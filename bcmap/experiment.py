
from . import data
from .profile import SequenceProfile

import numpy as np
import os

def translate(sequence):
    """
    Translate a nucleotide sequence into a protein sequence.  If there is a
    problem, write an "X" into the sequence.
    """

    try:
        return "".join([data.GENCODE[sequence[3*i:3*i+3]]
                        for i in range(len(sequence)//3)])
    except KeyError:
        out = []
        for i in range(len(sequence)//3):
            try:
                out.append(data.GENCODE[sequence[3*i:3*i+3]])
            except KeyError:
                out.append("X")
        return "".join(out)

class Experiment:
    """
    """

    def __init__(self,ref_file,cluster_file,key_size=25):

        self._ref_file = ref_file
        self._cluster_file = cluster_file
        self._key_size = key_size

        self._ref = SequenceProfile(self._ref_file)
        self._load_clusters()

    def _load_clusters(self):
        """
        Load a file containing barcodes placed into clusters.

        Returns:
            a dictionary with barcodes as keys and cluster indexes as values.
            a dictionary with clusters as keys and lists of barcodes as values.
        """

        cluster_counter = -1
        cluster_dict = {}
        inverse_dict = {}
        with open(self._cluster_file) as f:
            for l in f:
                col = l.split()
                cluster = int(col[0])
                if cluster == -1:
                    cluster = cluster_counter
                    cluster_counter = cluster_counter - 1
                cluster_dict[col[1].strip()[:self._key_size]] = cluster

                try:
                    inverse_dict[cluster].append(col[1].strip()[:self._key_size])
                except KeyError:
                    inverse_dict[cluster] = [col[1].strip()[:self._key_size]]

        self._cluster_dict = cluster_dict
        self._inverse_dict = inverse_dict

    def _call_consensus(self,n_probs,c_probs,signal_cutoff=0.3,call_cutoff=0.6):
        """
        Call the consensus between N-terminal and C-terminal reads.
        """

        num_mismatch = 0
        consensus = []
        support = []
        for i in range(n_probs.shape[0]):

            n_data = n_probs[i,:]
            c_data = c_probs[i,:]

            # Do n and c signals differ from base_freq_with_gap (i.e. no
            # signal at all)?
            n_has_signal = np.max(n_data) > signal_cutoff
            c_has_signal = np.max(c_data) > signal_cutoff

            # Average signals depending on whether they have signal
            if n_has_signal and c_has_signal:
                average_prob = (n_data + c_data)/2
            elif n_has_signal and not c_has_signal:
                average_prob = n_data
            elif not n_has_signal and c_has_signal:
                average_prob = c_data
            else:
                average_prob = (n_data + c_data)/2

            # See how many bases are above call_cutoff in their
            # probability
            above_cutoff_mask = average_prob >= call_cutoff

            # High-probability for one base
            if np.sum(above_cutoff_mask) == 1:
                base = np.argmax(average_prob)
                consensus.append(data.INDEX_TO_BASE[base])
                support.append(average_prob[base])

            # No good data
            elif np.sum(above_cutoff_mask) == 0:

                # If the value is low because of two, conflicting, votes
                # record as a mismatch.
                if np.max(n_data) > call_cutoff and np.max(c_data) > call_cutoff:
                    consensus.append("X")
                    support.append(np.max(average_prob))
                    num_mismatch += 1

                # Otherwise, it is just poor data.
                else:
                    consensus.append("N")
                    support.append(np.max(average_prob))

            # Conflicting data (won't ever happen if call_cutoff is > 0.5)
            else:
                consensus.append("X")
                support.append(np.max(average_prob))
                num_mismatch += 1

        return consensus, support, num_mismatch

    def map_from_cluster_files(self,cluster_directory,key_size=25):

        cluster_files = os.listdir(cluster_directory)
        cluster_files.sort()
        out_root = cluster_directory

        # Make output file names
        cluster_out_file = "{}.clusters".format(out_root)
        consensus_out_file = "{}_dna.fasta".format(out_root)
        prot_out_file = "{}_prot.fasta".format(out_root)

        # Make sure output files do not exist
        if os.path.isfile(cluster_out_file) or \
           os.path.isfile(consensus_out_file) or \
           os.path.isfile(prot_out_file):
           err = "output file(s) exist\n"
           raise FileExistsError(err)

        # Open output file pipes
        cluster_out = open(cluster_out_file,"w")
        consensus_out = open(consensus_out_file,"w")
        prot_out = open(prot_out_file,"w")

        for i, cf in enumerate(cluster_files):

            # cluster number
            cluster_number = cf.split("_")[1].split(".")[0]

            n_term_reads = []
            c_term_reads = []
            unique_bc = {}
            with open(os.path.join(cluster_directory,cf)) as f:
                for l in f:
                    col = l.split("|")

                    bc = col[0][:key_size].strip()
                    bases = col[2].strip()
                    quals = np.array([data.Q_DICT[k] for k in list(col[3].strip())])

                    ##### SOME SERIOUS HACKING HERE --> THIS SHOULD GO UP
                    # IN PROCESS READS FOR PRODUCTION XXXXX
                    if bases[-9:] == "GCCTAATAA":
                        c_term_reads.append((bases,quals))
                    else:
                        ## GET RID OF TAATAA FRAGMENT...
                        tmp_bases = bases[:-5]
                        tmp_quals = quals[:-5]
                        n_term_reads.append((tmp_bases,tmp_quals))
                    ##HACK HACK HACK

                    unique_bc[bc] = cluster_number

            # Make sure that both N- and C-terminal sequences are seen
            if len(n_term_reads) == 0 or len(c_term_reads) == 0:
                continue

            # Get N- and C-terminal sequences for this cluster
            n_term_probs = self._ref.calc_base_prob(n_term_reads,align_three_prime=False)
            c_term_probs = self._ref.calc_base_prob(c_term_reads,align_three_prime=True)

            # Get the consensus sequence defined by these
            consensus, support, num_mismatch = self._call_consensus(n_term_probs,
                                                                    c_term_probs)

            #Create human-readable support
            to_write = []
            for s in support:
                v = "{:.1f}".format(np.round(s))
                if v == "1.0": v = "0.9"

                to_write.append(v[-1])
            to_write = "".join(to_write)

            consensus = "".join(consensus)
            consensus_out.write(">cluster_{},{}\n{}\n{}\n".format(cluster_number,
                                                                  num_mismatch,
                                                                  consensus,
                                                                  to_write))

            protein = translate(consensus)
            prot_out.write(">cluster_{},{}\n{}\n".format(cluster_number,
                                                         num_mismatch,
                                                         protein))

            # Write out cluster numbers
            for bc in unique_bc.keys():
                # Write cluster/barcode pairs to a file
                cluster_out.write("{}:{}\n".format(bc,cluster_number))

            # Flush output files
            if i % 1000 == 0:
                cluster_out.flush()
                consensus_out.flush()
                prot_out.flush()
                print(i,"of",len(cluster_files))

        # Close output files
        cluster_out.close()
        consensus_out.close()
        prot_out.close()
