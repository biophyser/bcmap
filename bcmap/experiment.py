
from . import data
from . import util
from .profile import SequenceProfile

import numpy as np
import os, sys


class Experiment:
    """
    """

    def __init__(self,ref_file,key_size=25):

        self._ref_file = ref_file
        self._key_size = key_size
        self._ref = SequenceProfile(self._ref_file)

    def _call_consensus(self,n_probs,c_probs,signal_cutoff=0.3,call_cutoff=0.80):
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

    def _run_batch(self,batch):
        """
        Call the genotypes for each of the clusters in batch.
        """

        cons_buffer = []
        prot_buffer = []
        clust_buffer = []

        for b in batch:
            cons, prot, clust = self._call_genotype(b)
            if cons is not None:
                cons_buffer.append(cons)
                prot_buffer.append(prot)
                clust_buffer.append(clust)

        return cons_buffer, prot_buffer, clust_buffer


    def _call_genotype(self,cf):
        """
        Call the genotype of a given cluster. Store the outputs in the
        appropriate buffers.
        """

        # cluster number
        cluster_number = cf.split("_")[1].split(".")[0]

        n_term_reads = []
        c_term_reads = []
        unique_bc = {}
        with open(os.path.join(self._cluster_directory,cf.strip())) as f:
            for l in f:
                col = l.split("|")

                bc = col[0][:self._key_size].strip()
                bases = col[2].strip()
                quals = np.array([data.Q_DICT[k] for k in list(col[3].strip())])

                # HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
                # Hacking to deal with issues with N- and C-terminal
                # read assignment and truncation.  For production, this
                # should go into process_reads.
                if bases[-9:] == "GCCTAATAA":
                    c_term_reads.append((bases,quals))
                else:
                    if len(bases) > 345:
                        bases = bases[:345]
                        quals = quals[:345]

                    # Trim off any extra TAA
                    while bases[-3:] == "TAA":
                        bases = bases[:-3]
                        quals = quals[:-3]

                    n_term_reads.append((bases,quals))
                # HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK

                unique_bc[bc] = cluster_number

        num_n = len(n_term_reads)
        num_c = len(c_term_reads)

        # Make sure that both N- and C-terminal sequences are seen
        if num_n == 0 or num_c == 0:
            return None, None, None

        # somewhat cryptic call splits:
        # reads = [(s1,q1),(s2,q2),...(sn,qn)] into
        # (s1,s2,...,sn), (q1,q2,...,qn)
        n_bases, n_quals = zip(*n_term_reads)
        c_bases, c_quals = zip(*c_term_reads)

        # Get N- and C-terminal sequences for this cluster
        n_term_probs = self._ref.calc_base_prob(n_bases,n_quals,align_three_prime=False)
        c_term_probs = self._ref.calc_base_prob(c_bases,c_quals,align_three_prime=True)

        # Get the consensus sequence defined by these
        consensus, support, num_mismatch = self._call_consensus(n_term_probs,
                                                                c_term_probs)
        consensus = "".join(consensus)

        #Create human-readable support
        support_to_write = []
        for s in support:
            v = "{:.1f}".format(np.round(s))
            if v == "1.0": v = "0.9"
            support_to_write.append(v[-1])
        support_to_write = "".join(support_to_write)

        protein = self._ref.translate(consensus)
        protein = "".join(protein)

        allele = self._ref.get_allele(consensus)
        allele = "".join(allele)

        # Write outputs
        cons = []
        cons.append(">cluster_{},{},{},{},{}\n".format(cluster_number,
                                                       allele,
                                                       num_n,
                                                       num_c,
                                                       num_mismatch))
        cons.append("{}\n".format(consensus))
        cons.append("{}\n".format(support_to_write))
        cons.append("{}\n".format(util.prob_to_seq(n_term_probs)))
        cons.append("{}\n".format(util.prob_to_seq(c_term_probs)))

        prot = []
        prot.append(">cluster_{},{},{},{},{}\n".format(cluster_number,
                                                       allele,
                                                       num_n,
                                                       num_c,
                                                       num_mismatch))
        prot.append("{}\n".format(protein))

        # Write out cluster numbers
        clust = []
        for bc in unique_bc.keys():
            # Write cluster/barcode pairs to a file
            clust.append("{},{},{}\n".format(bc,cluster_number,allele))

        return "".join(cons), "".join(prot), "".join(clust)


    def map_from_cluster_files(self,list_of_cluster_files,cluster_directory,
                               out_root,key_size=25,batch_size=100):

        self._processes = []

        self._list_of_cluster_files = list_of_cluster_files
        self._cluster_directory = cluster_directory
        self._out_root = out_root
        self._key_size = key_size

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

        self._cluster_buffer = []
        self._consensus_buffer = []
        self._prot_buffer = []

        i = 0
        batch_queue = []
        with open(list_of_cluster_files) as list_pointer:
            for cf in list_pointer:

                if len(batch_queue) < batch_size:
                    batch_queue.append(cf)
                else:

                    # Dump this call onto a thread
                    cons, prot, clust = self._run_batch(batch_queue)

                    # Do this on main thread
                    consensus_out.write("".join(cons))
                    prot_out.write("".join(prot))
                    cluster_out.write("".join(clust))

                    batch_queue = []

                # Flush output files
                if i % batch_size == 0:

                    cluster_out.flush()
                    consensus_out.flush()
                    prot_out.flush()
                    print(i)
                    sys.stdout.flush()

                i += 1

        # Clear any remaining in the batch queue
        cons, prot, clust = self._run_batch(batch_queue)
        consensus_out.write("".join(cons))
        prot_out.write("".join(prot))
        cluster_out.write("".join(clust))

        # Close output files
        cluster_out.close()
        consensus_out.close()
        prot_out.close()
