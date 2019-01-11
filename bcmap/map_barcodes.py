import numpy as np
import copy, itertools, os

BASE_DICT = {"A":[1.00,0.00,0.00,0.00],
             "C":[0.00,1.00,0.00,0.00],
             "G":[0.00,0.00,1.00,0.00],
             "T":[0.00,0.00,0.00,1.00],
             "U":[0.00,0.00,0.00,1.00],
             "W":[0.50,0.00,0.00,0.50],
             "S":[0.00,0.50,0.50,0.00],
             "M":[0.50,0.50,0.00,0.00],
             "K":[0.00,0.00,0.50,0.50],
             "R":[0.50,0.00,0.50,0.00],
             "Y":[0.00,0.50,0.00,0.50],
             "B":[0.00,0.33,0.33,0.33],
             "D":[0.33,0.00,0.33,0.33],
             "H":[0.33,0.33,0.00,0.33],
             "V":[0.33,0.33,0.33,0.00],
             "N":[0.25,0.25,0.25,0.25],
             "-":[0.00,0.00,0.00,0.00]}

# Convert BASE_DICT values to numpy arrays, forcing normalization
# to 1.0.
for k in BASE_DICT.keys():
    value = np.array(BASE_DICT[k])
    if np.sum(value) != 0:
        value = value/np.sum(value)
    BASE_DICT[k] = value

INDEX_TO_BASE = list("ACGT-")
BASE_TO_INDEX = dict([(k,i) for i, k in enumerate(INDEX_TO_BASE)])
BASE_TO_INDEX["N"] = -1
BASE_FREQS = BASE_DICT["N"]
LN_BASE_FREQS = np.log(BASE_DICT["N"])

# genetic code
GENCODE = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
            '---':'-'}

# Dictionary mapping ascii value to Q-score
Q_DICT = {}
for i in range(33,76):
    Q_DICT[chr(i)] = 10**(-(i-33)/10)

def translate(sequence):
    """
    Translate a nucleotide sequence into a protein sequence.  If there is a
    problem, write an "X" into the sequence.
    """

    try:
        return "".join([GENCODE[sequence[3*i:3*i+3]]
                        for i in range(len(sequence)//3)])
    except KeyError:
        out = []
        for i in range(len(sequence)//3):
            try:
                out.append(GENCODE[sequence[3*i:3*i+3]])
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

    def _load_reads(self,read_file):
        """
        Create a dictionary where the keys are clusters in cluster_dict
        and the values are tuples of sequence/phred scores taken from
        read_file.
        """

        key_size = len(list(self._cluster_dict.keys())[0])

        out_dict = {}
        with open(read_file) as f:
            for l in f:
                col = l.split("|")

                bc = col[0][:key_size].strip()
                bases = col[2].strip()
                quals = np.array([Q_DICT[k] for k in list(col[3].strip())])

                cluster = self._cluster_dict[bc]

                try:
                    out_dict[cluster].append((bases,quals))
                except KeyError:
                    out_dict[cluster] = [(bases,quals)]

        return out_dict

    def _call_consensus(self,n_probs,c_probs,cutoff=0.5):

        diag_mask = np.eye(4,dtype=np.bool)

        num_mismatch = 0
        consensus = []
        support = []
        for i in range(n_probs.shape[0]):

            n_data = n_probs[i,:]
            c_data = c_probs[i,:]

            # Signal from N channel
            if np.sum(n_data) > 0:

                # Signal from C channel, average signals
                if np.sum(c_data) > cutoff:
                    average_prob = (n_data + c_data)/2

                # No signal from C channel, just take N channel
                else:
                    average_prob = n_data

            # No signal from N channel
            else:

                # Signal from C channel, use that one
                if np.sum(c_data) > 0:
                    average_prob = c_data

                # No signal from either channel
                else:
                    average_prob = np.zeros(n_data.shape,dtype=np.float)

            above_cutoff_mask = average_prob >= cutoff

            # High-probability for one base
            if np.sum(above_cutoff_mask) == 1:
                base = np.argmax(average_prob)
                consensus.append(INDEX_TO_BASE[base])
                support.append(average_prob[base])

            # No good data
            elif np.sum(above_cutoff_mask) == 0:
                consensus.append("N")
                support.append(np.max(average_prob))

            # Conflicting data
            else:
                consensus.append("X")
                support.append(np.max(average_prob))
                num_mismatch += 1

        return consensus, support, num_mismatch

    def map_from_cluster_files(self,cluster_directory,key_size=25):

        cluster_files = os.listdir(cluster_directory)
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

            n_term_reads = []
            c_term_reads = []
            with open(os.path.join(cluster_directory,cf)) as f:
                for l in f:
                    col = l.split("|")

                    bc = col[0][:key_size].strip()
                    bases = col[2].strip()
                    quals = np.array([Q_DICT[k] for k in list(col[3].strip())])

                    if bases[-6:] == "TAATAA":
                        c_term_reads.append((bases,quals))
                    else:
                        n_term_reads.append((bases,quals))

                    # Write cluster/barcode pairs to a file
                    cluster_out.write("{}:{}\n".format(bc,cf.split("_")[0]))

            # Make sure that both N- and C-terminal sequences are seen
            if len(n_term_reads) == 0 or len(c_term_reads) == 0:
                continue

            # cluster number
            c = cf.split("_")[1].split(".")[0]

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
            consensus_out.write(">cluster_{},{}\n{}\n{}\n".format(c,
                                                                  num_mismatch,
                                                                  consensus,
                                                                  to_write))

            protein = translate(consensus)
            prot_out.write(">cluster_{},{}\n{}\n".format(c,num_mismatch,protein))


            # Flush output files
            if i % 1000 == 0:
                cluster_out.flush()
                consensus_out.flush()
                prot_out.flush()
                print(i,"of",len(clusters_seen))
                break

        # Close output files
        cluster_out.close()
        consensus_out.close()
        prot_out.close()




    def create_map(self,n_term_file,c_term_file,out_root):

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

        # Load in sequences, mapping them to unique cluster numbers
        n_term_reads = self._load_reads(n_term_file)
        c_term_reads = self._load_reads(c_term_file)

        # Create list of all clusters seen
        clusters_seen = list(n_term_reads.keys())
        clusters_seen.extend(c_term_reads.keys())
        clusters_seen = list(set(clusters_seen))
        clusters_seen.sort()

        # Open output file pipes
        cluster_out = open(cluster_out_file,"w")
        consensus_out = open(consensus_out_file,"w")
        prot_out = open(prot_out_file,"w")

        # Go through every cluster
        for i, c in enumerate(clusters_seen):

            # Make sure that N- and C-termini have both been seen in this
            # cluster
            try:
                n_term_reads[c]
                c_term_reads[c]
            except KeyError:
                continue

            # Get N- and C-terminal sequences for this cluster
            n_term_probs = self._ref.calc_base_prob(n_term_reads[c],align_three_prime=False)
            c_term_probs = self._ref.calc_base_prob(c_term_reads[c],align_three_prime=True)

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
            consensus_out.write(">cluster_{},{}\n{}\n{}\n".format(c,
                                                                  num_mismatch,
                                                                  consensus,
                                                                  to_write))

            protein = translate(consensus)
            prot_out.write(">cluster_{},{}\n{}\n".format(c,num_mismatch,protein))

            # Write cluster/barcode pairs to a file
            for bc in self._inverse_dict[c]:
                cluster_out.write("{}:{}\n".format(bc,c))

            # Flush output files
            if i % 1000 == 0:
                cluster_out.flush()
                consensus_out.flush()
                prot_out.flush()
                print(i,"of",len(clusters_seen))

        # Close output files
        cluster_out.close()
        consensus_out.close()
        prot_out.close()


class SequenceProfile:

    def __init__(self,sequence_file=None):

        if sequence_file is not None:
            self.read_sequence_file(sequence_file)

    def _seq_to_prob(self,input_array):
        """
        Convert a string representation of a DNA sequence to a probability
        array.

        ACGT becomes:

        [[1,0,0,0],
         [0,1,0,0],
         [0,0,1,0],
         [0,0,0,1]]
        """

        out = np.zeros((len(input_array),4),dtype=np.float)
        for i in range(len(input_array)):
            out[i,:] = BASE_DICT[input_array[i]]

        return out

    def _array_diff(self,array_one,array_two,remove_gaps=True):
        """
        Take RMSD difference between two arrays of the same shape.
        """

        if array_one.shape != array_two.shape:
            err = "arrays must have the same shape\n"
            raise ValueError(err)

        array1 = np.copy(array_one)
        array2 = np.copy(array_two)

        # Do not compare gaps, where probability density is 0
        # by definition.
        if remove_gaps:
            array2[np.sum(array1,1) == 0,:] = 0.0
            array1[np.sum(array2,1) == 0,:] = 0.0

        return np.sum((array1 - array2)**2)/array1.shape[0]

    def read_sequence_file(self,sequence_file,merge_gaps=True,max_gap_merge=3):
        """
        First line should be a sequence.
        Second line should be sequence with all possible mutations.
        merge_gaps: boolean.  whether or not to merge sequential gaps when
                    making all possible gap combinations (default to True)
        max_gap_merge: maximum number of sequential gaps to merge.  defaults to
                       3, assuming gaps are codons.
        """

        f = open(sequence_file,"r")
        lines = [l.strip() for l in f.readlines()
                 if not l.startswith("#") and len(l.strip()) != 0]
        f.close()

        failed = True
        if len(lines) == 2:
            if len(lines[0]) == len(lines[1]):
                failed = False

        if failed:
            err = "The file must have two non-blank lines of the same length.\n"
            err += "The first line is the wildtype sequence, the second line\n"
            err += "should have the sequence containing the possible mutations.\n"
            err += "The sequences should be gapped appropriately relative to one\n"
            err += "another."
            raise ValueError(err)

        wt = lines[0]
        mut = lines[1]

        self._total_length = len(wt)
        self._b

        # Walk though sequence, recording non-gap sequence (real_wt, real_mut)
        # and position of each gap (gaps)
        real_wt = []
        real_mut = []
        base_profile = []
        gaps = []
        for i in range(len(wt)):

            # Don't care about gaps shared with wildtype and mutant
            if wt[i] == "-" and mut[i] == "-":
                self._total_length -= 1
                continue

            base_profile.append([])

            # Record wildtype sequence
            if wt[i] != "-":
                real_wt.append(wt[i])
                base_profile[-1].append(wt[i])

            # Insertion
            else:
                assigned = False

                # Should we try to merge gaps?
                if merge_gaps:

                    # Is there a previous gap to look at?
                    if len(gaps) > 0:

                        # Is previous gap an insertion?
                        if gaps[-1][2] == 1:

                            # Is the current gap part of a contiguous block of gaps
                            # with the previous gap?
                            if (gaps[-1][0] + gaps[-1][1]) == i:

                                # Is the previous gap short enough to be added to?
                                if gaps[-1][1] < max_gap_merge:

                                    # Attach to previous gap
                                    gaps[-1][1] += 1
                                    assigned = True

                # New gap
                if not assigned:
                    gaps.append([i,1,1])

            # Record mutant sequence
            if mut[i] != "-":
                real_mut.append(mut[i])
                base_profile[-1].append(mut[i])

            # Deletion
            else:
                assigned = False

                # Should we try to merge gaps?
                if merge_gaps:

                    # Is there a previous gap to look at?
                    if len(gaps) > 0:

                        # Is previous gap a deletion?
                        if gaps[-1][2] == -1:

                            # Is the current gap part of a contiguous block of gaps
                            # with the previous gap?
                            if (gaps[-1][0] + gaps[-1][1]) == i:

                                # Is the previous gap short enough to be added to?
                                if gaps[-1][1] < max_gap_merge:

                                    # Attach to previous gap
                                    gaps[-1][1] += 1
                                    assigned = True

                # New gap
                if not assigned:
                    gaps.append([i,1,-1])

        # Probability of seeing a given base at a given position
        self._prob_profile = np.zeros((self._total_length,4),dtype=np.float)
        for i in range(len(base_profile)):
            p = np.zeros(4,dtype=np.float)
            for base in base_profile[i]:
                p += BASE_DICT[base]
            p = p/len(base_profile[i])
            self._prob_profile[i,:] = p

        # Convert wildtype and mutant sequences into probability arrays
        self._real_wt = self._seq_to_prob(real_wt)
        self._real_mut = self._seq_to_prob(real_mut)

        # Make all possible combinations of the gaps
        gap_combos = []
        for i in range(len(gaps)+1):
            for c in itertools.combinations(gaps,i):
                gap_combos.append(c)

        # Go throuh every possible gap combination
        self._gap_masks = []
        self._gap_indexes = []
        for combo in gap_combos:

            # Convert gap pattern into a boolean mask
            gap_mask = np.zeros(self._total_length,dtype=np.bool)
            for gap in combo:
                for i in range(gap[1]):
                    gap_mask[gap[0] + i] = True

            # Create indexes flowing around gap pattern
            gap_indexes = -1*np.ones(self._total_length,dtype=np.int)
            counter = 0
            for i in range(self._total_length):
                if not gap_mask[i]:
                    gap_indexes[i] = counter
                    counter += 1

            # Record gap patterns
            self._gap_masks.append(gap_mask)
            self._gap_indexes.append(gap_indexes)


    def _align_sequence(self,prob_array,align_three_prime=False):

        result = None

        for i in range(len(self._gap_masks)):

            gap_mask = np.copy(self._gap_masks[i])
            gap_indexes = np.copy(self._gap_indexes[i])

            # Copy in probability array into to_align
            to_align = np.copy(prob_array)

            if not align_three_prime:

                # Chop the global indexes down to the length of to_align
                # (which will generally be shorter than the total_length)
                end_slice = np.where(gap_indexes == len(to_align))[0]
                if len(end_slice) == 0:
                    end_slice = len(gap_indexes)
                else:
                    end_slice = end_slice[0]

                gap_indexes = gap_indexes[:end_slice]

                # Make gap mask as long as final output
                tmp_gap_mask = gap_mask[:end_slice]
                gap_mask = np.ones(self._total_length,dtype=np.bool)
                gap_mask[:len(tmp_gap_mask)] = tmp_gap_mask

            # Chopping the 5' end
            else:

                # Figure out how much to shift to line up right-most position
                # in to_align with right-most position in global_index
                to_shift = len(to_align) - np.max(gap_indexes) - 1
                gap_indexes = gap_indexes + to_shift
                gap_indexes[gap_mask] = -1

                # Look for 0 in the set of indexes.  Values to the left of
                # that will be less than zero and should be trimmed out
                start_slice = np.where(gap_indexes == 0)[0]

                if len(start_slice) == 0:
                    start_slice = 0
                else:
                    start_slice = start_slice[0]

                tmp_indexes = gap_indexes[start_slice:]
                gap_indexes = np.zeros(self._total_length,dtype=np.int)
                gap_indexes[(self._total_length - len(tmp_indexes)):] = tmp_indexes

                tmp_gap_mask = gap_mask[start_slice:]
                gap_mask = np.ones(self._total_length,dtype=np.bool)
                gap_mask[(self._total_length - len(tmp_gap_mask)):] = tmp_gap_mask

            # Indexes in final alignment
            target_indexes = np.array(range(len(gap_indexes)))

            # Now align the input array to this particular gap pattern
            alignment = np.zeros((self._total_length,4),dtype=np.float)
            alignment[:,:] = BASE_FREQS
            alignment[target_indexes,:] = to_align[gap_indexes,:]

            # Make all gaps 0.0 --> no probability there
            alignment[gap_mask,:] = 0.0

            # Calculate the difference between this alignment and the probability
            # profile from the input file.
            score = self._array_diff(self._prob_profile,alignment)
            num_gaps = np.sum(gap_mask)

            # Decide whether this alignment is better than the best one
            # we've seen
            keep_result = False
            if result is None:
                keep_result = True
            else:
                # If scores are very close, take one with fewer gaps
                if np.isclose(score,result[0]):
                    if num_gaps < result[1]:
                        keep_result = True
                else:

                    # If score is better than previously stored result,
                    # keep that
                    if score < result[0]:
                        keep_result = True

            # Record result
            if keep_result:
                result = [score,num_gaps,alignment]

        return result[2]

    def calc_base_prob(self,reads,align_three_prime=False):
        """
        Calculate the probability that each position in an alignment has
        a particular identity.  Returns a total_length x 4 array where each
        row sums to 1, indicating the probability that position is a
        particular base.  If a row sums to 0.0, it is a gap.
        """

        # somewhat cryptic call splits:
        # reads = [(s1,q1),(s2,q2),...(sn,qn)] into
        # (s1,s2,...,sn), (q1,q2,...,qn)
        raw_bases, raw_quals = zip(*reads)

        # Load bases from reads and quality scores into stacks
        base_stack = -1*np.ones((len(raw_bases),self._total_length),dtype=np.int)
        qual_stack = np.zeros((len(raw_bases),self._total_length),dtype=np.float)
        for i in range(len(raw_bases)):
            N = len(raw_bases[i])
            for j in range(N):
                if not align_three_prime:
                    base_stack[i,j] = BASE_TO_INDEX[raw_bases[i][j]]
                    qual_stack[i,j] = raw_quals[i][j]
                else:
                    offset = self._total_length - N
                    base_stack[i,(offset + j)] = BASE_TO_INDEX[raw_bases[i][j]]
                    qual_stack[i,(offset + j)] = raw_quals[i][j]


        out = np.zeros((self._total_length,4),dtype=np.float)
        for i in range(self._total_length):

            has_data = base_stack[:,i] != -1

            # No data for this column; dump raw frequencies to indicate
            # uncertainty at this column.
            if np.sum(has_data) == 0:
                out[i,:] = BASE_FREQS
                continue

            last_column_with_data = i

            # Non-gap bases and quals
            bases = base_stack[:,i][has_data]
            quals = qual_stack[:,i][has_data]

            # Probability the call is correct given PHRED score
            ln_not_e = np.log(1 - quals)

            # Probability the call is wrong given PHRED score.
            # There are "degeneracy" ways to be wrong. The probability of being
            # wrong is split equally among those possibilities.  (This assumes the
            # base frequencies are equal... close enough for our purposes)
            ln_e = np.log(quals)/3

            # Calculate ln[P(base|vector_of_observed_bases)]
            for j in range(4):
                out[i,j] += np.sum(ln_not_e[bases == j]) + np.sum(ln_e[bases != j])

            # Normalize
            e_out = np.exp(out[i,:])
            out[i,:] = e_out/np.sum(e_out)

        # Only take data for columns that were actually seen.
        out = out[:(last_column_with_data + 1),:]
        aligned = self._align_sequence(out,align_three_prime)

        final_out = 0.20*np.ones((self._total_length,5),dtype=np.float)
        final_out[:,:4] = aligned


        return out


    @property
    def total_length(self):
        return self._total_length
