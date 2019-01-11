
from . import data

import numpy as np
import itertools

class SequenceProfile:

    def __init__(self,ref_file,merge_gaps=True,max_gap_merge=3):
        """
        ref_file: file describing the expected characteristics of the library
                  First line should be a sequence.
                  Second line should be sequence with all possible mutations.
        merge_gaps: boolean.  whether or not to merge sequential gaps when
                    making all possible gap combinations (default to True)
        max_gap_merge: maximum number of sequential gaps to merge.  defaults to
                       3, assuming gaps are codons.
        """"

        self._ref_file = ref_file
        self._merge_gaps = merge_gaps
        self._max_gap_merge = max_gap_merge
        self._read_ref_file()

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
            out[i,:] = data.BASE_DICT[input_array[i]]

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

    def _read_ref_file(self):
        """
        Read a reference input file.
        """

        f = open(self._ref_file,"r")
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

        # ------------------------------------------------------------------- #
        # Calcualate frequencies of bases seen in alignmnet
        # ------------------------------------------------------------------- #

        all_sequence_data = list(wt)
        all_sequence_data.extend(list(mut))
        num_gaps = sum([1 for i s in all_sequence_data if s == "-"])
        non_gaps = sum([1 for i s in all_sequence_data if s != "-"])
        total = num_gaps + non_gaps

        freq = np.zeros(4,dtype=np.float)
        for s in all_sequence_data:
            if s == "-":
                continue
            else:
                freq[:] += data.BASE_DICT[s]

        self._base_freq_no_gap = freq/non_gaps

        freq_with_gap = np.zeros(5,dtype=np.float)
        freq_with_gap[:4] = self._base_freq_no_gap*(non_gaps/total)
        freq_with_gap[4] = num_gaps/total

        self._base_freq_with_gap = freq_with_gap

        # ------------------------------------------------------------------- #
        # Parse gaps in the alignment
        # ------------------------------------------------------------------- #

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
                if self._merge_gaps:

                    # Is there a previous gap to look at?
                    if len(gaps) > 0:

                        # Is previous gap an insertion?
                        if gaps[-1][2] == 1:

                            # Is the current gap part of a contiguous block of gaps
                            # with the previous gap?
                            if (gaps[-1][0] + gaps[-1][1]) == i:

                                # Is the previous gap short enough to be added to?
                                if gaps[-1][1] < self._max_gap_merge:

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
                if self._merge_gaps:

                    # Is there a previous gap to look at?
                    if len(gaps) > 0:

                        # Is previous gap a deletion?
                        if gaps[-1][2] == -1:

                            # Is the current gap part of a contiguous block of gaps
                            # with the previous gap?
                            if (gaps[-1][0] + gaps[-1][1]) == i:

                                # Is the previous gap short enough to be added to?
                                if gaps[-1][1] < self._max_gap_merge:

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
                p += data.BASE_DICT[base]
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
        """
        Align the sequence stored in probaility array to the reference
        sequence data.  If align_three_prime, treat as a reverse read.
        """

        result = None

        # Try every gap combination
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
            alignment[:,:] = self._base_freq_no_gap
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
                    base_stack[i,j] = data.BASE_TO_INDEX[raw_bases[i][j]]
                    qual_stack[i,j] = raw_quals[i][j]
                else:
                    offset = self._total_length - N
                    base_stack[i,(offset + j)] = data.BASE_TO_INDEX[raw_bases[i][j]]
                    qual_stack[i,(offset + j)] = raw_quals[i][j]


        out = np.zeros((self._total_length,4),dtype=np.float)
        for i in range(self._total_length):

            has_data = base_stack[:,i] != -1

            # No data for this column; dump raw frequencies to indicate
            # uncertainty at this column.
            if np.sum(has_data) == 0:
                out[i,:] = self._base_freq_no_gap
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

        # Find best alignment
        aligned = self._align_sequence(out,align_three_prime)

        # Construct final output -- if we do not have data, set to bases at
        # frequencies observed in reference sequence
        final_out = np.zeros((self._total_length,5),dtype=np.float)
        final_out[:,:] = self._base_freq_with_gap

        # Load in aligned data
        final_out[:,:4] = aligned

        # Expand totally ambiguous positions to include the possibility of being
        # a gap
        ambig = np.sum(aligned[:] == self._base_freq_no_gap,1) == aligned.shape[1]
        final_out[ambig,:] = self._base_freq_with_gap

        # Set gaps
        gap_positions = np.sum(aligned,1) == 0
        final_out[gap_positions,:] = 0.0
        final_out[gap_positions,4] = 1.0

        return final_out

    @property
    def total_length(self):
        return self._total_length

    @property
    def base_freq_no_gap(self):
        return self._base_freq_no_gap

    @property
    def base_freq_with_gap(self):
        return self._base_freq_with_gap
