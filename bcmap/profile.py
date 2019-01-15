
from . import data
from . import util

import numpy as np
import itertools, copy

class SequenceProfile:

    def __init__(self,ref_file,merge_gaps=True,max_gap_merge=3):
        """
        ref_file: File describing the expected characteristics of the library.
                  The first line should be a sequence, the second line should be
                  a sequence with all possible mutations to the first.  The
                  two sequences should be gapped appropriately relative to each
                  other.  Gaps in the first sequence indicate insertions in the
                  mutant.  Gaps in the second sequence indicate deletions in the
                  mutant.  These sequences should use the standard nucleic acid
                  notation (ACGT-).  Degenerate positions are also allowed
                  (e.g.: N means any base, R means A or G, etc.)
        merge_gaps: boolean.  whether or not to merge sequential gaps when
                    making all possible gap combinations (default to True).
        max_gap_merge: Maximum number of sequential gaps to merge.  Defaults to
                       3, assuming gaps are codons.  This means six sequential
                       gaps would be treated as two blocks of three gaps.
        """

        self._ref_file = ref_file
        self._merge_gaps = merge_gaps
        self._max_gap_merge = max_gap_merge
        self._read_ref_file()


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

        # Filter out positions where both wt and mut are gaps
        tmp_wt = []
        tmp_mut = []
        for i in range(len(wt)):
            if wt[i] == "-" and mut[i] == "-":
                continue
            tmp_wt.append(wt[i])
            tmp_mut.append(mut[i])

        self._wt = tmp_wt[:]
        self._mut = tmp_mut[:]
        self._total_length = len(self._wt)

        # ------------------------------------------------------------------- #
        # Calculate total frequencies of bases seen in the alignment
        # ------------------------------------------------------------------- #

        all_sequence_data = list(self._wt)
        all_sequence_data.extend(list(self._mut))
        num_gaps = sum([1 for s in all_sequence_data if s == "-"])
        non_gaps = sum([1 for s in all_sequence_data if s != "-"])
        total = num_gaps + non_gaps

        freq = np.zeros(4,dtype=np.float)
        for s in all_sequence_data:
            if s == "-":
                continue
            else:
                freq[:] += data.BASE_DICT[s]

        self._base_freq_no_gap = freq/non_gaps
        self._ln_base_freq_no_gap = np.log(self._base_freq_no_gap)

        freq_with_gap = np.zeros(5,dtype=np.float)
        freq_with_gap[:4] = self._base_freq_no_gap*(non_gaps/total)
        freq_with_gap[4] = num_gaps/total

        self._base_freq_with_gap = freq_with_gap

        # ------------------------------------------------------------------- #
        # Parse gaps in the alignment
        # ------------------------------------------------------------------- #

        # Walk though sequence, recording non-gap sequence (wt_no_gap, mut_no_gap)
        # and position of each gap (gaps)
        wt_no_gap = []
        mut_no_gap = []
        base_profile = []
        gaps = []
        for i in range(len(self._wt)):

            base_profile.append([])

            # Record wildtype sequence
            if self._wt[i] != "-":
                wt_no_gap.append(self._wt[i])
                base_profile[-1].append(self._wt[i])

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
            if self._mut[i] != "-":
                mut_no_gap.append(self._mut[i])
                base_profile[-1].append(self._mut[i])

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
        self._wt_no_gap = util.seq_to_prob(wt_no_gap)
        self._mut_no_gap = util.seq_to_prob(mut_no_gap)

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
            score = util.array_diff(self._prob_profile,alignment)
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

        alignment = result[2]

        # Finally, everything after the last observation should be ambiguous,
        # not a gap.
        gaps = np.sum(alignment,1) == 0
        indexes = np.arange(alignment.shape[0],dtype=np.int)
        if not align_three_prime:
            # Set gaps to -1.  The biggest value in the index array will now be
            # the rightmost non-gap character.  Set everything after this to
            # the raw base frequencies
            indexes[gaps] = -1
            last_non_gap = np.max(indexes) + 1
            alignment[last_non_gap:,:] = self._base_freq_no_gap
        else:
            # Set gaps to array maximum.  The smallest value in index_array will
            # now be the leftmost non-gap character.  Set everything before this
            # to the raw base frequencies.
            indexes[gaps] = np.max(indexes)
            first_non_gap = np.argmin(indexes)
            alignment[:first_non_gap,:] = self._base_freq_no_gap

        return alignment


    def calc_base_prob(self,raw_bases,raw_quals,align_three_prime=False):
        """
        Calculate the probability that each position in an alignment has
        a particular identity.  Returns a total_length x 5 array where each
        row sums to 1, indicating the probability that position is a
        particular base.
        """

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

            # ln[P(base)] ...
            out[i,:] = self._ln_base_freq_no_gap

            # ... + ln[P(vector_of_observed_bases|base)]
            for j in range(4):
                out[i,j] += np.sum(ln_not_e[bases == j]) + np.sum(ln_e[bases != j])

            # Set best value to 0. This avoids numerical problems when taking
            # exponential.  The best is always 0 -> exp(0) = 1.  The worst could
            # be exp(-10000) -> 0, but safely b/c denominator is always at least
            # one
            out[i,:] = out[i,:] - np.max(out[i,:])

            # Normalize (sum(e_out) = P(vector_of_observed_bases))
            e_out = np.exp(out[i,:])
            out[i,:] = e_out/np.sum(e_out)

        # Only take data for columns that were actually seen.
        out = out[:(last_column_with_data + 1),:]

        # Find best alignment
        aligned = self._align_sequence(out,align_three_prime)

        # Construct final output -- if we do not have data, set to bases at
        # frequencies observed in reference sequence
        final_out = np.zeros((self._total_length,5),dtype=np.float)

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

    def translate(self,base_sequence):
        """
        Translate a base sequence, using the reference sequence to call codons
        with a missing base or two.
        """

        if len(base_sequence) != self._total_length:
            err = "input sequence length must be same as sequence profile\n"
            raise ValueError(err)

        translation = []
        for i in range(len(self.amino_acid_profile)):
            try:
                aa = self.amino_acid_profile[i][base_sequence[i*3:((i+1)*3)]]
            except KeyError:
                aa = "X"
            translation.append(aa)

        return translation

    def get_allele(self,base_sequence):
        """
        Get the allele given by a base sequence.  The wildtype amino acid is 0,
        other amino acids are 1, 2, 3 ... for other possibilities.
        """

        if len(base_sequence) != self._total_length:
            err = "input sequence length must be same as sequence profile\n"
            raise ValueError(err)

        allele = []
        translation = self.translate(base_sequence)
        for i, aa in enumerate(translation):
            if self.allele_call[i] is not None:

                try:
                    allele.append(self.allele_call[i][aa])
                except KeyError:
                    allele.append("?")

        return allele

    def allele_summary(self):
        """
        Summarize the mapping between amino acids and their summarized
        allele calls.
        """

        for i in range(len(self.allele_call)):

            if self.allele_call[i] is not None:

                to_write = []
                for k in self.allele_call[i].keys():
                    to_write.append((self.allele_call[i][k],k))
                to_write.sort()

                out = ["{}: ".format(i+1)]
                for j in range(len(to_write)):
                    out.append("{} = {}; ".format(to_write[j][0],to_write[j][1]))

                print("".join(out))


    @property
    def allele_call(self):
        """
        Dictionary returning the allele call for a given amino acid at a given
        position.
        """

        try:
            return self._allele_call
        except AttributeError:
            pass

        possible_alleles = "0123456789ABCDEFGHIJK"

        allele_call = []
        for i in range(len(self.amino_acid_profile)):
            aa = list(set(self.amino_acid_profile[i].values()))

            # Interesting call
            if len(aa) == 1:
                allele_call.append(None)
            else:
                call_dict = {}
                wt_aa = data.GENCODE["".join(self._wt[(i*3):((i+1)*3)]) ]
                call_dict[wt_aa] = "0"
                aa.remove(wt_aa)
                aa.sort()
                for j, allele in enumerate(aa):
                    call_dict[allele] = possible_alleles[j+1]
                allele_call.append(call_dict)

        self._allele_call = copy.deepcopy(allele_call)

        return self._allele_call

    @property
    def library_combos(self):
        """
        Possible amino acid combinations in the library.
        """

        try:
            return self._library_combos
        except AttributeError:
            pass

        to_expand = []
        for i in range(len(self.amino_acid_profile)):
            aa = set(self.amino_acid_profile[i].values())
            to_expand.append(aa)

        self._library_combos = list(itertools.product(*to_expand))

        return self._library_combos

    @property
    def amino_acid_profile(self):
        """
        Position-specific codon lookup table.
        """

        # Only calculate if you need it
        try:
            return self._amino_acid_profile
        except AttributeError:
            pass

        self._degenerate_table = util.codon_degeneracy(data.GENCODE)

        if self._total_length % 3 != 0:
            err = "sequence must be divisible by three (in frame)\n"
            raise ValueError(err)

        self._amino_acid_profile = []
        for i in range(0,self._total_length,3):

            possible_codons = []
            possible_aa = []

            # Set of three base positions in both wildtype and mutant
            triplet_bases = [self._wt[i:i+3]]
            if self._mut[i:i+3] != self._wt[i:i+3]:
                triplet_bases.append(self._mut[i:i+3])

            # For each of the base triplets
            for j in range(len(triplet_bases)):

                # Unpack into all possible bases at each position.  (This
                # deals with notation like R -> A or T).
                triplet = triplet_bases[j]
                possible_bases = []
                for k in range(3):
                    possible_bases.append(data.UNPACK_BASE[triplet[k]])

                # Construct a list of all possible codons given these triplets
                for c in itertools.product(*possible_bases):
                    possible_codons.append("".join(c))
                    possible_aa.append(data.GENCODE[possible_codons[-1]])

            # key a list of codons to each possible amino acid.  Include
            # codons with single "N" values in them.
            aa_dict = {}
            for j in range(len(possible_aa)):

                aa = possible_aa[j]

                try:
                    aa_dict[aa].append(possible_codons[j])
                except KeyError:
                    aa_dict[aa] = [possible_codons[j]]

                to_expand = [(possible_codons[j][0],"N"),
                             (possible_codons[j][1],"N"),
                             (possible_codons[j][2],"N")]

                for c in itertools.product(*to_expand):
                    aa_dict[aa].append("".join(c))

            # Convert lists of codons to sets
            all_aa_seen = list(aa_dict.keys())
            for aa in all_aa_seen:
                aa_dict[aa] = set(aa_dict[aa])

            # Go through all pairwise combinations of codons seen in for
            # different amino acids to find any codons (some with "N" at a
            # single position that overlap at this site)
            ambiguous = []
            for j in range(len(all_aa_seen)):
                set_j = aa_dict[all_aa_seen[j]]
                for k in range(j+1,len(all_aa_seen)):
                    set_k = aa_dict[all_aa_seen[k]]
                    overlap = set_j.intersection(set_k)
                    ambiguous.extend(overlap)

            final_codon_dict = {}
            for aa in all_aa_seen:
                for codon in aa_dict[aa]:
                    if codon not in ambiguous:
                        final_codon_dict[codon] = aa

            self._amino_acid_profile.append(copy.deepcopy(final_codon_dict))

        return self._amino_acid_profile

    @property
    def total_length(self):
        return self._total_length

    @property
    def base_freq_no_gap(self):
        return self._base_freq_no_gap

    @property
    def base_freq_with_gap(self):
        return self._base_freq_with_gap
