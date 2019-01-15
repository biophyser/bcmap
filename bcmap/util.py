__description__ = \
"""
Utility functions used throughout pacakge.
"""
__author__ = "Michael J. Harms"
__date__ = "2019-01-12"

from . import data
import numpy as np

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


def prob_to_seq(input_array,cutoff=0.8):

    out = []
    for i in range(input_array.shape[0]):

        if np.max(input_array[i,:])  > cutoff:
            base = np.argmax(input_array[i,:])
            base = data.INDEX_TO_BASE[base]
        else:
            base = "N"

        out.append(base)

    return "".join(out)



def seq_to_prob(sequence):
    """
    Convert a string representation of a DNA sequence to a probability
    array.

    For example, ACGT becomes:

    [[1,0,0,0],
     [0,1,0,0],
     [0,0,1,0],
     [0,0,0,1]]
    """

    out = np.zeros((len(sequence),4),dtype=np.float)
    for i in range(len(sequence)):
        out[i,:] = data.BASE_DICT[sequence[i]]

    return out

def array_diff(array_one,array_two,remove_gaps=True):
    """
    Take RMSD difference between two arrays of the same shape.
    """

    if array_one.shape != array_two.shape:
        err = "arrays must have the same shape\n"
        raise ValueError(err)

    # Do not compare gaps, where probability density is 0
    # by definition.
    if remove_gaps:
        array1 = np.copy(array_one)
        array2 = np.copy(array_two)
        array2[np.sum(array1,1) == 0,:] = 0.0
        array1[np.sum(array2,1) == 0,:] = 0.0
    else:
        array1 = array_one
        array2 = array_two

    return np.sum((array1 - array2)**2)/array1.shape[0]

def codon_degeneracy(codon_dict):
    """
    Create a dictionary mapping codons with a single ambiguous position to
    all possible amino acids given the bases known.
    """

    degenerate_codon_dict = {}

    codons = list(codon_dict.keys())
    for c in codons:

        # The original codon should be in here
        degenerate_codon_dict[c] = [codon_dict[c]]

        # Do not treat previous degenerate codons and gapped codons
        if len([pos for pos in c if pos in ["N","-"]]) > 0:
            continue

        # Go over all three possible positions
        for i in range(3):

            # Replace the position with "N"
            new_codon = list(c)
            new_codon[i] = "N"
            new_key = "".join(new_codon)

            # Make a new entry in the dengenerate_codon dictionary
            try:
                degenerate_codon_dict[new_key]
            except KeyError:
                degenerate_codon_dict[new_key] = []

            # Replace "N" with other bases and add possibilities to
            # degenerate_codon_dict
            tmp_codon = new_codon[:]
            for new_base in "ACGT":
                tmp_codon[i] = new_base
                tmp_key = "".join(tmp_codon)
                amino_acid = codon_dict[tmp_key]

                # Whatever amino acid is encoded by the temporary codon is
                # a possibility for a codon with "N" at this position.
                degenerate_codon_dict[new_key].append(amino_acid)

    # Get rid of duplicated amino acids for each degenerate codon
    for k in degenerate_codon_dict.keys():
        degenerate_codon_dict[k] = set(degenerate_codon_dict[k])

    return degenerate_codon_dict
