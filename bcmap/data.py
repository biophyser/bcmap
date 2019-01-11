__description__ = \
"""
Constants and data used across the package.
"""
__author__ = "Michael J. Harms"
__date__ = "2019-01-12"

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
            '---':'-',
            'TCN':'S', 'CTN':'L', 'CCN':'P', 'CGN':'R',
            'ACN':'T', 'GTN':'V', 'CGN':'A', 'GGN':'G'}

# Dictionary mapping ascii value to Q-score
Q_DICT = {}
for i in range(33,76):
    Q_DICT[chr(i)] = 10**(-(i-33)/10)
