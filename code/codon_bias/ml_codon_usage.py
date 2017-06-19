from collections import Counter
from itertools import islice
import numpy as np
import pandas as pd
from pyfaidx import Fasta
# import random
import re
from sklearn.svm import SVR, SVC
from sklearn.feature_selection import RFE, RFECV
from sklearn.cross_validation import StratifiedKFold
import sys


sys.path.insert(0, '.../influenza')
from count import actg_dict, count_kmer


# Plan.
'''
We already filtered some spliced proteins (those annotated with
brackets). Still our fasta contains duplicate accession numbers.
Which proteins are annotated in the NCBI flu db?
grep '>' FASTA.fa | awk -F'[_]' '{print $3}' - | sort | uniq -c
>gb|CY110306:1-2274_2_PB1
>gb|CY110306:95-367_2_PB1-F2  # this is the problematic case we filter

Then:
1. Sample interesting, clean subset (training) from lookup table.
2. Query fasta with pyfaidx.
3. Process and aggregate query, e.g. count unigrams (codons) etc.
4. ML.
5. Test on test set
'''


# Import fasta file and index it for quick random access.
file = 'iva_accession.fa'
fasta = Fasta(file, as_raw=True)  # we only want the sequence
# fasta['CAB95856'][:] or str(fasta['CAB95856'])


# Lookup table.
lu = pd.read_csv(
    'genomeset.map', sep='\t', index_col=None,
    keep_default_na=False, na_values=['_'])
proteinset = {'PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M1', 'NS1'}
lu = lu.loc[lu['protein_name'].isin(proteinset)]
# Counter(sub.protein_name)
# Counter({'HA': 23880,
#          'M1': 23880,
#          'NA': 23878,
#          'NP': 23879,
#          'NS1': 23879,
#          'PA': 23869,
#          'PB1': 23864,
#          'PB2': 23867})
# Clean enough.


# Sampling.
'''
Before sampling from each group randomly to get our training set,
we subset to an arbitrary of the 8 segments, because otherwise we
would not sample isolates but segments also. Note that we sample
groups of equal sizes, to not have one group dominate the prediction.
'''
subset = lu.loc[lu['protein_name'] == 'HA']
# subset = lu.loc[lu['protein_name'] == 'M1']
# subset = lu.loc[lu['protein_name'] == 'NS1']

size = 1000  # sample size
# Counter(lu.host)
# Counter({'Avian': 58632, 'Human': 109872, 'Swine': 22560})
# divide by 8
replace = False  # with replacement

# stackoverflow, 22472213
fn = lambda obj: obj.loc[np.random.choice(obj.index, size, replace), :]
sample = subset.groupby('host', as_index=False).apply(fn)

# train = lu.loc[lu['id'].isin(sample.id)]  # stackoverflow, 18250298
# test = lu.loc[~lu['id'].isin(sample.id)]


# Count.
'''
Include a checking condition to only include "full sets" of 8 segments.
'''


# Is HA sufficient?
# ---------------------------------------------------------------

train = subset.loc[subset['id'].isin(sample.id)]  # stackoverflow, 18250298
test = subset.loc[~subset['id'].isin(sample.id)]

length = 3
shift = 3

# train
l = []
# for key in islice(train.protein_accession, 2):
for key in train.protein_accession:
    s = str(fasta[key])
    c = count_kmer(s[1:], length, shift)  # frame shift: count_kmer(s[1:], 3, 3)
    # normalize
    sum_ = sum(c.values())
    for key in c.keys():
        c[key] = c[key] / sum_
    l.append(c)

# test
u = []
for key in test.protein_accession:
    s = str(fasta[key])
    c = count_kmer(s[1:], length, shift)
    # normalize
    sum_ = sum(c.values())
    for key in c.keys():
        c[key] = c[key] / sum_
    u.append(c)


X_train = pd.DataFrame.from_records(l)
y_train = np.array(train.host)
clf = SVC(decision_function_shape='ovo')  # multiclass classification
clf.fit(X_train, y_train)

X_test = pd.DataFrame.from_records(u)
y_test = clf.predict(X_test)
y_true = np.array(test.host)

pd.crosstab(
    y_true, y_test, rownames=['True'], colnames=['Predicted'], margins=True)

'''
> HA
Predicted  Avian  Human  Swine    All
True
Avian       5806      6     15   5827
Human        148  12014     71  12233
Swine         35     88   1197   1320
All         5989  12108   1283  19380

> M1
Predicted  Avian  Human  Swine    All
True
Avian       5801      6     20   5827
Human        129  11643    461  12233
Swine         22    166   1132   1320
All         5952  11815   1613  19380

Similar results: NS1, NP, NA, PB1, PB2, PA
'''

'''
> nucleotide count (1, 1), HA, n=1000
Predicted  Avian  Human  Swine    All
True
Avian       4014   1455    858   6327
Human        924   1720  10089  12733
Swine        133    102   1585   1820
All         5071   3277  12532  20880

> dinucleotide count (2, 1), HA, n=1000
Predicted  Avian  Human  Swine    All
True
Avian       4294   1048    985   6327
Human        173   8376   4184  12733
Swine         17    630   1173   1820
All         4484  10054   6342  20880

> codon count (3, 3), HA, n=1000
Predicted  Avian  Human  Swine    All
True
Avian       5524    227    576   6327
Human        168  11314   1251  12733
Swine         25    800    995   1820
All         5717  12341   2822  20880

> frameshift codon count (3, 3), HA, n=1000
Predicted  Avian  Human  Swine    All
True
Avian       5519    408    400   6327
Human        201   6152   6380  12733
Swine         34    557   1229   1820
All         5754   7117   8009  20880

> hexamer count (6, 3), HA, n=1000
Predicted  Avian  Human  Swine    All
True
Avian       6303      8     16   6327
Human        205  11302   1226  12733
Swine         61    329   1430   1820
All         6569  11639   2672  20880

> frameshift hexamer count (6, 3), HA, n=1000
Predicted  Avian  Human  Swine    All
True
Avian       6277     15     35   6327
Human        204  11302   1227  12733
Swine         38    343   1439   1820
All         6519  11660   2701  20880
'''

X_train = pd.DataFrame.from_records(l)
y_train = np.array(train.host)
estimator = SVC(kernel='linear', decision_function_shape='ovo')
selector = RFE(estimator, 5, step=1)  # prune to 5 remaining features
selector = selector.fit(X_train, y_train)

X_test = pd.DataFrame.from_records(u)
y_test = selector.predict(X_test)
y_true = np.array(test.host)

pd.crosstab(
    y_true, y_test, rownames=['True'], colnames=['Predicted'], margins=True)

selector.support_  # the number of selected features left after pruning by RFE
selector.ranking_

# X_train.columns[selector.support_]
# Index(['AAA', 'ACC', 'GGA', 'GTA', 'GTG'], dtype='object')
# X_train.columns[selector.ranking_]
'''
array([ 1, 38, 31, 54, 14,  1, 24, 47, 43,  6, 27, 21,  8, 20, 12, 37, 46,
       28,  3, 39, 44, 49, 17, 13, 30, 57, 48, 51, 42,  9, 53, 19, 26, 36,
       16, 25, 59,  4, 50, 18,  1, 33, 35, 45,  1, 10,  1, 11, 60,  5, 58,
       15,  2, 34, 23, 29, 56, 52, 55, 41, 22, 32,  7, 40])

Index(['AAC', 'GCG', 'CTT', 'TCG', 'ATG', 'AAC', 'CGA', 'GTT', 'GGT', 'ACG',
       'CGT', 'CCC', 'AGA', 'CCA', 'ATA', 'GCC', 'GTG', 'CTA', 'AAT', 'GCT',
       'GTA', 'TAC', 'CAC', 'ATC', 'CTG', 'TGC', 'TAA', 'TAT', 'GGG', 'AGC',
       'TCC', 'CAT', 'CGG', 'GCA', 'CAA', 'CGC', 'TGT', 'ACA', 'TAG', 'CAG',
       'AAC', 'GAC', 'GAT', 'GTC', 'AAC', 'AGG', 'AAC', 'AGT', 'TTA', 'ACC',
       'TGG', 'ATT', 'AAG', 'GAG', 'CCT', 'CTC', 'TGA', 'TCA', 'TCT', 'GGC',
       'CCG', 'GAA', 'ACT', 'GGA'],
      dtype='object')
'''
# Note, there are repetitions of features.

# Try to select the optimal number of feautures with CV.
# http://bit.ly/2c5c2ro

X_train = pd.DataFrame.from_records(l)
y_train = np.array(train.host)

estimator = SVC(kernel='linear', decision_function_shape='ovo')
selector = RFECV(estimator=estimator, step=1,
                 cv=StratifiedKFold(y_train, 2),
                 scoring='accuracy')
selector = selector.fit(X_train, y_train)
print('Optimal number of features : %d' % selector.n_features_)

X_test = pd.DataFrame.from_records(u)
y_test = selector.predict(X_test)
y_true = np.array(test.host)

pd.crosstab(
    y_true, y_test, rownames=['True'], colnames=['Predicted'], margins=True)
