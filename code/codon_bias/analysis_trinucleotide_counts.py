# TODO: Plot the number of complete sequences per IVA subtype

from Bio.SeqIO.FastaIO import FastaIterator
from itertools import product, repeat
import pandas as pd
import re


path_df = '.../genomeset.tsv'
df = pd.read_csv(path_df, header=None, index_col=False, delimiter='\t')
df.columns = ['accession host segment subtype country date length description age gender id'.split(' ')]

# How many entries are there for each unique IVA subtype?
df.subtype.value_counts()

# Eliminate all entries whose len(set(internal ID == x)) %% 8 != 0,
# i.e. all entries from which we do not have 8 sequences.

# In [40]: len(df)
# Out[40]: 225277
for i in df.id.unique():
    entries = df[df.id == i]
    n_entries = len(entries)
    if n_entries != 8:
        df = df.drop(entries.index)
# In [42]: len(df)
# Out[42]: 222568

df.subtype.value_counts() / 8
# Still not perferct, but it'll do.

# Include only subtypes with > 10 complete genome entries.
io = df.subtype.value_counts() / 8 > 10
include_subtype = set(io[io == True].index)
include_index = df.subtype.isin(include_subtype)
df = df[include_index]
# In [85]: len(df)
# Out[85]: 221159

# Filter host.
df = df[df.host.isin({'Human', 'Avian', 'Swine'})]
# In [89]: len(df)
# Out[89]: 215487


# ------------------------------------------------------------------
# TODO

path_fa = '.../influenza_genome_set.fa'
# TODO: adjust pattern, old pattern = '.*\|(.*?)\:.*'
# deprecated: pattern = '.*\s\d+\s.*\|(.*?)\:.*'
# deprecated: pattern = '.*\s\d+\s.*\|(.*?)\:(\d+)\-(\d+)'
pattern = '.*\s\d+\s.*\|(.*?)\:[<]?(\d+)[->]+(\d+)'
# a = 'cds:BAB39505 13383267 gb|AB049153:1-2280'
# re.search('.*\s\d+\s.*\|(.*?)\:.*', a).group(1)
# 'AB049153'
# shit to deal with: weird signs around range numbers
# 'cds:AAK51352 13991402 gb|AF320067:1->1032'
# 'cds:BAK24091 333778340 gb|AB629711:<11-2284'


include_accession = set(df.accession.unique())


def filter_influenza_fa(in_fasta, out_fasta, pattern, accession_set):
    '''
    accession_set .. a set of accession IDs that we query
    the fasta header against

    l, count = filter_influenza_fa(path_fa, pattern, include_accession)
    '''
    cache_previous = ()
    count, l = 0, []

    with open(in_fasta) as handle, open(out_fasta, 'a+') as out:
        for record in FastaIterator(handle):
            # [^1]
            if '(' not in list(record.description):
                cache_current = re.search(
                    pattern, record.description).group(1, 2, 3)
                if cache_current[0] in cache_previous:
                    # [^2]
                    continue
                acc = cache_current[0]
                cache_previous = cache_current
                if acc in accession_set:
                    count += 1
                    l.append(acc)
                    out.write('>' + acc + '\n')
                    out.write(str(record.seq) + '\n')
    return(count)  # l could be returned also





# ------------------------------------------------------------------
# deprecated

# pattern = re.compile('\(.*\|.*?\:.*')  # see comment below

# count = 0
# l = []
# include_accession = set(df.accession.unique())
# with open(path_fa) as handle:
#     for record in FastaIterator(handle):
#         accession = re.search('\w\:(.*)', record.id).group(1)
#         if accession in include_accession:
#             count += 1
#             l.append(accession)

# count = 0
# l = []
# include_accession = set(df.accession.unique())

# outfile = '/Users/pi/tmp/result.tsv'
# with open(outfile, 'a+') as out:
#     with open(path_fa) as handle:
#         for record in FastaIterator(handle):
#             accession = re.search('.*\|(.*?)\:.*', record.id).group(1)
#             if pattern.match(record.id):
#                 continue
#             if accession in include_accession:
#                 count += 1
#                 l.append(accession)
#                 print(record.id)

#                 result = count_tri(str(record.seq))
#                 out.write(accession + '\t' + '\t'.join(result) + '\n')

# In [124]: len(l)
# Out[124]: 302072
# Duplicates?
# In [125]: len(set(l))
# Out[125]: 215406
# Yes.
# Added loop to eliminate () duplicates. However, some remain.
# data$ grep "CY041248" influenza.cds
# >gb|CY041248:15-2291|Influenza A virus (A/mallard/Netherlands/3/2005(H3N8))
#   segment 2, complete sequence
# >gb|CY041248:109-381|Influenza A virus (A/mallard/Netherlands/3/2005(H3N8))
#   segment 2, complete sequence
# Because we have to get on, I'll leave this for the moment at that.

# ------------------------------------------------------------------


# Sampling.

df_avian = df[df.host == 'Avian']
# sample 200 internal IDs

sample_id_avian = set(df_avian.id.sample(200))
# internal ID guarantees that we get 8 segments per ID

sample_df_avian = df_avian[df_avian.id.isin(sample_id_avian)]
# Because each ID represents 8 segments, we should find 8 * 200 =
# 1600 rows in this df.
# In [260]: len(sample_df_avian)
# Out[260]: 1584
# Close enough.

# TODO: repeat for human and swine
# ALso: just do all, no sampling required because file small.

# Counting.

# Now count triplets etc., of the form "host | count 1st codon combination
# | 2nd ..."
# Better: have a separate file for each host

# human_3.tsv
# human_2.tsv
# numbers refer to triplets etc.


# Normalizing.
# by total count per sample
# use log scale?

# This goes into SVM and ggplot.



import numpy as np
from sklearn.decomposition import PCA
X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
pca = PCA(n_components=2)  # n_components=’mle’, TP Minka
pca.fit(X)  # gets vectors
pca.fit_transform(X)  # gets transformed datapoints
# In [9]: pca.fit_transform(X)
# Out[9]:
# array([[-1.38340578,  0.2935787 ],
#        [-2.22189802, -0.25133484],
#        [-3.6053038 ,  0.04224385],
#        [ 1.38340578, -0.2935787 ],
#        [ 2.22189802,  0.25133484],
#        [ 3.6053038 , -0.04224385]])

# write to file as host,x,y and ggplot aes(x=x, y=y, color=host)


# SVM
# http://scikit-learn.org/stable/modules/svm.html


# Is the CDS divisable by 3 as logic dictates?

path_set = '/Users/pi/projects/influenza/data/influenza_genome_set.fa'
l = []
with open(path_set) as handle:
    for record in FastaIterator(handle):
        # count
        acc = record.id
        seq = str(record.seq)
        l.append(len(seq) % 3)
ll = [0 if i==0 else 1 for i in l]
# In [391]: Counter(ll)
# Out[391]: Counter({0: 317618, 1: 283})
# for subset
# Out[395]: Counter({0: 209050, 1: 1})
# Goodie.


# Counting.

path_subset = '/Users/pi/projects/influenza/data/influenza_genome_subset.fa'
with open(path_subset) as handle, open('/Users/pi/tmp/fooluu.tsv', 'a+') as out:
    for record in FastaIterator(handle):
        # count
        acc = record.id
        seq = str(record.seq)
        d = count_tri(seq)
        # What about triplets like GRT, ANC etc.?
        header = [key for key, value in sorted(d.items())]
        result = [str(value) for key, value in sorted(d.items())]
        out.write(acc + '\t' + '\t'.join(result) + '\n')
# Don't normalize at this step, as it would throw away the information
# about the total number of counts.
# output: AB049162  8   9   8   16  7 ...
header = header


# Read in count data and link to df.

df2 = pd.read_csv('/Users/pi/tmp/fooluu.tsv',
    delimiter='\t', header=None, index_col=None)
df2.columns = ['accession'] + header


dfdf2 = pd.merge(df, df2, how='inner', on=['accession'])
grouped = dfdf2.groupby('id')
sum_ = grouped.sum()
annotation = grouped.first()
# Because the first entry for most columns is equal to the other 7
# in the same ID group.
# In [451]: len(grouped.first())
# Out[451]: 26141
# In [452]: len(grouped.sum())
# Out[452]: 26141
# Note that the only columns we cannot use are
# * segment
# * accession
# as they only shows the first genome fragments for all rows.
sum_.drop(['segment', 'length'], axis=1, inplace=True)


df_count = pd.concat([annotation.host, annotation.subtype, sum_], axis=1)
# sum_['index'] = sum_.index

df_count.to_csv(
    '/Users/pi/projects/influenza/results/count_trinucleotide_annotated.tsv',
    sep='\t', header=True, index=None)


# Create dataframe for SVM.
# 1. Select required data.
# 2. Split in train and test.
# 3. Fit.

ah = annotation.host
# # x = df3.icol(list(range(8, 137)))  # icol is deprecated
# x = df3.iloc[:, 8:137]  # like working with dataframes in R

# Do the indices match for x (here sum_) and y?
# In [552]: Counter(sum_.index == y.index)
# Out[552]: Counter({True: 26141})

# from sklearn import svm
# X = sum_.iloc[0:10, :]
# y = ah.iloc[0:10]
# clf = svm.SVC()
# clf.fit(X, y)

# In [569]: clf.predict(X.iloc[1,:])
# > DeprecationWarning: Passing 1d arrays as data is deprecated in 0.17
# and willraise ValueError in 0.19. Reshape your data either using
# X.reshape(-1, 1) if your data has a single feature or X.reshape(1, -1)
# if it contains a single sample.

# clf.predict(X.iloc[1,:].reshape(-1, 64))
# Out[571]: array(['Human'], dtype=object)


# Sample.
# http://stackoverflow.com/questions/35346421/pandas-random-sample-with-ration-11-of-specific-column-entry
# TODO: write a function that automatizes this,
# subsample(proportions=[0.1, 0.4, 0.5])
# needs to match number of labels and sum to 1

# Create views of label.
human = ah.loc[ah == 'Human']
avian = ah.loc[ah == 'Avian']
swine = ah.loc[ah == 'Swine']

human_train = human.sample(5000, replace=False)
avian_train = avian.sample(5000, replace=False)
swine_train = swine.sample(5000, replace=True)  # len(swine) 2841

host_train = pd.concat([human_train, avian_train, swine_train], axis=0)
sum_train = sum_.ix[host_train.index]

sa = set(ah.index)  # len 26141
sb = set(host_train.index)  # len 12364, remember, swine with replacement
test_index = sa.difference(sb)  # len 13777

host_test = ah.ix[test_index]
sum_test = sum_.ix[test_index]


# ML.

from sklearn import svm
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.metrics import precision_score, recall_score, f1_score

# CLassify.

X = sum_train
y = host_train
clf = svm.SVC(decision_function_shape='ovo')  # multiclass classification
clf.fit(X, y)

clf.predict(sum_test)

y_true = host_test
y_pred = clf.predict(sum_test)
pd.crosstab(y_true, y_pred,
    rownames=['True'], colnames=['Predicted'], margins=True)
# Predicted  Avian  Human  Swine    All
# True
# Avian       4547      7      6   4560
# Human        435   8275     30   8740
# Swine        171      8    298    477
# All         5153   8290    334  13777


# Evaluate some more.
confusion_matrix(y_true=y_true, y_pred=y_pred,
    labels=['Human', 'Avian', 'Swine'])
accuracy_score(y_true=y_true, y_pred=y_pred)  # 0.95

# precision_score(y_true=y_true, y_pred=y_pred)
# DeprecationWarning: The default `weighted` averaging is deprecated, and
# from version 0.18, use of precision, recall or F-score with multiclass or#
# multilabel data or pos_label=None will result in an exception. Please set an
# explicit value for `average`, one of (None, 'micro', 'macro', 'weighted',
# 'samples'). In cross validation use, for instance, scoring="f1_weighted"
# instead of scoring="f1".


recall_score(y_true=y_true, y_pred=y_pred,
    average=None, labels=['Avian', 'Human', 'Swine'])
# [ 0.99714912,  0.94679634,  0.62473795]
precision_score(y_true=y_true, y_pred=y_pred,
    average=None, labels=['Avian', 'Human', 'Swine'])
# [ 0.8823986 ,  0.99819059,  0.89221557]
f1_score(y_true=y_true, y_pred=y_pred,
    average=None, labels=['Avian', 'Human', 'Swine'])
# [ 0.93627098,  0.97181445,  0.73489519]


from sklearn.decomposition import PCA

X = sum_
pca = PCA(n_components=2)  # 'mle'
pca.fit(X)
pca_df = pd.DataFrame(
    pca.fit_transform(X)
    )
pca_df['host'] = list(annotation.host)
pca_df['subtype'] = list(annotation.subtype)

pca_df.to_csv('/Users/pi/projects/influenza/results/pca_trinucleotide.tsv',
    sep='\t', header=True, index=None)
