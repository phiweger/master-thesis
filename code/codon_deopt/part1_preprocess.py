'''
Aim: Identify "interesting positions" to feed to ML in multiple sequence
alignment.

1. Search database w/ filter applied, select appropriate documents.
2. Multiple sequence alignment, respecting codons.
3. From MSA, create train and test set.
4. For training data, select positions with high entropy.
5. Create feature matrix from positions for train and test set and encode.
6. Learn.
7. Optimize.
8. Assess.
9. Benchmark.

Appendix: Vis.

TODO

- true train, test
- xgboost benchmark (metrics?)
- repeat both to get error estimates + vis
- look at:

clf.feature_importances_  # overall
clf.estimators_[1][1].feature_importances_  # individual trees
# clf.estimators_.shape
# (1000, 3)
# more info: stackoverflow, 17057139
clf.get_params()

- Bayesian model? what we do is classify with SNPs, surely done before
    - specifically: PSD model and scaling via ADVI
    - the SNPs come from GWAS
- repreat for individual segments and whole genome and other viruses
- how overlapping is our entropy based position selection and SNP calling?
- does one of these approaches have a better predictive performance?
- anything known about positions identified (they are of the form: "a
C at this position", "no A at this position")
- better performance with protein sequence (lots more feats)?
- how did our competition do it?
- include all positions and compare to preselecting them
- can we work around MSA prerequisite? maybe something kmer based?
- reduce feautures to how little? i.e. motivate threshold with some
investigation of cost
- tune hyperparams, n of learners?
- at som point entropy as filter might give way to really a position
being important for being this very position, so entropy really is more
of a filter
- mine this material concerning xgboost/ skleran gradient ...:



http://machinelearningmastery.com/data-preparation-gradient-boosting-xgboost-python/

- from sklearn.ensemble import GradientBoostingClassifier
- from xgboost import XGBClassifier

http://datascience.stackexchange.com/questions/10943/why-is-xgboost-so-much-faster-than-sklearn-gradientboostingclassifier

XGBClassifier(base_score=0.5, colsample_bylevel=1, colsample_bytree=1,
gamma=0, learning_rate=0.05, max_delta_step=0, max_depth=10,
min_child_weight=1, missing=None, n_estimators=500, nthread=-1,
objective='binary:logistic', reg_alpha=0, reg_lambda=1,
scale_pos_weight=1, seed=0, silent=True, subsample=1)

GradientBoostingClassifier(init=None, learning_rate=0.05, loss='deviance',
max_depth=10, max_features=None, max_leaf_nodes=None,
min_samples_leaf=1, min_samples_split=2,
min_weight_fraction_leaf=0.0, n_estimators=10,
presort='auto', random_state=None, subsample=1.0, verbose=0,
warm_start=False)


https://www.kaggle.com/c/higgs-boson/discussion/10335

> Indeed, these are very important questions for reproducibility concerns.
XGBoost is a really great piece of software, but some subtle things may or
may not makes it fully comparable with respect to other implementations.
In particular, we have been trying to compare XGBoost and
GradientBoostingClassifier and it happens that the trees that are built are
quite often very different -- which shouldn't be the case if they both are
properly implementing what is called "Gradient Boosted Decision Trees". It
seems many less nodes are often built in XGBoost, as if construction often
terminates early. In many cases, impurity improvements appear to be close
to 0.

> Feature importance is already supported in python version, see get_fscore.

try both, as always, compare speed and accuracy

install:
https://www.kaggle.com/c/melbourne-university-seizure-prediction/discussion/23591

brew update
brew tap homebrew/boneyard
brew install clang-omp  # takes forever, 12:45 - x
pip install xgboost  # in venv lab3, i.e. Pytjon 3.5
# Successfully installed xgboost-0.6a2
# // No edits to .bashrc necessary.


tune parameters:
https://www.analyticsvidhya.com/blog/2016/02/complete-guide-parameter-tuning-gradient-boosting-gbm-python/
'''


from collections import Counter
from itertools import chain
import numpy as np
import os
import pandas as pd
import progressbar
from pymongo import MongoClient
from skbio import TabularMSA, DNA
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
import sys
from xgboost import XGBClassifier


from zoo.utils import random_draw, iter_seq, to_fasta
from zoo.align import msa_subset
from zoo.learn import entropy_select


sys.path.append('/Users/pi/commonplace/compbio/protocols/utils')
from one_hot_encoding import ohe_matrix
from shannon import shannon
from chunk import chunk


'''
0. Database connection, functions, constants etc.
'''


client = MongoClient("localhost:27017")
db = client["zoo"]
collection = db.get_collection('influenza_a_virus')


HOSTS = 'Human Avian Swine'.split(' ')
YEAR = 2000
ANNO = 'NA'  # No described reassortment, see E. Holmes' book.
THRESHOLD = 1.2  # entropy threshold for position selection
SIZE = 0.33
SEED = 42


'''
1. Search database w/ filter applied, select appropriate documents,
write to fasta.
'''


# How many entries exist for each tag (e.g. host)?
match = {
    'annotation.name': ANNO,
    'metadata.host': {'$in': HOSTS},
    'metadata.date.y': {'$gte': YEAR}
    }

group = {
    '_id': '$metadata.host',  # note the "$" prefix
    'cnt': {'$sum': 1}
    }

pipeline = [{'$match': match}, {'$group': group}]
q = collection.aggregate(pipeline)
d = {}
for i in q:
    d[i['_id']] = i['cnt']
# PA: {'Avian': 9945, 'Human': 14791, 'Swine': 3099}


n = min(d.values())  # for a balanced sample
l = []
for h in HOSTS:
    match = {
        'annotation.name': ANNO,
        'metadata.host': h,
        'metadata.date.y': {'$gte': YEAR}
    }
    query = random_draw(collection, match, n)
    l.append(query)
'''
# If run, will consume cursors.
for i in l:
    u = []
    for j in i:
        u.append(j['metadata']['host'])
    print(Counter(u))

# Counter({'Swine': 3099})
# Counter({'Avian': 3099})
# Counter({'Human': 3099})
'''


query = chain(*l)
# counter = 0
# for i in query:
#     if len(i['sequence']) % 3 != 0:
#         counter += 1


'''
# If run, will consume cursors.
counter = 0
for i in query:
    counter += 1

query = get_id(collection, list(X))
'''

fp = '/Users/pi/data/influenza/mod/msa/subset_NA.fa'
to_fasta(fp, iter_seq(query, annotation=ANNO))  # save coding sequence


'''
2. Multiple sequence alignment, respecting codons.
'''


''' shell
python ~/commonplace/compbio/protocols/alignment/msa/codon_msa.py subset.fa
# creates codon.fa and protein.fa

# or if codon alignment not necessary
# TODO: make this happen in temporary file or create run-specific folder
# "analysis_NA_2017-03-13"
mafft --thread -1 subset_NA.fa > subset_NA.mafft.fa
'''


'''
3. From MSA, create train and test set.
'''

fp_msa = '/Users/pi/data/influenza/mod/msa/NA/subset_NA.mafft.fa'
# fp_msa = '/Users/pi/data/influenza/mod/msa/codon.fa'
msa = TabularMSA.read(
    fp_msa, format='fasta', constructor=DNA, lowercase=True)
# msa.shape
# Shape(sequence=9297, position=2151)

# Create a dataframe with ID column and column(s) to stratify the IDs by.
ids = [msa[i].metadata['id'] for i in range(len(msa))]
query = collection.find({'_id': {'$in': ids}})
gen = ((i['_id'], i['metadata']['host']) for i in query)
df = pd.DataFrame.from_records(gen, columns=['id', 'tag'])
# Counter(df['tag'])
# Counter({'Avian': 3099, 'Human': 3099, 'Swine': 3099})


# Split into train and test set.
X_train, X_test, y_train, y_test = train_test_split(
    df['id'], df['tag'],
    test_size=SIZE,
    random_state=SEED,
    stratify=df['tag']
    )
# Counter(y_train)
# Counter({'Avian': 2076, 'Human': 2076, 'Swine': 2076})


'''
We reassign X, y here to make sure their correct mapping is carried.
'''
M, X_train, y_train = msa_subset(msa, X_train, y_train)
N, X_test, y_test = msa_subset(msa, X_test, y_test)


'''
4. For training data, select positions with high entropy.

Why high entropy? The idea is a bit that those are the positions that evolution
can "play with", i.e. where nt substitutions are likely synonymous. So when an
influenza a virus "instance" (=virion) jumps to a new host, these sites are
the buttons that nature can turn in order to make adaptation possible, like the
wheels on a number lock.

Calculate column-wise Shannon entropy, ignoring non-canonical nt.
Write to file for vis w/ R. Note only M, i.e. the training set is used
for this. However it does not make much of a difference, i.e. high
entropy positions in train tend to be high entropy in test.
'''

# training data only (matrix M)

actg = ['A', 'C', 'T', 'G']
l = []

C = np.apply_along_axis(Counter, axis=0, arr=M)
counter = 0
# fp_out = '/Users/pi/data/influenza/mod/msa/conserved_train.tsv'
fp_out = '/Users/pi/data/influenza/mod/msa/NA/conserved_train.tsv'


with open(fp_out, 'w+') as outfile:
    for c in C:
        counts = [value for (key, value) in c.items() if key in actg]
        entropy = shannon(np.array(counts))
        l.append(entropy)
        outfile.write('{}\t{}\t{}\n'.format(
            counter, entropy, sum(counts)
            ))
        counter += 1

# to compare, test data (matrix N)
C = np.apply_along_axis(Counter, axis=0, arr=N)
counter = 0
fp_out = '/Users/pi/data/influenza/mod/msa/conserved_test.tsv'
with open(fp_out, 'w+') as outfile:
    for c in C:
        counts = [value for (key, value) in c.items() if key in actg]
        entropy = shannon(np.array(counts))
        outfile.write('{}\t{}\t{}\n'.format(
            counter, entropy, sum(counts)
            ))
        counter += 1


'''
see R script below for vis of entropy
'''


'''
5. Create feature matrix from positions for train and test set and encode.
'''
entropy = np.array(l)
# selection = np.where(entropy > THRESHOLD)[0]

TOP = 69
selection = entropy_select(entropy, n=50)  # 0.05
# random
selection = random.sample(range(len(entropy)), k=100)
# all
selection = entropy_select(entropy, n=1000)
# really all
selection = np.array(range(len(entropy)))


''' test
a = np.where(entropy > THRESHOLD)[0]
b = entropy_select(entropy, len(a))  # Will select top 69 elements.
len(set(a).intersection(set(b)))  # 69
'''

# random for negative control
# import random
# selection = random.sample(range(len(entropy)), k=69)


# another negative control
# selection = np.where(entropy < 0.1)[0][:69]
# selection = random.sample(set(np.where(entropy < 0.05)[0]), 69)
# selection = random.sample(set(np.where(entropy == 0)[0]), 69)

X_train = ohe_matrix(M[:, selection])
X_test = ohe_matrix(N[:, selection])


'''
6. Learn.
'''


# Fit training data.

# shrinkage and regularization
# http://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_regularization.html#sphx-glr-auto-examples-ensemble-plot-gradient-boosting-regularization-py

# n_estimators
# http://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_oob.html#sphx-glr-auto-examples-ensemble-plot-gradient-boosting-oob-py

# feature transformations
# http://scikit-learn.org/stable/auto_examples/ensemble/plot_feature_transformation.html#sphx-glr-auto-examples-ensemble-plot-feature-transformation-py


NE = 3000
clf = GradientBoostingClassifier(
    n_estimators=NE,
    learning_rate=0.1,
    subsample=0.5,
    random_state=SEED,
    verbose=True
    ).fit(X_train, y_train)


# os.system('afplay /System/Library/Sounds/Sosumi.aiff')


# Predict test data.
y_true = np.array(y_test)
y_pred = clf.predict(X_test)


# assess
pd.crosstab(
    pd.Series(y_true),
    pd.Series(y_pred),
    rownames=['True'],
    colnames=['Predicted'],
    margins=True)

'''
Predicted  Avian  Human  Swine   All
True
Avian       1006      5     12  1023
Human         11    991     21  1023
Swine          9    119    895  1023
All         1026   1115    928  3069
'''


'''
How stable?
'''

# selection = random.sample(range(len(entropy)), k=100)
selection = np.array(range(len(entropy)))
X_train = ohe_matrix(M[:, selection])
X_test = ohe_matrix(N[:, selection])

# NE = 3000
# l = []
# for i in range(20):
#     clf = GradientBoostingClassifier(
#         n_estimators=NE,
#         learning_rate=0.1,
#         subsample=0.5,
#         # random_state=SEED,
#         verbose=True
#         ).fit(X_train, y_train)
#     l.append(clf.feature_importances_)


'''
# scaling unnecessary, feature importance already scaled
# stackoverflow, 26414913
from sklearn import preprocessing
x = df.values #returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
df = pd.DataFrame(x_scaled)
'''


''' R
library(ggplot2)

fp <- '/Users/pi/data/influenza/mod/feature_importance/PA_random100_repeat20.csv'
df <- read.table('~/tmp/fi_rep.csv', header=F, stringsAsFactors=F, sep=',')
names(df) <- c('pos', 'rep', 'val')

p <-
ggplot(df, aes(x=as.factor(rep), y=val)) +
    geom_boxplot(outlier.size=NA) +
    scale_y_log10() +
    theme_default() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        ) +
    ylab('feature importance') +
    xlab('position')
    # scale_y_continuous(limits=c(0, 0.1))

fp_img = '/Users/pi/projects/influenza/img/feature_importance/PA_random100_repeat20_log.png'
ggsave(fp_img, p, width=20, height=5, units='cm')
'''


'''
7. Optimize.

- see 28:48 here: https://www.youtube.com/watch?v=IXZKgIsZRm0
- more
'''


# Boosted trees
# about 100 parameter combinations
param_grid = {
    'learning_rate': [0.1, 0.05, 0.02, 0.01],
    'max_depth': [4, 6],
    'min_samples_leaf': [3, 5, 9, 17],
    'max_features': [1.0, 0.3, 0.1]
}
est = GradientBoostingClassifier(n_estimators=1000)
%time gs_cv = GridSearchCV(est, param_grid).fit(X_train, y_train)  # %time
'''
CPU times: user 7h 33min 57s, sys: 2min 55s, total: 7h 36min 53s
Wall time: 7h 46min 5s

gs_cv.best_params_
{'learning_rate': 0.02,
 'max_depth': 4,
 'max_features': 0.1,
 'min_samples_leaf': 3}
'''

m1 = GradientBoostingClassifier(n_estimators=1000).fit(X_train, y_train)
m2 = GradientBoostingClassifier(
    n_estimators=1000,
    learning_rate=0.02,
    max_depth=4,
    max_features=0.1,
    min_samples_leaf=3).fit(X_train, y_train)
# clf = gs_cv.best_estimator_
m3 = gs_cv.best_estimator_.fit(X_train, y_train)

m1.score(X_test, y_test)  # 0.927
m2.score(X_test, y_test)  # 0.927
m3.score(X_test, y_test)  # 0.925


'''
8. Assess.
'''

# http://machinelearningmastery.com/feature-importance-and-feature-selection-with-xgboost-in-python/


# clf.fit(X_train, y_train)
# clf.feature_importances_
# np.argsort(clf.feature_importances_)[::-1]
# y_true = y_test
# y_pred = clf.predict(X_test)


clf.feature_importances_  # overall
chunks = chunk(4, clf.feature_importances_)
feat_imp = ['%.3f' % round(sum(i), 3) for i in chunks]  # stackoverflow, 56820
feat_ent = entropy[selection]

fp_out = '/Users/pi/data/influenza/mod/msa/feat_importance_entropy.tsv'
with open(fp_out, 'w+') as outfile:
    outfile.write('{}\t{}\n'.format('entropy', 'importance'))
    for i in zip(feat_ent, feat_imp):
        outfile.write('{}\t{}\n'.format(i[0], i[1]))







clf.estimators_[1][1].feature_importances_  # individual trees
# clf.estimators_.shape
# (1000, 3)
# more info: stackoverflow, 17057139
clf.get_params()

# from sklearn.metrics import classification_report
# from sklearn.metrics import roc_curve, auc

# stackoverflow, 29148355
# overall accuracy
acc = clf.score(X_test, y_test)  # 0.94



# get roc/auc info
# Y_score = clf.decision_function(X_test)
# fpr = dict()
# tpr = dict()
# fpr, tpr, _ = roc_curve(y_test.transpose(), Y_score)

# roc_auc = dict()
# roc_auc = auc(fpr, tpr)


'''
# header: number of selected positions, accuracy
# Note that a suffix like 69r means "random draw", i.e. negative control.
# rl .. random and low entropy
69    0.943
69r   0.920
69rl  0.358  threshold < 0.05
69rl  0.331  threshold == 0

54    0.931    1000 learner
      0.929    500
      0.928    250
      0.924    125
50    0.928

'''

'''
check if correlation feature importance and entropy when learning on all
feats

hypothesis: by simply changing highly entropic sites, virus cannot react
(in the codon optimization field)
'''

len(clf.feature_importances_) == X_train.shape[1]


'''
sum feature importance accross 4 one hot encoding sites, then compare
visually to entropy
'''


NE = 125
fp = ''
with open(fp, 'w+') as outfile:

    selection = entropy_select(entropy, n=0.05)
    X_train = ohe_matrix(M[:, selection])
    X_test = ohe_matrix(N[:, selection])
    # Fit training data.
    clf = GradientBoostingClassifier(n_estimators=NE)
    clf.fit(X_train, y_train)
    # Predict test data.
    y_true = np.array(y_test)
    y_pred = clf.predict(X_test)


'''
7. Benchmark.

- http://xgboost.readthedocs.io/en/latest/model.html
'''


model = XGBClassifier()
model.fit(X_train, y_train)
y_true = np.array(y_test)
y_pred = model.predict(X_test)
pd.crosstab(
    pd.Series(y_true),
    pd.Series(y_pred),
    rownames=['True'],
    colnames=['Predicted'],
    margins=True)

'''
Predicted  Avian  Human  Swine   All
True
Avian       1011      3      9  1023
Human         12    992     19  1023
Swine         15    155    853  1023
All         1038   1150    881  3069
'''


'''
Appendix: Vis.
'''


''' R

library(ggplot2)

fp_train <- '~/data/influenza/mod/msa/conserved_train.tsv'
df_train <- read.table(fp_train, sep='\t', stringsAsFactors=F, header=F)
names(df_train) = c('index', 'entropy', 'sum')
df_train$label = 'training'

fp_test <- '~/data/influenza/mod/msa/conserved_test.tsv'
df_test <- read.table(fp_test, sep='\t', stringsAsFactors=F, header=F)
names(df_test) = c('index', 'entropy', 'sum')
df_test$label = 'test'

df = rbind(df_train, df_test)


theme_default <- function(
    base_size = getOption("bayesplot.base_size", 12),
    base_family = getOption("bayesplot.base_family", "serif")) {

        theme_bw(base_family = base_family, base_size = base_size) +
        theme(
            plot.background = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(size = 0.4),
            axis.ticks = element_line(size = 0.3),
            strip.background = element_blank(),
            strip.text = element_text(size = rel(0.9)),
            strip.placement = "outside",
            # strip.background = element_rect(fill = "gray95", color = NA),
            panel.spacing = unit(1.5, "lines"),
            legend.position = "right",
            legend.background = element_blank(),
            legend.text = element_text(size = 13),
            legend.text.align = 0,
            legend.key = element_blank()
        )
    }


# entropy vis
ggplot(df[df$label == 'training',], aes(x=entropy)) +
    geom_freqpoly() +
    theme_default()
fp_out = '~/projects/influenza/img/entropy_PA_freqpoly.pdf'
ggsave(fp_out, height=5, width=5, unit='cm')

# ggplot(df, aes(x=entropy, color=label)) +
#     geom_freqpoly() +
#     scale_color_manual(
#         name="split",
#         values=c(test="#DCBCBC", training="#8F2727")) +
#     theme_default()


ggplot(df[df$label == 'training',], aes(x=index, y=entropy)) +
    geom_point(size=0.1) +
    theme_default()
fp_out = '~/projects/influenza/img/entropy_PA_train_test.pdf'
ggsave(fp_out, height=5, width=7.5, unit='cm')


delta <- as.data.frame(df_test$entropy - df_train$entropy)
names(delta) <- "delta"
ggplot(delta, aes(x=delta)) +
    geom_histogram(fill='white', color='black') +
    theme_default()
fp_out = '~/projects/influenza/img/entropy_PA_delta.pdf'
ggsave(fp_out, height=5, width=7.5, unit='cm')


# feature importance
fp <- '~/data/influenza/mod/msa/feat_importance_entropy.tsv'
df <- read.table(fp, sep='\t', stringsAsFactors=F, header=T)
ggplot(df, aes(x=entropy, y=importance)) +
    geom_point(size=0.5) +
    scale_y_log10() +
    theme_default()
fp_out = '~/projects/influenza/img/entropy_PA_feat_importance.pdf'
ggsave(fp_out, height=5, width=5, unit='cm')

'''






# ----------------------------------------------------------------------------

fp_msa = '/Users/pi/data/influenza/mod/msa/deprecated/train/codon.fa'
msa = TabularMSA.read(
    fp_msa, format='fasta', constructor=DNA, lowercase=True)


'''
Create a lookup table mapping an ID (accession number) to a host.
'''
ids = [msa[i].metadata['id'] for i in range(len(msa))]
hosts = []
for i in ids:
    cursor = collection.find({'_id': i})
    hosts.append(next(cursor)['metadata']['host'])
lookup = {key: value for (key, value) in zip(ids, hosts)}


'''
Check out same wird samples: CY159618, CY135240, CY159920, CY135194, ...
They have large stretches of "N" characters in their sequences. This
prevents us from using this:
'''

conserved = msa.conservation(gap_mode='nan', degenerate_mode='nan')
consensus = msa.consensus()

'''
Manually coding this w/o the introduction of all the nans.
'''

actg = ['A', 'C', 'T', 'G']
counter = 0
fp_out = '/Users/pi/data/influenza/mod/msa/deprecated/train/conserved.tsv'
bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)


with open(fp_out, 'w+') as outfile:
    outfile.write('index\tshannon\tsum\n')  # header

    for i in msa.iter_positions(ignore_metadata=True):
        counts = [
            value for (key, value) in Counter(str(i)).items() if key in actg
            ]
        outfile.write('{}\t{}\t{}\n'.format(
            counter, shannon(np.array(counts)), sum(counts)
            ))
        counter += 1
        bar.update(counter)


''' R

library(ggplot2)

fp_train <- '~/data/influenza/mod/msa/conserved_train.tsv'
df_train <- read.table(fp_train, sep='\t', stringsAsFactors=F, header=F)
names(df_train) = c('index', 'entropy', 'sum')
df_train$label = 'training'

fp_test <- '~/data/influenza/mod/msa/conserved_test.tsv'
df_test <- read.table(fp_test, sep='\t', stringsAsFactors=F, header=F)
names(df_test) = c('index', 'entropy', 'sum')
df_test$label = 'test'

df = rbind(df_train, df_test)



theme_default <- function(
    base_size = getOption("bayesplot.base_size", 12),
    base_family = getOption("bayesplot.base_family", "serif")) {

        theme_bw(base_family = base_family, base_size = base_size) +
        theme(
            plot.background = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(size = 0.4),
            axis.ticks = element_line(size = 0.3),
            strip.background = element_blank(),
            strip.text = element_text(size = rel(0.9)),
            strip.placement = "outside",
            # strip.background = element_rect(fill = "gray95", color = NA),
            panel.spacing = unit(1.5, "lines"),
            legend.position = "right",
            legend.background = element_blank(),
            legend.text = element_text(size = 13),
            legend.text.align = 0,
            legend.key = element_blank()
        )
    }

ggplot(df[df$entropy>0,], aes(x=index, y=entropy, color=entropy)) +
    #geom_hline(yintercept=1, color='#8F2727') +
    geom_point(size=0.5) +
    scale_color_gradient(high='#8F2727', low='#DCBCBC') +
    theme_default()

fp_out = '~/projects/influenza/img/entropy_PA.pdf'
ggsave(fp_out, height=5, width=10, unit='cm')


ggplot(df[df$entropy>0,], aes(x=entropy)) +
    geom_freqpoly() +
    theme_default()

fp_out = '~/projects/influenza/img/entropy_PA_freqpoly.pdf'
ggsave(fp_out, height=5, width=5, unit='cm')


ggplot(df, aes(x=index, y=entropy)) +
    #geom_hline(yintercept=1, color='#8F2727') +
    geom_point(size=0.5, color='#DCBCBC') +
    geom_point(
        data=df[df$entropy>1.2 & df$label=='test',],
        color='#8F2727', size=0.5) +
    geom_point(
        data=df[df$entropy>1.2 & df$label=='training',],
        color='#B97C7C', size=0.5) +
    #scale_color_gradient(high='#8F2727', low='#DCBCBC') +
    theme_default()

fp_out = '~/projects/influenza/img/entropy_PA_train_test.pdf'
ggsave(fp_out, height=5, width=10, unit='cm')
'''


# Interesting thresholds are 0.5 and 1.2.
THRESHOLD = 1.2


'''
We want a feature matrix X

columns are features = positions in MSA
rows are samples = nt sequences indexed by GenBank accession id

sample/ pos 1 2 3
a A T G
b G T G
c A T C

Furthermore we need a target vector y, holding a label for each sample

sample label
a Avian
b Human
c Avian
'''


'''
X
'''


fp_conserved = '/Users/pi/data/influenza/mod/msa/deprecated/train/conserved.tsv'
entropy = []
with open(fp_conserved, 'r+') as infile:
    file = infile.readlines()[1:]  # skip header
    for line in file:
        e = float(line.split('\t')[1])  # entropy val, cast to float
        entropy.append(e)

entropy = np.array(entropy)
selection = np.where(entropy > THRESHOLD)[0]


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
d = {}
index = 0
for i in msa.iter_positions(ignore_metadata=True):
    if index in selection:
        d[index] = list(str(i))  # pos index is key, value are nt at pos
    index += 1
    bar.update(index)
# len(d.keys()) == len(selection)  # True


# X = pd.DataFrame.from_dict(d)
X = np.array(list(d.values())).transpose()


'''
y
'''


y = np.array(hosts)


'''
Label encode y. One hot encode X.
'''
# from sklearn.preprocessing import LabelEncoder, OneHotEncoder

from sklearn.model_selection import train_test_split
# from sklearn.grid_search import GridSearchCV
from sklearn.ensemble import GradientBoostingClassifier


# label_encoder = LabelEncoder()
# label_encoder = label_encoder.fit(y)
# label_encoded_y = label_encoder.transform(y)


'''
We have non-canonocale, i.e. non {ACTG} entries. We could filter them.
However, XGBoost, which we will be using, does not care. Since we have to
one hot encode the nucleotides, an "A" will become [1, 0, 0, 0], while
R/Y/N/foobar will get [0, 0, 0, 0].

sources:

    - https://github.com/dmlc/xgboost/issues/21
    - https://arxiv.org/abs/1603.02754#
'''




X = ohe_matrix(X)

n_estimators = 1000
# Create test and train set.
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.33, random_state=42)
# ML
clf = GradientBoostingClassifier(n_estimators=n_estimators)
clf.fit(X_train, y_train)
# clf.fit(X, y)

y_true = y_test
y_pred = clf.predict(X_test)

pd.crosstab(
    pd.Series(y_true),
    pd.Series(y_pred),
    rownames=['True'],
    colnames=['Predicted'],
    margins=True)


'''
Predicted  Avian  Human  Swine   All
True
Avian        691      0      5   696
Human          6    656     11   673
Swine          5     91    591   687
All          702    747    607  2056
'''


'''
The real test data.
'''

fp_msa = '/Users/pi/data/influenza/mod/msa/test/codon.fa'
msa = TabularMSA.read(
    fp_msa, format='fasta', constructor=DNA, lowercase=True)


ids = [msa[i].metadata['id'] for i in range(len(msa))]
hosts = []
for i in ids:
    cursor = collection.find({'_id': i})
    hosts.append(next(cursor)['metadata']['host'])


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
d = {}
index = 0
for i in msa.iter_positions(ignore_metadata=True):
    if index in selection:
        d[index] = list(str(i))  # pos index is key, value are nt at pos
    index += 1
    bar.update(index)
# len(d.keys()) == len(selection)  # True


X_test = np.array(list(d.values())).transpose()
X_test = ohe_matrix(X_test)
y_test = np.array(hosts)


y_true = y_test
y_pred = clf.predict(X_test)


pd.crosstab(
    pd.Series(y_true),
    pd.Series(y_pred),
    rownames=['True'],
    colnames=['Predicted'],
    margins=True)


