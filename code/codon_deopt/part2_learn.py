'''
transition matrix, codon, codon bias
'''


import sys
sys.path.append('/Users/pi/commonplace/compbio/protocols/utils')
from chunk import chunk

# sum(clf.feature_importances_ > 0.01)


def feature_to_transmat(seq, feature_importance, threshold, cds=True, shift=3):
    '''
    Backtranslate the feature importances of the one hot encoded
    DNA to a transition matrix.

    returns transition matrix, codon changes

    codon changes as dict: {'ACT': 'CTG'}, argument "shift" to allow
    dinucleotide an bicodon bias investigation?

    build on this function: is change synonymous? codon bias change

    note that we have to deal w/ case where neiboring positions in same codon
    best use some sliding window approach, asking every time is this position
    > theshold
    '''
    pass

'''
Write out feature importance.
'''
a = list(chunk(4, clf.feature_importances_))
df = pd.DataFrame.from_records(a)
df.columns = 'A C T G'.split(' ')
df.unstack().reset_index().to_csv(
    '/Users/pi/data/influenza/mod/feature_importance/fi_NA.csv',
    header=None, index=False)


'''
Write out entropy.
'''
s = pd.Series(entropy)
s.to_csv(
    '/Users/pi/data/influenza/mod/feature_importance/entropy_NA.csv',
    index=True, header=None)





# for superheat, see below
# pd.DataFrame.from_records(a).to_csv(
#     '/Users/pi/tmp/fi_squares.csv', header=None, index=False)

''' R
library(ggplot2)
# library(superheat)  # https://github.com/rlbarter/superheat

fp <- '/Users/pi/data/influenza/mod/feature_importance/fi_square_PA.csv'
df <- read.table(fp, header=F, stringsAsFactors=F, sep=',')
names(df) <- c('nt', 'pos', 'val')

p <-
ggplot(df, aes(x=pos, y=nt, fill=log(val))) +
    geom_tile() +
    theme_default() +
    theme(legend.key.size = unit(0.2, "cm")) +
    xlab('position') +
    ylab('nt') +
    scale_fill_continuous(
        na.value = 'white',
        low='white',
        high='black',
        guide = guide_legend(
            title = 'FI',
            direction = "horizontal",
            title.position = "top",
            label.position="bottom",
            label.hjust = 0.5,
            label.vjust = 0.5,
            label.theme = element_text(angle = 90)
            )
        )

q <-
ggplot(df, aes(x=pos, y=sqrt(val))) +
    geom_point(size=0.2) +
    facet_wrap(~ nt, ncol=1) +
    theme_default()

r <-
ggplot(df, aes(x=pos, y=val, color=as.factor(nt))) +
    geom_point(size=0.5) +
    theme_default() +
    xlab('position') +
    ylab('FI') +
    #ylab(expression(sqrt(FI))) +
    scale_color_discrete(guide = guide_legend(title = 'nt'))


fp_out = '~/projects/influenza/img/feature_importance/points_PA.pdf'
ggsave(fp_out, r, height=5, width=10, unit='cm')
'''


'''
How stale prediction depending on learning rate (i.e. regularization)?
'''


NE = 1000
l = []
for lr in [0.5, 0.25, 0.1, 0.05, 0.025, 0.01]:
    print('lr:', lr)
    clf = GradientBoostingClassifier(
        n_estimators=NE,
        learning_rate=lr,
        subsample=0.5,
        # random_state=SEED,
        verbose=True
        ).fit(X_train, y_train)
    l.append(clf.feature_importances_)


fp = '/Users/pi/data/influenza/mod/feature_importance/PA_regularized.csv'
df = pd.DataFrame.from_records(l)

# df.transpose().unstack().reset_index().to_csv(
#     fp, index=False, header=None)


rates = []
lr = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01]
counter = 0
for i in df.iterrows():
    a = list(chunk(4, i[1]))
    b = pd.DataFrame.from_records(a)
    b.columns = 'A C T G'.split(' ')
    c = b.unstack().reset_index()
    c['lr'] = lr[counter]
    counter += 1
    rates.append(c)

result = pd.concat(rates)
result.to_csv(fp, header='nt pos val lr'.split(' '), index=False)


''' R
library(ggplot2)

fp <- '/Users/pi/data/influenza/mod/feature_importance/PA_regularized.csv'
df <- read.table(fp, header=T, stringsAsFactors=F, sep=',')

p <-
ggplot(df, aes(x=pos, y=val, color=as.factor(nt))) +
    geom_point(size=0.5) +
    theme_default() +
    xlab('position') +
    ylab('FI') +
    #ylab(expression(sqrt(FI))) +
    scale_color_discrete(guide = guide_legend(title = 'nt')) +
    facet_wrap(~ as.factor(lr), ncol=1, scale='free')


fp_out = '~/projects/influenza/img/feature_importance/points_PA_regularized.pdf'
ggsave(fp_out, p, height=20, width=8, unit='cm')
'''


''' R
# Plot entropy and feature importance for NA protein.

library(ggplot2)
library(dplyr)


fp_f <- '/Users/pi/data/influenza/mod/feature_importance/fi_NA.csv'
fp_h <- '/Users/pi/data/influenza/mod/feature_importance/entropy_NA.csv'
fp_m <- '/Users/pi/data/influenza/mod/feature_importance/mismatch_DL_NA.csv'

df_h <- read.table(fp_h, header=F, stringsAsFactors=F, sep=',')
names(df_h) <- c('pos', 'h')
df_f <- read.table(fp_f, header=F, stringsAsFactors=F, sep=',')
names(df_f) <- c('nt', 'pos', 'fi')
df_m <- read.table(fp_m, header=F, stringsAsFactors=F, sep=',')
 names(df_m) <- c('pos', 'mismatch')

# x: position, y: feature importance, color: entropy, facet: nucleotide


p <-
ggplot(df_f, aes(x=pos, y=sqrt(fi), color=as.factor(nt))) +
    geom_point(size=0.5) +
    theme_default() +
    xlab('i') +
    ylab('F') +
    #ylab(expression(sqrt(FI))) +
    scale_color_discrete(guide = guide_legend(title = 'nt'))

q <-
ggplot(df_h, aes(x=pos, y=h)) +
    geom_point(size=0.5) +
    theme_default() +
    xlab('i') +
    ylab('H')

df <- inner_join(df_h, df_f, by='pos')

r <-
ggplot(df, aes(x=pos, y=sqrt(fi), color=h)) +
    geom_point(size=0.3) +
    theme_default() +
    facet_wrap(~ nt, scale='fixed', ncol=1) +
    xlab('i') +
    ylab(expression(sqrt(F))) +
    scale_color_continuous(
        guide = guide_legend(title = 'H'),
        low='#DCBCBC', high='#8F2727'
        )

fp_out = '~/projects/influenza/img/feature_importance/NA_fi_entropy_1.pdf'
ggsave(fp_out, r, height=15, width=8, unit='cm')


# x: position, y: entropy, color: feature importance (threshold)


s <-
ggplot(df, aes(x=pos, y=h)) +
    geom_point(size=0.3, color='#DCBCBC') +
    geom_point(data=df[df$fi > 0.005,], color='#8F2727', size=0.3) +
    theme_default() +
    xlab('i') +
    ylab('H')

fp_out = '~/projects/influenza/img/feature_importance/NA_fi_entropy_2.pdf'
ggsave(fp_out, s, height=5, width=5, unit='cm')


# x: entropy, y: feature importance


t <-
ggplot(df, aes(x=h, y=sqrt(fi))) +
    geom_point(size=0.3, alpha=0.5) +
    theme_default() +
    xlab('H') +
    ylab(expression(sqrt(F)))

fp_out = '~/projects/influenza/img/feature_importance/NA_fi_entropy_3.pdf'
ggsave(fp_out, t, height=5, width=5, unit='cm')


# x: position, y: feature importance, facet: nt


u <-
ggplot(df_f, aes(x=pos, y=sqrt(fi))) +
    geom_point(size=0.2) +
    theme_default() +
    facet_wrap(~ nt, scale='fixed', ncol=2) +
    xlab('i') +
    ylab(expression(sqrt(F)))

fp_out = '~/projects/influenza/img/feature_importance/NA_fi_entropy_4.pdf'
ggsave(fp_out, u, height=16, width=5, unit='cm')


# There is a position w/ really low entropy. What is going on there?
subdf <- df[df$fi > 0.005,]
#       pos           h nt          fi
# ...
# 2382  595 1.813309122  C 0.007918206
# 2836  708 0.844838785  G 0.016908831
# 2878  719 0.004495383  C 0.009828745 -
# 2880  719 0.004495383  G 0.009611707 -
# 3354  838 1.555032482  C 0.012731394
# 3433  858 0.953769945  A 0.008158394
# 3436  858 0.953769945  G 0.009635674
# ...


# enter mismatch


df <- inner_join(df, df_m, by='pos')


# x: position, y: feature importance, facet: nt


v <-
ggplot(
        df[df$fi > 0.0002,],
        aes(x=pos, y=sqrt(fi), color=as.factor(mismatch))
        ) +
    geom_point(size=0.2) +
    theme_default() +
    facet_wrap(~ nt, scale='fixed', ncol=1) +
    xlab('i') +
    ylab(expression(sqrt(F))) +
    scale_color_manual(
        values=c('#DCBCBC', '#8F2727'),
        guide = guide_legend(title = 'modified')
        )


fp_out = '~/projects/influenza/img/feature_importance/NA_fi_deopt_DL.pdf'
ggsave(fp_out, v, height=16, width=7, unit='cm')

'''



