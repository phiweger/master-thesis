library(dplyr)
library(ggplot2)
library('RColorBrewer')
source('.../theme_minimal.R')


df = read.table('pca_trinucleotide.tsv', header=T, sep='\t')
names(df) = c('PC1', 'PC2', 'Host', 'Subtype')
df2 <- dplyr::select(df, -Host, -Subtype)


# TODO: add axis


# host
p <- ggplot(df, aes(x=PC1, y=PC2)) + 
    geom_point(data=df2, color='grey70', size=.2) +
    geom_point(aes(color=Host), size=.2, alpha=.3) +
    coord_fixed(ratio=1) +
    facet_wrap(~Host, ncol=3) +
    theme_minimal() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
    scale_y_continuous(breaks=c(-20, 0, 20)) +
    scale_x_continuous(breaks=c(-20, 0, 20)) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size=2)))
    # http://stackoverflow.com/questions/5290003/how-to-set-legend-alpha-with-ggplot2

    # scale_colour_brewer(palette = "Set1")
ggsave('.../pca_5b08d.png', 
    p, width=18, height=6, units='cm')
ggsave('.../pca_5b08d.pdf', 
    p, width=18, height=6, units='cm')


# subtype
# invert host and subtype
q <- ggplot(df, aes(x=PC1, y=PC2, color=Subtype)) + 
    geom_point(data=df2, color='grey70', size=.2) +
    geom_point(size=.2, alpha=.3) +
    coord_fixed(ratio=1) +
    facet_wrap(~Host) +
    theme_minimal() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
    scale_y_continuous(breaks=c(-20, 20)) +
    scale_x_continuous(breaks=c(-20, 20)) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size=2)))

# color host, facet subtype
ggsave('.../pca_0aaeb.png', 
    q, width=25, height=25, units='cm')
ggsave('.../pca_0aaeb.pdf', 
    q, width=25, height=25, units='cm')
# color subtype, facet host
ggsave('.../pca_0dcec.png', 
    q, width=25, height=15, units='cm')
ggsave('.../pca_0dcec.pdf', 
    q, width=25, height=15, units='cm')


r <- ggplot(df, aes(x=PC1, y=PC2, color=Subtype)) + 
    #geom_point(data=df2, color='grey70', size=.2) +
    geom_point(size=.2, alpha=.3) +
    coord_fixed(ratio=1) +
    # facet_wrap(~Subtype) +
    theme_minimal() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
    scale_y_continuous(breaks=c(-20, 20)) +
    scale_x_continuous(breaks=c(-20, 20)) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size=2)))

ggsave('.../pca_4566a.png', 
    r, width=15, height=15, units='cm')
ggsave('.../pca_4566a.pdf', 
    r, width=15, height=15, units='cm')

# > t(table(dplyr::select(df, Host, Subtype)))
#        Host
# Subtype Avian Human Swine
#   H10N1    32     0     0
#   H10N2    14     0     0
#   H10N3    84     0     0
#   H10N4    29     0     0
#   H10N5    19     0     1
#   H10N6    52     0     0
#   H10N7   388     0     0
#   H10N8    63     4     0
#   H10N9    16     0     0
#   H11N1    24     0     0
#   H11N2    74     0     0
#   H11N3    23     0     0
#   H11N8     9     0     0
#   H11N9   294     0     0
#   H12N5   106     0     0
#   H13N6    26     0     0
#   H13N8    12     0     0
#   H16N3    97     0     0
#   H1N1    336  6948  1225
#   H1N2     33    21   634
#   H1N3     29     0     0
#   H1N6     36     0     0
#   H1N8     14     0     0
#   H1N9     45     0     0
#   H2N1     32     0     0
#   H2N2     47    73     0
#   H2N3    202     0     2
#   H2N5     17     0     0
#   H2N7     22     0     0
#   H2N8     13     0     0
#   H2N9     41     0     0
#   H3N1     28     0    15
#   H3N2    213  6551   914
#   H3N3     14     0     2
#   H3N5     18     0     0
#   H3N6    130     0     1
#   H3N8    817     0     2
#   H3N9     11     0     0
#   H4N1     10     0     1
#   H4N2    122     0     0
#   H4N3     15     0     0
#   H4N4     16     0     0
#   H4N6    878     0     1
#   H4N8    167     0     1
#   H4N9     20     0     0
#   H5N1    950    60    20
#   H5N2    442     0     0
#   H5N3     43     0     0
#   H5N5     37     0     0
#   H5N6     84     3     2
#   H5N8    152     0     0
#   H5N9     11     0     0
#   H6N1    252     0     0
#   H6N2    373     0     0
#   H6N3     13     0     0
#   H6N5     51     0     0
#   H6N6    169     0     1
#   H6N8    141     0     0
#   H7N1     73     0     0
#   H7N2    123     0     1
#   H7N3    398     2     0
#   H7N4     18     0     0
#   H7N6     26     0     0
#   H7N7    133     2     0
#   H7N8     13     0     0
#   H7N9    409    72     0
#   H8N4     91     0     0
#   H9N1     16     0     0
#   H9N2    843     4    18
#   H9N5     11     0     0


dfc = read.table(
    '.../count_trinucleotide_annotated.tsv', 
    header=T, sep='\t')
dfc2 <- dplyr::select(dfc, -host)

# wide2long
# http://stackoverflow.com/questions/2185252/reshaping-data-frame-from-wide-to-long-format
dfc_long <- reshape2::melt(dfc, id=c('host', 'subtype'), 
    variable_name='triplet')

s <- ggplot(dfc_long, aes(x=variable, y=value)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~host, scale='fixed', nrow=3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('codon') + ylab('count')

ggsave('.../boxplot_016c0.png', 
    s, width=30, height=20, units='cm')
ggsave('.../boxplot_016c0.pdf', 
    s, width=30, height=20, units='cm')



