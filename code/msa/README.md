## Multiple alignment

For virus genomes, the "linsi" parameter set works well. This is because even if viruses have widely divergent genomes, this configuration will retrieve common motivs nonetheless.

~~~
mafft --localpair --maxiterate 1000 input_file > output_file
# or
linsi input_file > output_file

# e.g. MSA .. multisequence alignment for hemagglutinin protein
grep '_HA_' -A1 iva_cds_singleline.fa | grep -v "\-\-" - | head -n 100 | linsi - > ../tmp/msa_ha.fa
# output: aligned sequence in fasta format
# >gb|AB049159:1-1683_cds:BAB39511_13383279_HA_4_A/parakeet/Chiba/1/97_H9N2_Avian
# atggaaacaatatc------------actactaactatactactag------------ta
# gtaacagcaagcaatgcagataaaatctgcatcggccaccagtcaacaaactccacagaa
# [...]
~~~

>Iterative refinement methods:
* FFT-NS-i (Slow; iterative refinement method)
* E-INS-i (Very slow; recommended for <200 sequences with multiple conserved domains and long gaps) Help  Updated (2015/Jun)
* L-INS-i (Very slow; recommended for <200 sequences with one conserved domain and long gaps) Help
* G-INS-i (Very slow; recommended for <200 sequences with global homology) Help
* Q-INS-i (Extremely slow; secondary structure of RNA is considered; recommended for a global alignment of highly divergent ncRNAs with <200 sequences × <1,000 nucleotides; the number of iterative cycles is restricted to two, 2016/May) 

[help](mafft.cbrc.jp/alignment/software/algorithms/algorithms.html)


# Vis

[MSAViewer](http://msa.biojs.net/) wrapped in [an R package](https://github.com/bene200/msaR)

~~~ R
install.packages('devtools')
install.packages('htmlwidgets')
devtools::install_github('bene200/msaR')
library(msaR)
msaR('/Users/pi/projects/influenza/data/tmp/msa_ha')
~~~


## How to "glue" sequences together (if need be)?

Markus puts 100 "N" characters between segments, done.