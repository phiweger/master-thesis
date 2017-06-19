## Mafft

- how to [append](http://mafft.cbrc.jp/alignment/software/addsequences.html) sequence to existing alignment

```
# codon based alignment
python ~/commonplace/compbio/protocols/alignment/msa/codon_msa.py subset.fa
# creates codon.fa and protein.fa

# or if codon alignment not necessary: mafft default
# TODO: make this happen in temporary file or create run-specific folder
# "analysis_NA_2017-03-13"
mafft --thread -1 subset_NA.fa > subset_NA.mafft.fa

# append sequence
mafft --add ~/data/influenza/external/dusan_deoptimized_influenza/H1N1_NA_wildtype.fa --reorder subset_NA.mafft.fa > subset_NA_wt.mafft.fa

grep ">" subset_NA.mafft.fa | wc -l
# 9333
grep ">" subset_NA_wt.mafft.fa | wc -l
# 9334
```
