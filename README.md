
Notice that the reference genome used for the mapping is different for the one [here](ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz)
Some region in Y chromosome seems to have been masked (only affect a single read)
The mitochondrial chromosome have a different name. To correct:
```sed 's/>MT/>M/' DATA/human_g1k_v37.fasta > DATA/human_g1k_v37_M.fasta```

Fast computation of the number of mismatches using the cigar string operation and NM tag
```nim compile --run get_mismatch_statistics.nim --bam=DATA/C1000.bam  --ref=DATA/human_g1k_v37_M.fasta```

Compute the type of mismatches and create a sam file with th '=' and 'X' cigar string operation
```nim compile --run update_bam.nim --bam=DATA/C1000.bam  --new-sam=gg.sam --ref=DATA/human_g1k_v37_M.fasta```