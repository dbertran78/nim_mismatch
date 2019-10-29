

Fast computation of the number of mismatches using the cigar string operation and NM tag
```nim compile --run get_mismatch_statistics.nim --bam=DATA/C1000.bam  --ref=DATA/human_g1k_v37_M.fasta```

Compute the type of mismatches and create a sim file with th '=' and 'X' cigar string operation
```nim compile --run update_bam.nim --bam=DATA/C1000.bam  --new-sam=gg.sam --ref=DATA/human_g1k_v37_M.fasta```