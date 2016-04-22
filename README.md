# htseq-transcriptome
htseq script to count reads multiply mapping to a transcriptome assembly

Because a transcriptome assembly includes collections of isoforms
from the same gene, any reads that would have mapped uniquely to 
a genome assembly - where there are multiple isoforms - would in 
a transcriptome assembly map multiply - once to each isoform.

Because these reads are discarded in the classic implementation of 
HTSeq, a different implementation needs to be developed to account
for the legitimately multiply-mapped reads in a transcriptome 
alignment.
