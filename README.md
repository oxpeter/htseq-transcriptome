# htseq-transcriptome
htseq script to count reads multiply mapping to a transcriptome assembly

Because a transcriptome assembly includes collections of isoforms
from the same gene, any reads that would have mapped uniquely to
a genome assembly - where there are multiple isoforms - would in
a transcriptome assembly map multiply - once to each isoform.

Because these reads are discarded in the classic implementation of
HTSeq-count, a different implementation needs to be developed to account
for the legitimately multiply-mapped reads in a transcriptome
assembly.

HTSeq-transcriptome provides this functionality. It is built with the
same framework as HTSeq-count, and the same parameters are available
(or will be added soon!). It therefore should be very easy to use for
anyone who is familiar with using HTSeq-count.

The difference is that HTSeq-transcriptome has a unique way of
assessing read bundles. Rather than discard all reads that are
multiply-mapped, it first looks to see if the multiply-mapped loci all
have the same gene_id. In this way, all reads mapping to different
isoforms of the same gene will not get discarded, but are counted as
a SINGLE count toward the gene's expression.