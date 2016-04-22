#!/usr/bin/env python

import argparse
import collections
import itertools
import os

import HTSeq as hts


def assess_bundle(bundle, features):
    counts = collections.Counter()
    for almnt in bundle:
        if not almnt.aligned:
            count[ "_unmapped" ] += 1
            continue
        gene_ids = set()
        for iv, val in features[ almnt.iv ].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            counts[ list(gene_ids)[0] ] += 1
        elif len(gene_ids) == 0:
            counts[ "_no_feature" ] += 1
        else:
            counts[ "_ambiguous" ] += 1
    return counts

# define some test files:
samfile = '/home/antqueen/booster/PRO_Odontomachus/trinity_denovo_normalized_camponotus/Star/Cplan_Q2_16Aligned.out.sam'
gtffile = '/home/antqueen/genomics/experiments/analyses/PRO20160405_camponotus/trinity_denovo_normalized_camponotus/Transdecoder_ss/merge_genesets/Cpla_td_gff.Apr21_11.15.families.gtf'

# create gtf iterator
gtf = hts.GFF_Reader(gtffile)

# create genomic array and populate with exon features (transcripts and genes)
exons = hts.GenomicArrayOfSets( "auto", stranded=True)
for feature in gtf:
    if feature.type == "exon":
        exons[ feature.iv ] += feature.name

# create Reader class for samfile:
sam_file = hts.SAM_Reader(samfile)

# classic alignment counter for ungapped single end reads:
counts = collections.Counter( )
for almnt in itertools.islice( sam_file, 100):
    if not almnt.aligned:
        count[ "_unmapped" ] += 1
        continue
    gene_ids = set()
    for iv, val in exons[ almnt.iv ].steps():
        gene_ids |= val
    if len(gene_ids) == 1:
        gene_id = list(gene_ids)[0]
        counts[ gene_id ] += 1
    elif len(gene_ids) == 0:
        counts[ "_no_feature" ] += 1
    else:
        counts[ "_ambiguous" ] += 1

for g, c in counts.items():
    print "%-10s %d" % (g, c)


read_pair_bundles = hts.pair_SAM_alignments( sam_file , bundle=True)
