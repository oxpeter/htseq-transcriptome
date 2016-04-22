#!/usr/bin/env python

import argparse
import collections
import itertools
import os
import sys

import HTSeq as hts


def assess_bundle(bundle, features):
    """
    Takes a bundle of (potentially multiply mapped) paired-end read pairs, and looks to
    see how many features map to each aligned pair.
    """
    counts = collections.Counter()
    for (p1,p2) in bundle:
        if not p1 or not p2: # ie, one of the mate pairs is missing
            counts[ "_unmapped" ] += 1
            continue
        elif not (p1.aligned and p2.aligned): # mate pairs present, but unaligned
            counts[ "_unmapped" ] += 1
            continue

        # collect all genes that map to this alignment
        gene_ids = set()
        for iv, val in features[ p1.iv ].steps():
            gene_ids |= val
        for iv, val in features[ p2.iv ].steps():
            gene_ids |= val

        # evaluate:
        if len(gene_ids) == 1:
            counts[ list(gene_ids)[0] ] += 1
        elif len(gene_ids) == 0:  # TODO: test that mate pairs are matching features!
            counts[ "_no_feature" ] += 1
        else:
            counts[ "_ambiguous" ] += 1
    return counts

def ungapped_se_counter_demo(sam_reader, feature_array):
    """
    classic alignment counter for ungapped single end reads
    """
    counts = collections.Counter( )
    for almnt in itertools.islice( sam_reader, 100):
        if not almnt.aligned:
            count[ "_unmapped" ] += 1
            continue
        gene_ids = set()
        for iv, val in feature_array[ almnt.iv ].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[ gene_id ] += 1
        elif len(gene_ids) == 0:
            counts[ "_no_feature" ] += 1
        else:
            counts[ "_ambiguous" ] += 1

    return counts

def ungapped_pe_counter_demo(sam_reader, feature_array):
    counts = collections.Counter( )
    pair_iterator = hts.pair_SAM_alignments( sam_reader, bundle=True )
    # bundle puts all multiply-mapped pairs together.
    for bundle in itertools.islice( pair_iterator, 100 ):
        bcounts = assess_bundle(bundle, feature_array)

        """
        To evaluate the multiply mapped bundles, each pair in a bundle must still ALWAYS
        and ONLY map to a single feature. Thus, every aligned pair has come from the same
        feature (gene), and this bundle counts as evidence of one read for this gene.

        If any of the read pairs maps to a different gene, or no gene, or multiple genes,
        then the bundle is considered ambiguous.

        If all pairs in a bundle map as _no_feature, _unmapped or _ambiguous, then the
        bundle counts as one count towards this feature type. (ie, it is passed on to
        the final counter to increment by 1).
        """

        if len(bcounts) > 1: # ie, is a multiply mapped feature with multiple gene mappings
            counts[ "_ambiguous" ] += 1
            continue
        elif len(bcounts) == 0:  # uh oh! There is an error somewhere.
            print "#" * 40
            print "Error! bundle was not assigned any status"
            print "Contents of bundle:"
            print bundle
            continue
        else:
            counts[ bcounts.keys()[0] ] += 1

    return counts

# define some test files:
samfile = '/home/antqueen/booster/PRO_Odontomachus/trinity_denovo_normalized_camponotus/Star/Cplan_Q2_16Aligned.out.sam'
gtffile = '/home/antqueen/genomics/experiments/analyses/PRO20160405_camponotus/trinity_denovo_normalized_camponotus/Transdecoder_ss/merge_genesets/Cpla_td_gff.Apr21_11.15.families.gtf'

# create gtf iterator
print "Reading gtf file...",
gtf = hts.GFF_Reader(gtffile)
print " done."

# create genomic array and populate with exon features (transcripts and genes)
print "Populating genomic array with GTF features...",
sys.stdout.flush()
exons = hts.GenomicArrayOfSets( "auto", stranded=True)
for feature in gtf:
    if feature.type == "exon":
        exons[ feature.iv ] += feature.name
print "done."

# create Reader class for samfile:
sam_file = hts.SAM_Reader(samfile)

# get counts assuming ungapped SE reads :
counts = ungapped_se_counter_demo(sam_file, exons)
print "\nUngapped SE counts:"
for g, c in sorted(counts.items()):
    print "%-10s %d" % (g, c)

# get counts assuming ungapped PE reads :
counts = ungapped_pe_counter_demo(sam_file, exons)

print "\nUngapped PE counts:"
for g, c in sorted(counts.items()):
    print "%-10s %d" % (g, c)

