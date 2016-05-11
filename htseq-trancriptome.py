#!/usr/bin/env python

import argparse
import collections
import datetime
import itertools
import logging
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

def ungapped_se_counter(sam_reader, feature_array):
    """
    classic alignment counter for ungapped single end reads
    """
    counts = collections.Counter( )
    for almnt in sam_reader:
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

def ungapped_pe_counter(sam_reader, feature_array):
    counts = collections.Counter( )
    pair_iterator = hts.pair_SAM_alignments( sam_reader, bundle=True )
    # bundle puts all multiply-mapped pairs together.

    t0 = datetime.datetime.now()
    for ic, bundle in enumerate(pair_iterator):

        # report progress (to prove that it is still alive):
        if ic % 1000000 == 0:
            t1 = datetime.datetime.now()
            print "\r%d read bundles counted in %s\r" % (ic, t1-t0)
            sys.stdout.flush()

        if bundle == []: # first bundle for some reason is always an empty list
            continue

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


parser = argparse.ArgumentParser(description=
        """An HTSeq script to count paired-end reads that multiply map, for making
        inferences about gene counts in a transcriptome assembly""")
### input options ###
parser.add_argument("-o", "--output", type=str, default='htseq_transcriptome',
                    help="""specify the filename base to save results to [default is
                    'htseq_transcriptome_counts'""")
parser.add_argument("-d", "--directory", type=str, default=os.getcwd(),
                    help="""specify the directory to save results to [default is current
                    working directory""")


# data file options:
parser.add_argument('alignment_file', nargs=1, type=str,
                    help="""The <alignment_file> contains the aligned reads in the SAM
                    format. (Note that the SAMtools contain Perl scripts to convert most
                    alignment formats to SAM.)""")
parser.add_argument('gtf_file', nargs=1, type=str,
                    help="""The <gtf_file> contains the features in the GFF format.""")

parser.add_argument("-f", "--format", type=str, choices=['sam', 'bam'],
                    help="""Format of the input data. Possible values are sam
                    (for text SAM files) or bam (for binary BAM files). Default is sam""")
parser.add_argument("-r", "--order", type=str,
                    help="""
                    For paired-end data, the alignment have to be sorted either by
                    read name or by alignment position. If your data is not sorted, use
                    the samtools sort function of samtools to sort it. Use this option,
                    with name or pos for <order> to indicate how the input data has been
                    sorted. The default is name.

                    If name is indicated, htseq-count expects all the alignments for the
                    reads of a given read pair to appear in adjacent records in the input
                    data. For pos, this is not expected; rather, read alignments whose
                    mate alignment have not yet been seen are kept in a buffer in memory
                    until the mate is found. While, strictly speaking, the latter will
                    also work with unsorted data, sorting ensures that most alignment
                    mates appear close to each other in the data and hence the buffer
                    is much less likely to overflow.""")
parser.add_argument("-s", "--stranded", choices=['yes', 'no', 'reverse'], type=str,
                    default="yes",
                    help="""
                    whether the data is from a strand-specific assay (default: yes)

                    For stranded=no, a read is considered overlapping with a feature
                    regardless of whether it is mapped to the same or the opposite strand
                    as the feature. For stranded=yes and single-end reads, the read has to
                    be mapped to the same strand as the feature. For paired-end reads,
                    the first read has to be on the same strand and the second read on
                    the opposite strand. For stranded=reverse, these rules are reversed.
                    """)
parser.add_argument("-t", "--type", type=str, default='exon',
                    help="""feature type (3rd column in GFF file) to be used, all features
                    of other type are ignored (default, suitable for RNA-Seq analysis
                    using an Ensembl GTF file: exon)""")
parser.add_argument("-i", "--idattr", type=str, default='gene_id',
                    help="""GFF attribute to be used as feature ID. Several GFF lines with
                    the same feature ID will be considered as parts of the same feature.
                    The feature ID is used to identity the counts in the output table.
                    The default, suitable for RNA-Seq analysis using an Ensembl GTF file,
                    is gene_id.""")
parser.add_argument("-y", "--read_type", choices=['single_end', 'paired_end'], type=str,
                    default="paired_end",
                    help="""
                    whether the data is from a single-end read library, or a paired-end
                    library. Default is paired-end.
                    """)
args = parser.parse_args()



# define some test files:
samfile = '/home/antqueen/booster/PRO_Odontomachus/trinity_denovo_normalized_camponotus/Star/Cplan_Q2_16Aligned.out.sam'
gtffile = '/home/antqueen/genomics/experiments/analyses/PRO20160405_camponotus/trinity_denovo_normalized_camponotus/Transdecoder_ss/merge_genesets/Cpla_td_gff.Apr21_11.15.families.gtf'

# create gtf iterator
print "\nReading gtf file %s..." % (args.gtf_file[0]),
gtf = hts.GFF_Reader(args.gtf_file[0])
print " done."

# create genomic array and populate with exon features (transcripts and genes)
print "Populating genomic array with GTF features...",
sys.stdout.flush()

if args.stranded == 'yes':
    feature_array = hts.GenomicArrayOfSets( "auto", stranded=True)
elif args.stranded == 'no':
    feature_array = hts.GenomicArrayOfSets( "auto", stranded=False)

for feature in gtf:
    if feature.type == args.type:
        feature_array[ feature.iv ] += feature.name

print "done.\n\n"

# create Reader class for samfile:
if args.format == 'sam':
    alnmt_file = hts.SAM_Reader(args.alignment_file[0])
else:
    alnmt_file = hts.BAM_Reader(args.alignment_file[0])


# count reads:
print "Counting reads..."

if args.read_type == 'single_end':
    counts = ungapped_se_counter(alnmt_file, feature_array)

    print "\nSample output for ungapped SE counts:"
    countlist = sorted(counts.items())
    for g, c in countlist[-10:]:
        print "%-10s %d" % (g, c)
else:
    counts = ungapped_pe_counter(alnmt_file, feature_array)

    print "\nSample output for ungapped PE counts:"
    countlist = sorted(counts.items())
    for g, c in countlist[-10:]:
        print "%-10s %d" % (g, c)


# output counts to file:
print "Saving gene counts to file...",
sys.stdout.flush()
outfile = args.output + "_counts.txt"
handle = open(outfile, 'w')
for g, c in countlist:
    handle.write("%-20s %d\n" % (g, c))
handle.close()
print "done."

print "All done."

