from __future__ import division
from optparse import OptionParser
import sys, os, re
import ntpath
import math
import csv

import pysam
import pybedtools
import numpy as np
from scipy import stats


def get_depth(options):
    """Get the number of mapped reads in both t and n bam files in a genomic region (e.g. X:1000000-2000000)
       Then get the number of mapped reads within the CNV region
    """
    tumour_name = get_sample(options.tumour)
    normal_name = get_sample(options.normal)

    print("Window size set to: %s" % options.window)
    print("Sampling within region %s" % options.region)
    chromosomes = ["2L", "2R", "3L", "3R", "4", "X", "Y"]

    # Count reads in both bam files
    t_reads_by_chrom, tumour_mapped = count_reads(options.tumour, chromosomes)
    n_reads_by_chrom, normal_mapped = count_reads(options.normal, chromosomes)
    print("Total reads in %s: %s" % (tumour_name, tumour_mapped))
    print("Total reads in %s: %s" % (normal_name, normal_mapped))

    # Calculate ratio to correct normal read counts by
    mapped_ratio = round(tumour_mapped / normal_mapped, 3)
    print("Using read count ratio [%s] to normalise reads between bams" % (mapped_ratio))

    chrom, start, end = extract_loci(options.region)

    with open('window_size.bed', 'w') as region_out: region_out.write('\t'.join([chrom, str(0), str(options.window)]))
    a = pybedtools.BedTool('window_size.bed')

    with open('incl.bed', 'w') as out: out.write('\t'.join([chrom, str(start), str(end)]))

    norms = []
    tums = []

    if options.debug:
        t_out = tumour_name + '.bed'
        if os.path.exists(t_out): os.remove(t_out)

    for i in range(1, int(options.sample)):
        c, s, e = str(a.shuffle(genome={'3R':(0,32079331)}, incl='incl.bed')).split()
        s = int(s)
        e = int(e)
        if i % 10 == 0: print("Sampled %s regions" % (i))
        tum_count = region_depth(options.tumour, c, s, e, options)
        norm_count = region_depth(options.normal, c, s, e, options)
        tums.append(tum_count)
        norms.append(norm_count)
        if options.debug:
            t_out = tumour_name + '.bed'
            info = "T:" + str(tum_count) + "/N:" + str(norm_count)
            with open(t_out, mode='a') as counts_file:
                counts_writer = csv.writer(counts_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                counts_writer.writerow([c, s, e, info])


    norm_adj = []
    for n in norms:
        norm_adj.append(int(n * mapped_ratio))

    print("average read count per window:")
    print(" o tumour (%s): %s" % (tumour_name, round(np.average(tums), 2)))
    print(" o normal (%s): %s" % (normal_name, round(np.average(norms), 2)))
    print(" p adjusted normal: %s" % (round(np.average(norm_adj), 2)))

    fc, log2FC = calc_ratio(tums, norm_adj)
    print("Average ratio: %s [log2: %s]" % (fc, log2FC))

    t, p = stats.ttest_ind(tums, norm_adj)

    print(p)

    sig_val = is_sig(p, fc, log2FC)

    info = [tumour_name, round(np.average(tums))]
    # info.append(tumour_name, round(np.average(tums)))

    with open('stats.csv', mode='a') as stats_file:

        info_writer = csv.writer(stats_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # info_writer.writerow(["sample", 'tissue', 'total mapped reads', 'region', 'n samples', 'average read count in window', 'fc', 'log2FC' ])
        info_writer.writerow([tumour_name, 'tumour', tumour_mapped, options.region, options.sample, round(np.average(tums), 2), fc, log2FC, p, sig_val])
        info_writer.writerow([normal_name, 'normal', normal_mapped, options.region, options.sample, round(np.average(norm_adj), 2), fc, log2FC, p, sig_val])

    return True


def get_sample(b):
    return ntpath.basename(b).split(".")[0]


def is_sig(p, fc, log2FC):
    sig_val = "-"
    if p > 0.01:
        print(" -> No significant difference (p: %s) found: fc ratio: %s [log2: %s]" % (round(p, 5), fc, log2FC))
        return sig_val
    sig_val = "*"
    if p <= 0.001: sig_val = "**"
    if p <= 0.0001: sig_val = "***"

    if log2FC <= 0:
        print(" -> *** Sig (p: %s) contamination found in tumour: fc ratio: %s [log2: %s]" % (round(p, 5), fc, log2FC))
        return sig_val
    print(" -> *** (p: %s) contamination found in normal: fc ratio: %s [log2: %s]" % (round(p, 5), fc, log2FC))
    return sig_val


def count_reads(bamfile, chromosomes):
    """Count the total number of mapped reads in a BAM file, filtering
    for chromosomes in `chromosomes`
    """

    lines = pysam.idxstats(bamfile)

    if type(lines) is str:
        lines = lines.strip().split('\n')

    total_mapped = 0
    total_mapped_chrom = {}

    for line in lines:
        chrom, len, mapped, unmapped = line.split('\t')
        if chrom in chromosomes:
            total_mapped_chrom[chrom] = int(mapped)
            total_mapped += int(mapped)

    return total_mapped_chrom, total_mapped


def calc_ratio(t, n):
    r = []
    for i,j in zip(t,n):
        i = i + 0.01
        j = j + 0.01
        r.append(i/j)

    fc = round(np.average(r), 2)
    log2fc = round(math.log2(np.average(r)), 3)

    return fc, log2fc


def region_depth(bamfile, chrom, start, end, options):
    """Count the total number of mapped reads in a genomic region"""

    samfile = pysam.Samfile(bamfile, "rb")
    count = 0

    for read in samfile.fetch(chrom, start, end):
        if read.is_unmapped:
            continue
        elif read.mapq < 3:
            continue

        count += 1

    return count


def extract_loci(region):
    try:
        chrom, start, end = re.split('[:\\-]', region)
    except ValueError:
        sys.exit("[!] Must provide a valid region [ chrom:start:end ].Exiting.")

    return chrom, int(start), int(end)


def get_args():
    parser = OptionParser()

    parser.add_option("-t", "--tumour", dest="tumour", action="store", help="Bam file for tumour sample")
    parser.add_option("-n", "--normal", dest="normal", action="store", help="Bam file for normal sample")
    parser.add_option("-s", "--sample", dest="sample", action="store", help="Number of times to sample [Default: 10]")
    parser.add_option("-r", "--region", dest="region", action="store", help="Region to compare. E.g. X:100000-200000")
    parser.add_option("-w", "--window", dest="window", action="store", help="Window in bps to sample. [Default: 10000] ")
    parser.add_option("-o", "--out_file", dest="out_file", action="store", help="File to write annotated vars to")
    parser.add_option("--debug", dest="debug", action="store_true", help="Debug by writing bed files for sampled regions")

    parser.set_defaults(window=10000,
                        sample=10,
                        tumour='/Volumes/perso/Analysis/Temp/L342-01.tagged.filt.SC.RG.bam',
                        normal='/Volumes/perso/Analysis/Temp/L342-02.tagged.filt.SC.RG.bam')

    options, args = parser.parse_args()

    if not options.tumour or not options.normal:
        parser.print_help()
        print()
        sys.exit("[!] Must provide a --tumour and --normal bam files. Exiting.")
    elif not options.region:
        parser.print_help()
        print()
        sys.exit("[!] Must specify a genomic --region. Exiting.")
    else:
        return options, args


def main():
    options, args = get_args()

    print("deepy")
    get_depth(options)


if __name__ == "__main__":
    sys.exit(main())
