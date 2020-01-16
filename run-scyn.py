import argparse
import scyn
import os

parser = argparse.ArgumentParser(description="SCYN: Single cell CNV profiling method usingdynamic programming efficiently andeffectively", add_help=False, usage="python %(prog)s [-h] [options] -i input_bams_dir", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")
required.add_argument("-i", "--indir", type=str, required=True,
                    help="<str> the input bams directory", metavar="\b")
optional.add_argument("-o", "--outdir", type=str, default="./",
                    help="<str> the output directory", metavar="\b")
optional.add_argument("--seq", type=str,
                    default="single-end", help="<str> the reads type: single-end or paired-end.", metavar="\b")
optional.add_argument("--bin_len", type=int,
                    default="500", help="<int> the bin length, default is 500K.", metavar="\b")
optional.add_argument("--ref", type=str,
                      default="hg19", help="<str> the reference genome version: hg19 or hg38.", metavar="\b")
optional.add_argument("--reg", type=str,
                      default="*.bam", help="<str> the regular expression to match all BAM files in your input directory. For example, \".bam\" will match all BAM files ended with '.bam'.", metavar="\b")
optional.add_argument("--mapq", type=int,
                      default="40", help="<int> the mapping quality cutoff when calculating the reads coverage.", metavar="\b")
optional.add_argument("--K", type=int,
                      default="10", help="<int> the changepoints number for each chromosome.", metavar="\b")
optional.add_argument("--verbose", type=int,
                      default="1", help="<int> If > 0, print log messages.", metavar="\b")
optional.add_argument("-h", "--help", action="help")
args = vars(parser.parse_args())

# create SCYN object
scyn_operator = scyn.SCYN(seq=args.seq, bin_len=args.bin_len,
                          ref=args.ref, reg=args.reg, mapq=args.mapq, K=args.K, verbose=args.verbose)

# call cnv
# bam_dir is the input bam directory and output_dir is the output directory
scyn_operator.call(args.indir, args.outdir)

# store cnv matrix to a csv file
scyn_operator.cnv.to_csv(os.path.join(args.outdir, 'cnv.csv'))
scyn_operator.cnv.T.to_csv(os.path.join(args.outdir, 'cnv_T.csv'))
scyn_operator.cnv.T.to_csv(os.path.join(args.outdir, 'cnv_T.csv'))
scyn_operator.meta.T.to_csv(os.path.join(args.outdir, 'meta.csv'))
