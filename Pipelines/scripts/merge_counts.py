# this script merges the results from featureCounts or Salmon in to a single gene X cell counts matrix

"""
Merge count tables for genes for multiple samples into a gene X single cell/sample counts matrix

Input files as a list of comma-separated file paths

Default output is stdout, redirect using the --output flag
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
import re
import logging

# set up the command line parser
parser = argparse.ArgumentParser(description='Merge gene count matrices together')
parser.add_argument("--output", dest="output", type=str,
                    default=sys.stdout,
                    help="Output file path")

parser.add_argument("--input-directory", dest="input_dir", type=str,
                    help="The directory containing the input files")

parser.add_argument("--file-regex", dest="file_regex", type=str,
                    help="Regular expression to match to input files. Accepts "
                    "normal regular expression syntax")

parser.add_argument("--input-format", dest="in_format", type=str,
                    choices=["salmon", "feature"],
                    default="feature",
                    help="input file format, should be either feature or salmon")

parser.add_argument("--log", dest="logfile", type=str,
                    help="destination file for script logging",
                    default=sys.stderr)

args = parser.parse_args()

# setup the logger
logging.basicConfig(level=logging.INFO,
                        filename=args.logfile)

# can't take thousands of input files because posix has a limit to the
# character limit of an argument to a script
# need to pass an input directory and a glob/regex instead
# could be a problem for snakemake to handle?
reg_compile = re.compile(args.file_regex)
if args.in_format == "feature":
    found_files = [ft for ft in os.listdir(args.input_dir) if re.search(reg_compile, ft)]
    input_files = [os.path.join(args.input_dir, fx) for fx in found_files]
elif args.in_format == "salmon":
    found_files = [ft for ft in os.listdir(args.input_dir) if re.search(reg_compile, os.path.join(ft, "quant.sf"))]
    input_files = [os.path.join(args.input_dir, fx, "quant.sf/quant.sf") for fx in found_files]
else:
    raise ValueError("Did not recognize input format.  Should be either 'feature' or 'salmon'")
#input_files = sys.argv[-1].split(",")

logging.info("Found {} count files to merge".format(len(input_files)))

# this WILL take a long time for many thousands of files
# sequentially open and merge files based on the following columns:
# Geneid, Chr, Start, End, Strand, Length
# parse the file name input filename

n_files = 0
for cfile in input_files:
    if args.in_format == "feature":
        fname = cfile.split("/")[-1]
        samp_name = fname.split(".")[:-2]
    elif args.in_format == "salmon":
        fname = cfile.split("/")[-3]
        samp_name = ".".join(fname.split("-")[:-2])
    if n_files == 0:
        count_df = pd.read_table(cfile, sep="\t", header=0, index_col=None, 
                                 comment='#')
        if args.in_format == "salmon":
            count_df = count_df.loc[:, ["Name", "Length", "EffectiveLength", "NumReads"]]
        columns = list(count_df.columns[:-1])
        columns.append(samp_name[0])
        count_df.columns = columns
        n_files += 1
    else:
        _df = pd.read_table(cfile, sep="\t", header=0, index_col=None, 
                                 comment='#')
        if args.in_format == "salmon":
            _df = _df.loc[:, ["Name", "Length", "NumReads"]]
        columns = list(_df.columns[:-1])
        columns.append(samp_name[0])
        _df.columns = columns
        n_files += 1
        # merge
        if args.in_format == "feature":
            count_df = pd.merge(count_df, _df,
                                left_on=["Geneid", "Chr", "Start", "End", "Strand", "Length"],
                                right_on=["Geneid", "Chr", "Start", "End", "Strand", "Length"])
        elif args.in_format == "salmon":
            count_df = pd.merge(count_df, _df,
                                left_on=["Name", "Length"],
                                right_on=["Name", "Length"])
        if not n_files % 100:
            logging.info("Processed {} count files".format(n_files))

logging.info("Writing merged gene counts to file: {}".format(args.output))
count_df.to_csv(args.output, sep="\t",
                index=None)
