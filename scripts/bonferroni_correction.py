#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import re
import logging
from scipy import stats
from argparse import ArgumentParser


def list_files(path_pattern):
    """List files recursively with specific pattern"""
    file_list = []
    for name in glob.glob(path_pattern):
        # base = os.path.basename(name)
        file_list.append(name)
    return file_list


def files_to_dict(results_file):
    """Returns a dictionary with gene name as key and
    results files as values"""
    results = dict()
    for x in results_file:
        gene = os.path.basename(x).split("_")[0]
        if gene not in results.keys():
            results[gene] = []
        results[gene].append(x)
    return results



def lnl(codeml_in):
    """Extracts Ln Likelihood values from output files created by codeml"""
    with open(codeml_in, "r", newline=None) as file_in:
        if os.stat(codeml_in).st_size > 0:  # if file is not empty
            for line in file_in:
                if line.startswith("lnL"):
                    vals = re.findall(r"[-]?\d*\.\d+|\d+", line)
                    for x in vals:
                        if float(x) < 0:
                            log = float(x)
                    break
        else:
            log = None
    return log


def chisqprob(chisq, df):
    value = stats.chi2.sf(chisq, df)
    return value


def lrt(lnl1, lnl0, df=1):
    """Likelihood Ratio Test calculator.
    Degrees of freedom (df, default = 1)"""
    if lnl1 and lnl0:
        val = 2 * (float(lnl1) - (float(lnl0)))
        p_val = chisqprob(val, df)
    else:
        p_val = None
    return p_val


def if_significant(alt_file, null_file, alpha_value=0.05):
    """Returns true when a gene analyzed
    has a pval < alpha according to LRT"""
    ln_alt = lnl(alt_file)
    ln_null = lnl(null_file)
    p_val = lrt(ln_alt, ln_null)
    if p_val < alpha_value:
        return True
    else:
        return False

def bonferroni_alpha(alternative_files):
    """Calculates new level of significance using
    multiple branch testing Bonferroni's correction"""
    n = len(alternative_files)
    if n == 0:
        logging.error("No GWideCodeMl results files found.")
        exit()
    alpha_value = 0.05 / n
    return alpha_value



def parse_arguments():
    parser = ArgumentParser(description="Parameters for multiple branch tetsing with Bonferroni's correction")
    parser.add_argument_group('required arguments')
    parser.add_argument("-model", dest="mode", help="model testing, BM or BS", type=str, default="BS")
    parser.parse_args()
    arguments = parser.parse_args()
    return arguments


# Main steps of the program
def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s :: %(levelname)s :: %(message)s')
    args = parse_arguments()
    path = os.path.abspath(os.getcwd())
    if args.mode not in ["BS", "BM"]:
        logging.error("Multiple branch testing correction can be only applied for multiple branch testing of"
                      "branch (BM) or branch-site (BS) models")
        exit()
    pattern = "/*/*_" + args.mode + "_alt.txt"
    alt_files = list_files(path + pattern)
    logging.info("{} hypotheses will be corrected.".format(len(alt_files)))
    alpha = bonferroni_alpha(alt_files)
    logging.info("New level of significance: {}".format(alpha))

    tests = files_to_dict(alt_files)
    logging.info("Applying multiple testing correction to {} genes".format(len(tests.keys())))
    significants = []
    for k in tests.keys():
        for v in tests[k]:
            alt = v
            null = alt.replace("_alt.txt", "_null.txt")
            if if_significant(alt, null, alpha):
                significants.append(os.path.basename(alt))

    # Corrected hypothesis
    if len(significants) != 0:
        with open("multiple_branch_correction.txt", "w") as outfile:
            outfile.write("\n".join(str(item) for item in significants))

        logging.info("Null hypothesis rejected in {} cases. "
                     "Check multiple_branch_correction.txt file.".format(len(significants)))

    elif len(significants) == 0:
        logging.info("Alternative hypothesis not accepted in any case.")
        exit()






if __name__ == "__main__":
    main()
