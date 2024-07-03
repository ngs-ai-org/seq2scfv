#!/usr/bin/env python3

"""
This script processes scFv annotation data to filter sequences based on 
various criteria such as linker quality, domain order, productivity, and chain length.
"""

# Import necessary libraries
import pandas as pd 
from optparse import OptionParser
import logging

                
# Parse input arguments
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)

parser.add_option("--scFv", dest="hits", default=None, type="string", action="store", help="TSV file containing annotations of all in-frame scFv.")
parser.add_option("--linker-min-id", dest="idpcnt", default="90", type="string", action="store", help="Flag scFvs with a linker score (%) lower than the selected threshold. Default: 90.")
parser.add_option("--mismatches", dest="mm", default="2", type="string", action="store", help="Flag scFvs with a linker presenting more mismatches than the selected threshold. Default: 2.")
parser.add_option("--linker-overhang", dest="lo", default=None, type="string", action="store", help="Linker overhang (in number of aa) allowed relative to the reference linker. Default: None.")
parser.add_option("--chain-length", dest="cl", default=None, type="string", action="store", help="Minimum V-domain length in aa. Default: None.")
parser.add_option("--log", dest="loglevel", type="string", default="info")

(options, args) = parser.parse_args()

input = options.hits
score_pcnt = options.idpcnt
mismatches = options.mm
linkeroverhang = options.lo
chainlen = options.cl
loglevel = options.loglevel

# Set logging level
numeric_level = getattr(logging, loglevel.upper(), None)
if not isinstance(numeric_level, int):    
    raise ValueError('Invalid log level: %s' % loglevel)
logging.basicConfig(level=numeric_level)

logging.debug(options)

# Validate essential options
if input is None:
    raise RuntimeError("Error! no scFv dataframe provided (--scFv)")

def read_scfv_table(input): 
    """
    Reads the scFv annotation data from a TSV file.

    Parameters:
    input (str): Path to the input TSV file.

    Returns:
    pd.DataFrame: DataFrame containing the scFv annotation data.
    """
    annot = pd.read_csv(input, sep = "\t", header = 0, index_col = [0])
    return(annot)

def check_domain_order(VH, VL): 
    """
    Checks if the VH domain appears before the VL domain.

    Parameters:
    VH (str): Position of the VH domain.
    VL (str): Position of the VL domain.

    Returns:
    int: 1 if VH appears before VL, else 0.
    """
    if int(VH) < int(VL) : 
        return 1
    else: 
        return 0

def check_linker_qual(pcnt, mm, overhang, score_pcnt, mismatches, linkeroverhang):
    """
    Checks if the linker meets the quality criteria based on score percentage,
    mismatches, and overhang.

    Parameters:
    pcnt (str): Percentage score of the linker.
    mm (str): Number of mismatches in the linker.
    overhang (str): Number of amino acids overhanging in the linker.
    score_pcnt (str): Threshold score percentage.
    mismatches (str): Allowed mismatches.
    linkeroverhang (str): Allowed linker overhang.

    Returns:
    int: 1 if the linker passes the quality check, else 0.
    """
    if pd.notna(pcnt) and pd.notna(mm) and pd.notna(overhang):
        if int(pcnt) >= int(score_pcnt) and int(mm) <= int(mismatches):
            if linkeroverhang == "undefined":
                return 1
            elif int(overhang) <= int(linkeroverhang):
                return 1
            else:
                    return 0
        else:
                return 0
    else:
            return 0

def check_productive_scfv(VH, VL): 
    """
    Checks if both VH and VL domains are productive.
    """
    if VH == "T" and VL =="T": 
        return 1
    else: 
        return 0

def check_complete_scfv(VH, VL): 
    """
    Checks if both VH and VL domains are complete.
    """
    if VH == "T" and VL == "T": 
        return 1
    else: 
        return 0
    
def check_chain_length(VH, VL, chainlen): 
    """
    Checks if both VH and VL sequences meet the minimum chain length requirement.
    """
    if (len(VH) >= int(chainlen)) and (len(VL) >= int(chainlen)):
        return 1
    else: 
        return 0
        
def main(input, score_pcnt, mismatches, linkeroverhang, chainlen): 
    df = read_scfv_table(input)
    df["VH_IGH"] = df.apply(
        lambda row: check_domain_order(
            row["VH_v_sequence_start"], 
            row["VL_v_sequence_start"]
            ),
        axis=1)
    df["linker_pass"] = df.apply(
        lambda row: check_linker_qual(
            row["score_pcnt"], 
            row["mismatches"], 
            row["linker_overhang"], 
            score_pcnt, 
            mismatches, 
            linkeroverhang), 
        axis = 1)
    df["productive_scfv"] = df.apply(
        lambda row: check_productive_scfv(
            row["VH_productive"], 
            row["VL_productive"]), 
        axis =  1)
    df["complete_scfv"] = df.apply(
        lambda row: check_productive_scfv(
            row["VH_complete_vdj"], 
            row["VL_complete_vdj"]), 
        axis =  1)
    if chainlen != "undefined":
        df["min_chain_aa_len"] = df.apply(
            lambda row: check_chain_length(
                row["VH_sequence_alignment_aa"], 
                row["VL_sequence_alignment_aa"], 
                chainlen), 
            axis = 1)
    else: 
        df["min_chain_aa_len"] = "NA"
    df.to_csv(
        "in_frame_igBLAST_paired_delim_linker_scored_flags.tsv", 
        sep = "\t"
        )

if __name__ == '__main__':
    main(input, score_pcnt, mismatches, linkeroverhang, chainlen)