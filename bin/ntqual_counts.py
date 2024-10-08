#!/usr/bin/env python3

"""
This script combines data from the file containing scFv annotations and the file containing 
correspondence information about uniquely-tagged libraries. It counts occurrences of each 
library from a specified list of focus libraries within the correspondence file data.
"""

# Import necessary libraries
from optparse import OptionParser
import logging
import pandas as pd 

# Parse input arguments
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("--scFv", dest="hits", default=None, type="string", action="store", help="TSV file containing annotations of all in-frame scFv returned by flagging process.")
parser.add_option("--correspondence", dest="corr", default=None, type="string", action="store", help="Correspondence file generated by merging of uniquely-tagged libraries.")
parser.add_option("--focus-libraries", dest="focuslibs", default=90, type="string", action="store", help="Text file with the basename of all provided libraries.")
parser.add_option("--log", dest="loglevel", type="string", default="info")

(options, args) = parser.parse_args()

input = options.hits
correspondence = options.corr
libraries = options.focuslibs
loglevel = options.loglevel

# Set logging level
numeric_level = getattr(logging, loglevel.upper(), None)
if not isinstance(numeric_level, int):    
    raise ValueError('Invalid log level: %s' % loglevel)
logging.basicConfig(level=numeric_level)

logging.debug(options)

# Validate essential options
if input is None:
    raise RuntimeError("Error! No TSV file containing scfv_flags output given (--scFv)")
elif correspondence is None:
    raise RuntimeError("Error! No correspondence file of uniquely-tagged libraries given (--correspondence)")
elif libraries is None:
    raise RuntimeError("Error! No .txt file with list of focus libraries given (--focus-libraries)")

def parse_focuslibs(libs): 
    """
    Parse the list of focus libraries from a text file.

    Parameters:
    libs (str): Path to the text file containing focus libraries.

    Returns:
    list: List of focus libraries.
    """
    with open(libs, "r") as focuslibs:
        lines = focuslibs.readlines()
        formatted_lines = [line.strip() for line in lines]
        return formatted_lines

def count_and_map(row, given_libs):
    """
    Count occurrences of each focus library in a DataFrame row and map the counts.

    Parameters:
    row (pd.Series): Row of a DataFrame containing library information.
    given_libs (list): List of focus libraries.

    Returns:
    pd.Series: Updated row with counts mapped to corresponding focus libraries.
    """
    if row['library'] in given_libs: 
        row[row['library']] = 1
        return row
    else: 
        return row

def count_table(corrtable, libs): 
    """
    Count occurrences of each focus library in a correspondence table.

    Parameters:
    corrtable (str): Path to the correspondence table file.
    libs (str): Path to the text file containing focus libraries.

    Returns:
    pd.DataFrame: DataFrame summarizing counts of focus libraries per correspondence category.
    """
    given_libs = parse_focuslibs(libs)
    data = pd.read_csv(corrtable, 
                       header = None, 
                       sep = "\t")       
    df = pd.DataFrame(data)
    df.columns =['nt_qual', 'ccs', 'library']
    df[[given_libs]] = 0
    df = df.apply(lambda row: count_and_map(row, given_libs),axis=1)
    df = df.drop(['ccs', 'library'], axis = 1)
    df_sum = df.groupby(['nt_qual']).sum()
    df_sum["total_nt_qual"] = df_sum.sum(axis = 1)
    df_sum = df_sum.reset_index()
    df_sum = df_sum.rename({'nt_qual':'sequence_id'}, axis = 1)
    return(df_sum)

def read_scfv_table(input): 
    """
    Read the scFv annotation table from a TSV file.
    """
    annot = pd.read_csv(input, sep = "\t", header = 0, index_col = [0])
    return(annot)

def main(input, correspondence, libraries): 
    df_corr = count_table(correspondence, libraries)
    df_scfv = read_scfv_table(input)
    df_full = pd.merge(df_scfv, df_corr, how = "left")
    df_full.to_csv("in_frame_igBLAST_paired_delim_linker_scored_flags_counts.tsv", sep = "\t")


if __name__ == '__main__':
    main(input, correspondence, libraries)
    
    

