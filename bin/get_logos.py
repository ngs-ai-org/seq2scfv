#!/usr/bin/env python3

# Import necessary libraries
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from optparse import OptionParser
import logging
import logomaker as lm

pd.options.mode.chained_assignment = None  # default='warn'

# Parse input arguments
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)

parser.add_option("--nt-aln", 
                  dest="nt_input", 
                  default=None, 
                  type="string", 
                  action="store", 
                  help="nt linker alignment to generate logo (.aln)")


parser.add_option("--aa-aln", 
                  dest="aa_input", 
                  default=None, 
                  type="string", 
                  action="store", 
                  help="aa linker alignment to generate logo (.aln)")

parser.add_option("--check-linker", 
                  dest="linkcheck", 
                  default=None, 
                  type="string", 
                  action="store", 
                  help="User-specified amino-acid linker sequence (string) to compare to inferred consensus linker. If not known: 'undefined'")

parser.add_option("--log", dest="loglevel", type="string", default="info")

(options, args) = parser.parse_args()

nt = options.nt_input
aa = options.aa_input
linkercheck = options.linkcheck
loglevel = options.loglevel

# Set logging level
numeric_level = getattr(logging, loglevel.upper(), None)
if not isinstance(numeric_level, int):    
    raise ValueError('Invalid log level: %s' % loglevel)
logging.basicConfig(level=numeric_level)

logging.debug(options)

# Validate essential options
if nt is None:
    raise RuntimeError("Error! no nt linker alignment given to process (--nt-aln)")
if aa is None:
    raise RuntimeError("Error! no aa linker alignment givengiven to process (--aa-aln)")
if linkercheck is None:
    raise RuntimeError("Error! no linker sequence given to process (--check-linker)")

def plot_logo(input, output):
    """Generate logos from nt & aa MSA of subsampled linkers    
    Args:
        input: .aln file of nt or aa linker MSA
        output: file name (.pdf) for saving the figure
    """
    # parse input file 
    aln = []
    for record in SeqIO.parse(input, "fasta"): 
        seq = str(record.seq)
        aln.append(seq)

    total_seqs = len(aln)
    
    # create count matrix
    counts_df = lm.alignment_to_matrix(sequences=aln, 
                                        to_type='counts', 
                                        characters_to_ignore='.-X')

    # filter positions based on counts
    # plot only positions with > 50% of sequences present
    num_seqs = counts_df.sum(axis=1)
    pos_to_keep = num_seqs > len(aln)/2
    filt_counts_df = counts_df[pos_to_keep]
    filt_counts_df.reset_index(drop=True, inplace=True)
    
    # generate the sequence logo & save
    logo = lm.Logo(filt_counts_df)
    plt.savefig(output)  # e.g. aa_linker_logo.pdf
    
    counts_df['Weight'] = round(num_seqs / total_seqs, 4)

    return counts_df


def get_consensus(logo):
    """Detect the nt & aa consensus linker sequence from logodata output
    
    This function uses the weights, defined here as the fraction of non-gap symbols,
    to select the most significant positions in the alignment. Then, for each position, 
    the consensus is defined by the most frequent nucleotide. 

    Args:
        logo: dataframe returned by plot_logo function
    """
    logodf_080 = logo[logo["Weight"] > 0.80]
    logodf_080["maxcount"] = logodf_080.idxmax(axis=1)
    logoseq = ''.join(logodf_080["maxcount"]).replace(" ", "")
    return logoseq

def write_fasta(sequence, name): 
    """Write fasta files
    
    Args: 
        sequence (str): linker sequence deduced by get_consensus function or provided by user with check_linker arg
        name (str): identifier to use for the fasta header and filename
    """ 
    SeqIO.write(
        SeqRecord(
            seq = Seq(sequence), 
            id = name, 
            name="", 
            description=""), 
        (name + ".fasta"), 
        "fasta"
        )    

def main(nt_linkers, aa_linkers, linkercheck): 
    """Detect the nt & aa consensus linker sequence from LogoMaker output
    
    Verification of the correspondence between the inferred linker sequence and the 
    user-provided linker sequence (at the amino acid level) is optional and only issues a warning. 
     
    Args:
        logo: .tsv file returned by weblogo3 (after trimming of annotation lines ('#'))
        linkercheck (string): user-provided amino acid linker sequence (optional)
    """
    nt_df = plot_logo(nt_linkers, "nt_linkers_logo.pdf")
    aa_df = plot_logo(aa_linkers, "aa_linkers_logo.pdf")
    
    nt_df.to_csv("nt_linkers_logo.tsv", sep = "\t")
    aa_df.to_csv("aa_linkers_logo.tsv", sep = "\t")
    
    nt_linker_consensus = get_consensus(nt_df)
    aa_linker_consensus = get_consensus(aa_df)
    
    write_fasta(nt_linker_consensus, "nt_inferred_consensus_linker")
    write_fasta(aa_linker_consensus, "aa_inferred_consensus_linker") 
    
   
    if (str(linkercheck) == "undefined"): 
        write_fasta(aa_linker_consensus, "aa_reference_linker")
        notification = (
            "No amino acid linker sequence was provided by the user. \n"
            "The inferred consensus linker sequence has been determined as the reference."
        )
    else:
        write_fasta(linkercheck, "aa_reference_linker")
        if str(Seq(aa_linker_consensus)) == str(linkercheck):
            notification = "The inferred and the provided amino acid linker sequences are identical. \n"
        else: 
            notification = (
                "The inferred and the provided amino acid linker sequences are not identical.\n"
                f"Provided linker sequence: {str(linkercheck)}\n"
                f"Inferred consensus linker sequence: {str(aa_linker_consensus)}\n"
            )
    with open("inferred_vs_provided_linker_evaluation.txt", "w") as file:
        file.write(notification) 

    
if __name__ == '__main__':
    main(nt, aa, linkercheck)