#!/usr/bin/env python3

# Import necessary libraries
import pandas as pd 
from optparse import OptionParser
import logging
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
                
# Parse input arguments
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("--scFv", dest="hits", default=None, type="string", action="store", help="TSV file containing annotations of all in-frame scFv.")
parser.add_option("--linker", dest="linker", default=None, type="string", action="store",help="Reference linker sequence.")
parser.add_option("--log", dest="loglevel", type="string", default="info")

(options, args) = parser.parse_args()

scfvs = options.hits
deduced_linker = options.linker
loglevel = options.loglevel

# Set logging level
numeric_level = getattr(logging, loglevel.upper(), None)
if not isinstance(numeric_level, int):    
    raise ValueError('Invalid log level: %s' % loglevel)
logging.basicConfig(level=numeric_level)

logging.debug(options)

# Validate essential options
if scfvs is None:
    raise RuntimeError("Error! no scFv dataframe provided (--scFv)")
if deduced_linker is None:
    raise RuntimeError("Error! no reference aa linker sequence provided (--linker)")

def read_scfv_table(input): 
    """
    Reads the scFv annotation table from a tsv file.

    Args:
    input (str): Path to the input tsv file.

    Returns:
    pd.DataFrame: DataFrame containing the scFv annotations.
    """
    annot = pd.read_csv(input, sep = "\t", header = 0, index_col = [0])
    return(annot)

def scoring(df, consensus, aligner):
    """
    Scores the alignment of linker sequences against the reference sequence.

    Args:
    df (pd.DataFrame): DataFrame containing scFv annotations with linker sequences.
    consensus (str): Consensus sequence to align against.
    aligner (Align.PairwiseAligner): Aligner object for scoring and aligning sequences.

    Returns:
    pd.DataFrame: Updated DataFrame with scores, mismatches, and linker overhangs.
    """
    scores = {}
    mismatches = {}
    overhangs = {}
    for i, row in df.iterrows():
        if pd.notna(row["aa_linker"]): 
            aa_linker = row["aa_linker"]
            if aa_linker not in scores:
                raw_score = aligner.score(consensus, aa_linker)
                if raw_score == 0:
                    score = raw_score
                    mismatch = "NA"
                    linker_overhang = "NA"
                else: 
                    score = (raw_score/len(consensus))*100
                    aligned_query = (aligner.align(consensus, aa_linker)[0]).aligned[1]
                    aligned_query_length = aligned_query[0, 1] - aligned_query[0, 0]
                    mismatch = aligned_query_length - raw_score
                    linker_overhang = len(aa_linker) - aligned_query_length   
                scores[aa_linker] = score
                mismatches[aa_linker] = mismatch
                overhangs[aa_linker] = linker_overhang
            else:
                score = scores[aa_linker]
                mismatch = mismatches[aa_linker]
                linker_overhang = overhangs[aa_linker]
        else:
            score = "NA"
            mismatch = "NA"
            linker_overhang = "NA"
        df.at[i, "score_pcnt"] = score
        df.at[i, "mismatches"] = mismatch
        df.at[i, "linker_overhang"] = linker_overhang
    return df

def define_aligner(): 
    """
    Defines and configures a pairwise aligner for local alignment.

    Returns:
    Align.PairwiseAligner: Configured pairwise aligner.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local' 
    aligner.mismatch_score = 0
    aligner.query_internal_open_gap_score = -10
    aligner.query_internal_extend_gap_score = -0.5
    aligner.target_internal_open_gap_score = -10
    aligner.target_internal_extend_gap_score = -0.5
    aligner.target_left_open_gap_score = -10
    aligner.target_left_extend_gap_score = -0.5
    aligner.target_right_open_gap_score = -10
    aligner.target_right_extend_gap_score = -0.5  
    aligner.query_left_open_gap_score = -10 
    aligner.query_left_extend_gap_score = -0.5
    aligner.query_right_open_gap_score = -10
    aligner.query_right_extend_gap_score = -0.5 
    return(aligner)


def main(input, linker):
    """
    Args:
    input (str): Path to the input scFv annotation table.
    linker (str): Path to the FASTA file containing the linker consensus sequence.

    Returns:
    pd.DataFrame: DataFrame with updated scFv annotations and scores.
    """
    consensus = str(SeqIO.read(linker, "fasta").seq)
    annot = read_scfv_table(input)
    aligner = define_aligner()
    annotscore = scoring(annot, consensus, aligner)
    annotscore.to_csv("in_frame_igBLAST_paired_delim_linker_scored.tsv", sep = "\t")
    return(annotscore)
    
    
if __name__ == '__main__':
    main(scfvs, deduced_linker)