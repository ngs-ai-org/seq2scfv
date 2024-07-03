#!/usr/bin/env python3

"""
scFv identification and characterization. 
"""

import logging
import multiprocessing as mp
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import scFv_tools
from antpack import SingleChainAnnotator
import re

logger = logging.getLogger(__name__)

# Parse input arguments
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("--igBLAST-hits", dest="hits", default=None, type="string", action="store", help="Pre-filtered IG hits stemming from V-scan.")
parser.add_option("--domain-type", dest="domain", default=None, type="string", action="store", help="Domain system to be used for segment annotation: 'imgt' or 'kabat'")
parser.add_option("--log", dest="loglevel", type="string", default="info")

(options, args) = parser.parse_args()

input = options.hits
domain = options.domain
loglevel = options.loglevel

# Set logging level
numeric_level = getattr(logging, loglevel.upper(), None)
if not isinstance(numeric_level, int):    
    raise ValueError('Invalid log level: %s' % loglevel)
logging.basicConfig(level=numeric_level)

logging.debug(options)

# Validate essential options
if input is None:
    raise RuntimeError("Error! no IgBLAST result file given to process (--query)")
if domain is None:
    raise RuntimeError("Error! no domain system specified (--domain-type)")

# compile all information in df
def annotator(paired_V):
    """
    Adds scFv and linker sequences to the dataframe.

    Args:
    paired_V (pd.DataFrame): DataFrame containing paired VH and VL sequences.

    Returns:
    pd.DataFrame: Updated DataFrame with scFv and linker sequences annotated.
    """
    logging.info("Adding scFv and linker sequences")
    # filter out sequences with stops 
    logging.debug(f"Shape of scFv table after linkers: {paired_V.shape}")
    with mp.Pool(mp.cpu_count()) as pool:
        logging.info("Running MP row annotations")
        paired_V[[ 
            "nt_scFv_id",
            "nt_scFv",
            "nt_scFv_start",
            "nt_scFv_end",
            "nt_linker",
            "nt_linker_start",
            "nt_linker_end",
            "nt_VH", 
            "nt_VL",
            "aa_sequence",
            "frame",
            "aa_scFv_id",
            "aa_scFv",
            "aa_scFv_start",
            "aa_scFv_end",
            "aa_linker",
            "aa_linker_start",
            "aa_linker_end"
        ]] = pool.map(
            scFv_tools.annotate_row, # scFv_tools.annotate_row
            paired_V[[
                "sequence",
                "VH_sequence_alignment_aa",
                "VL_sequence_alignment_aa",
            ]].itertuples(index=False, name=None))
    logging.debug(f"Shape of scFv table after annotations: {paired_V.shape}")
    return paired_V

# Use Antpack to number chains
def number_seq(aa):
    """
    Numbers amino acid sequences using the Antpack annotator.

    Args:
    aa (str): Amino acid sequence to be numbered.

    Returns:
    pd.Series: Series containing the numbering, percent identity, chain type, and error message.
    """
    numbering, percent_identity, chain_type, err_message = numbering_annotator.analyze_seq(aa)
    return pd.Series({
        'numbering': numbering,
        'percent_identity': percent_identity,
        'chain_type': chain_type,
        'err_message': err_message
    })

def number_df(df):
    """
    Numbers amino acid sequences for all VH and VL sequences in the DataFrame.

    Args:
    df (pd.DataFrame): Input DataFrame with VH and VL sequences.

    Returns:
    pd.DataFrame: Updated DataFrame with numbered sequences.
    """
    logger.info('Numbering aa with anpack')
    with mp.Pool(mp.cpu_count()) as pool:
        logging.info('Running MP row annotations')
        df[['VH_numbering', 
            'VH_percent_identity', 
            'VH_chain_type', 
            'VH_err_message']]  = pool.map(
                number_seq, 
                df['VH_sequence_alignment_aa'])
    with mp.Pool(mp.cpu_count()) as pool:
        logging.info('Running MP row annotations')
        df[['VL_numbering', 
            'VL_percent_identity', 
            'VL_chain_type', 
            'VL_err_message']]  = pool.map(
                number_seq, 
                df['VL_sequence_alignment_aa'])
    logging.debug(f'Shape of table after numbering: {df.shape}')
    return df

# write subsequences to csv
def subSeq(
    annot_vhvl: pd.DataFrame
) -> None:
    """Helper function to write out the FASTA files from the annotation table
    
    The annotation table, 'annot_vhvl', is a Pandas dataframe holding the
    information on VH, VL and linker parts. Each row corresponds to a
    'nt_qual' identifier - deduplicated high-quality nt sequences containing
    a VH and a VL. Therefore, different rows can contain the same scFv (either
    nt or aa). This function writes out 4 fasta files:
    
    1. nt_scFv: Deduplicated scFv nt sequences (nt_qual ids are in description)
    2. aa_scFv: Deduplicated scFv aa sequences (nt_qual ids are in description)
    3. nt_linker: Linker nt sequences (not dedup) from deduplicated scFv nt
    4. aa_linker: Linker aa sequences (not dedup) from deduplicated scFv aa
    
    Args:
        annot_vhvl: A Pandas dataframe returned by 'annotator' function
    """
    logging.info("Writing fasta files")
    df_filt = annot_vhvl[annot_vhvl["frame"] != "NA"]
    df = df_filt[
        ["sequence_id",
         "frame",
         "nt_scFv_id",
         "nt_scFv",
         "aa_scFv_id",
         "aa_scFv",
         "nt_linker",
         "aa_linker"]
    ]
    ntscfv, aascfv, ntlinker, aalinker = {}, {}, [], []
    for row in df.itertuples(index=False):
        # deduped at scFv nt level
        if not row.nt_scFv_id in ntscfv:
            ntscfv[row.nt_scFv_id] = SeqRecord(
                seq=Seq(row.nt_scFv),
                id=row.nt_scFv_id,
                name="",
                description=row.sequence_id
            )
            if row.nt_linker: 
                ntlinker.append(
                    SeqRecord(
                        seq=Seq(row.nt_linker),
                        id=f"{row.nt_scFv_id}_linker",
                        name="",
                        description=""
                    )
                )
        else:
            ntscfv[row.nt_scFv_id].description += f";{row.sequence_id}"
        # deduped at scFv aa level
        if not row.aa_scFv_id in aascfv:
            aascfv[row.aa_scFv_id] = SeqRecord(
                seq=Seq(row.aa_scFv),
                id=row.aa_scFv_id,
                name="",
                description=row.sequence_id
            )
            if row.aa_linker:
                aalinker.append(
                    SeqRecord(
                        seq=Seq(row.aa_linker),
                        id=f"{row.aa_scFv_id}_linker",
                        name="",
                        description=""
                    )
                )
        else:
            aascfv[row.aa_scFv_id].description += f";{row.sequence_id}"
    SeqIO.write(ntscfv.values(), "nt_inferred_scFv.fasta", "fasta")
    SeqIO.write(aascfv.values(), "aa_inferred_scFv.fasta", "fasta")
    SeqIO.write(ntlinker, "nt_inferred_linkers.fasta", "fasta")
    SeqIO.write(aalinker, "aa_inferred_linkers.fasta", "fasta")



# Prep the antibody annotator
numbering_annotator = SingleChainAnnotator(['H', 'K', 'L'], scheme = domain)
# Prep aa regions
aa_chain_regions = ['VH_fwr1_aa', 'VH_cdr1_aa', 'VH_fwr2_aa', 'VH_cdr2_aa', 'VH_fwr3_aa', 'VH_cdr3_aa', 'VH_fwr4_aa', 
                    'VL_fwr1_aa', 'VL_cdr1_aa', 'VL_fwr2_aa', 'VL_cdr2_aa', 'VL_fwr3_aa', 'VL_cdr3_aa', 'VL_fwr4_aa']


def main(input):
    """
    Main function to process the igBLAST hits and annotate scFv sequences.

    Args:
    input (str): Path to the input igBLAST file.
    """
    # Read input
    df = scFv_tools.igBLAST_to_df(input)
    # Get dataframe with sequences containing 1 VH & 1 VL; write "orphan" chains to tsv
    paired_V_regions = scFv_tools.select_V_pairs(df)
    # Keep only sequences without * in any of the chains 
    no_stop_paired_V = paired_V_regions[~(paired_V_regions['VH_sequence_alignment_aa'].str.contains('\*') | paired_V_regions['VL_sequence_alignment_aa'].str.contains('\*'))]
    # Annotate scFv & linker at the nt and aa levels
    scFv_delimited = annotator(no_stop_paired_V)
    # Only sequences that are in frame are considered from here on
    # Filter out of frame if frame == NA or if linker contains * | some corner cases cannot be caught with annotator
    scFv_in_frame = scFv_delimited[(scFv_delimited["frame"] != "NA") & (~scFv_delimited["aa_linker"].str.contains('\*'))]
    # Delimit FWR / CDR regions at the aa level, relative to each chain
    region_aa = scFv_tools.set_aa_boundaries(scFv_in_frame.dropna(subset = aa_chain_regions, how = 'any')) 
    scFv_in_frame_aa_pos = scFv_in_frame.merge(region_aa, how = 'left', on = 'sequence_id')
    # Run numbering with Antpack
    scFv_numbered = number_df(scFv_in_frame_aa_pos)
    # Compare IgBlast annotation vs Antpack numbering
    # 4 columns are created: VH & VL length and start position consistency between alignment / numbering
    # T / F annotation, with possibility of "Inferred incomplete sequencing of chain" for VH_start_consistency
    df_consistent = scFv_tools.eval_consistency(scFv_numbered[(scFv_numbered['VH_err_message'] != 'Invalid sequence supplied -- nonstandard AAs') |
                                                   (scFv_numbered['VL_err_message'] != 'Invalid sequence supplied -- nonstandard AAs')])
    # Update CDR / FWR region start / end at the nt level to be relative to the chain (not the read)
    df_consistent_nt_coord = df_consistent.apply(scFv_tools.update_coordinates, axis=1)
    # Add information of numbering system
    df_consistent_nt_coord.rename(columns={'VH_numbering': str('VH_' + domain + '_numbering')}, inplace=True)
    df_consistent_nt_coord.rename(columns={'VL_numbering': str('VL_' + domain + '_numbering')}, inplace=True)
    # Write fasta files 
    subSeq(df_consistent_nt_coord)
    # Write tables 
    scFv_numbered[(scFv_numbered['VH_err_message'] == 'Invalid sequence supplied -- nonstandard AAs') |
                  (scFv_numbered['VL_err_message'] == 'Invalid sequence supplied -- nonstandard AAs')].to_csv('in_frame_igBLAST_delim_nonstandard_aas.tsv', sep='\t', index=False)
    df_consistent_nt_coord.to_csv("in_frame_igBLAST_paired_delim.tsv", sep="\t")    
    out_of_frame1 = paired_V_regions[(paired_V_regions['VH_sequence_alignment_aa'].str.contains('\*') | paired_V_regions['VL_sequence_alignment_aa'].str.contains('\*'))]
    out_of_frame2 = scFv_delimited[(scFv_delimited["frame"] == "NA") | (scFv_delimited["aa_linker"].str.contains('\*'))]
    pd.concat([out_of_frame1, out_of_frame2], ignore_index=True, sort=False).to_csv('out_of_frame_igBLAST_paired_delim.tsv', sep='\t', index=False)
    
if __name__ == '__main__':
    main(input)
