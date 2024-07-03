#!/usr/bin/env python3

# Import necessary libraries
import logging
import warnings
from Bio.Seq import Seq
import pandas as pd
import re
import multiprocessing as mp
import getID

logger = logging.getLogger(__name__)

id_generator = getID.SeqIDGenerator("seguid")
    
def igBLAST_to_df(input):
    """
    Reads IgBLAST results from a tsv file and converts them into a pandas DataFrame.
    """    
    logger.info("Reading igBLAST table")
    igtsv = pd.read_csv(input, sep = "\t", header = 0, index_col = [0])
    igtsv["sequence"] = igtsv["sequence"].astype(str)
    return(igtsv)


def select_V_region(df, V, prefix):
    """
    Annotates columns with VH or VL prefixes based on the locus type (IGH or IGK).
    
    Args:
    df (pd.DataFrame): Input DataFrame.
    V (str): The locus type (IGH or IGK).
    prefix (str): Prefix to add (VH_ or VL_).
    
    Returns:
    pd.DataFrame: DataFrame with prefixed columns.
    """
    logger.info("Annotating VH/VL")
    V = df[df["locus"] == V].add_prefix(prefix)
    V.columns = [col.replace(prefix, '') if i < 2 else col for i, col in enumerate(V.columns)]
    return(V)

def select_V_pairs(df):
    """
    Filters sequences to obtain those containing exactly 1 VH and 1 VL, writes excluded sequences to a CSV.
    
    Args:
    df (pd.DataFrame): Input DataFrame.
    
    Returns:
    pd.DataFrame: DataFrame containing sequences with 1 VH and 1 VL.
    """

    logger.info("Filtering 1 VH + 1 VL")
    delimit = df[df["v_sequence_start"].notna() & df["j_sequence_end"].notna()]
    delimitpair = delimit.groupby("sequence_id").filter(lambda x : len(x) == 2)
    VHVL1 = pd.merge(select_V_region(delimitpair, "IGH", "VH_"), select_V_region(delimitpair, "IGK", "VL_"), on=['sequence_id', 'sequence'], how='inner')
    VHVL2 = pd.merge(select_V_region(delimitpair, "IGH", "VH_"), select_V_region(delimitpair, "IGL", "VL_"), on=['sequence_id', 'sequence'], how='inner')
    VHVL_combined = pd.concat([VHVL1, VHVL2]).reset_index()
    nonVHVL = delimit[~delimit['sequence_id'].isin(VHVL_combined['sequence_id'])]
    nonVHVL.to_csv("igBLAST_non_paired_nondelim.tsv", sep = "\t")
    return(VHVL_combined)

def find_forward_translation(seq, aa_vh, aa_vl):
    """
    Identifies the correct translation frame for the full sequence containing VH and VL amino acid sequences.
    
    Args:
    seq (str): Nucleotide sequence.
    aa_vh (str): VH amino acid sequence.
    aa_vl (str): VL amino acid sequence.
    
    Returns:
    tuple: Translated amino acid sequence and the translation frame.
    """
    nuc_record, VH, VL = Seq(seq), Seq(str(aa_vh)), Seq(str(aa_vl))
    frame = []
    numb = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(3): 
            fi = str(nuc_record[i:].translate(table = "Standard"))
            if (str(VH) in fi) and (str(VL) in fi):  
                frame.append(fi)
                numb.append(i)
    if (len(frame) == 1): 
        aa = str(frame[0])
        aa_num = int(numb[0])
    elif (len(frame) == 0):
        aa = "not_in_frame"
        aa_num = "NA"
    elif(len(frame) > 1): 
        aa = "incongruent_frames"
        aa_num = "NA"
    return aa, aa_num



def get_coords(seq, aa_sequence, trans_frame, aa_vh, aa_vl):
    """
    Identifies the start and end coordinates of VH and VL regions and annotates scFv and linker regions.
    
    Args:
    seq (str): Nucleotide sequence.
    aa_sequence (str): Amino acid sequence.
    trans_frame (int): Translation frame.
    aa_vh (str): VH amino acid sequence.
    aa_vl (str): VL amino acid sequence.
    
    Returns:
    tuple: Coordinates and sequences of scFv and linker regions in both nucleotide and amino acid forms.
    """
    # Find the VH & VL stat-end positions
    VH = re.search(aa_vh, aa_sequence)
    VL = re.search(aa_vl, aa_sequence)
    if VH and VL: 
        # Write the VH & VL stat-end positions
        aa_VH_start, aa_VH_end = VH.start(), VH.end()
        aa_VL_start, aa_VL_end = VL.start(), VL.end()
        # Write aa scFv delimitations 
        aa_scFv_start = min(aa_VH_start, aa_VL_start)
        aa_scFv_end = max(aa_VH_end, aa_VL_end)
        aa_linker_start = min(aa_VH_end, aa_VL_end) 
        aa_linker_end = max(aa_VH_start, aa_VL_start)
        # Write aa scFv & linker 
        aa_scFv = aa_sequence[aa_scFv_start:aa_scFv_end]
        aa_linker = aa_sequence[aa_linker_start:aa_linker_end]
        # Write nt scFv delimitations
        nt_scFv_start = (int(aa_scFv_start)*3) + trans_frame
        nt_scFv_end = (int(aa_scFv_end)*3) + trans_frame
        nt_linker_start = (int(aa_linker_start)*3) + trans_frame
        nt_linker_end = (int(aa_linker_end)*3) + trans_frame
        # Write nt scFv & linker 
        nt_scFv = seq[nt_scFv_start:nt_scFv_end]
        nt_linker = seq[nt_linker_start:nt_linker_end] 
        # get VH & VL inferred from aa sequence
        nt_VH = seq[((int(aa_VH_start)*3) + trans_frame):((int(aa_VH_end)*3) + trans_frame)]
        nt_VL = seq[((int(aa_VL_start)*3) + trans_frame):((int(aa_VL_end)*3) + trans_frame)]
        # for reporting in AIRR format: 1-based, closed intervals
        aa_scFv_start += 1 
        aa_linker_start += 1 
        nt_scFv_start += 1 
        nt_linker_start += 1 
    else: 
        nt_scFv, nt_linker, nt_scFv_start, nt_scFv_end, \
        nt_linker_start, nt_linker_end, \
        aa_scFv, aa_linker, \
        aa_scFv_start, aa_scFv_end, \
        aa_linker_start, aa_linker_end, \
        nt_VH, nt_VL = "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
    # Get scFv and linkers
    return  (nt_scFv, 
             nt_scFv_start, nt_scFv_end,
             nt_linker,
             nt_linker_start, nt_linker_end,
             nt_VH, nt_VL,
             aa_scFv, 
             aa_scFv_start, aa_scFv_end, 
             aa_linker,
             aa_linker_start, aa_linker_end)

def annotate_row(row_tuple):
    """
    Annotates a row with aa read & scFv and linker sequence and coordinates.
    
    Args:
    row_tuple (tuple): A tuple containing sequence and VH/VL amino acid alignments.
    
    Returns:
    tuple: Annotated data with scFv and linker coordinates.
    """
    (
        sequence,
        VH_sequence_alignment_aa,
        VL_sequence_alignment_aa,
    ) = row_tuple
    aa_read, frame = find_forward_translation(
        seq=sequence,
        aa_vh=str(VH_sequence_alignment_aa),
        aa_vl=str(VL_sequence_alignment_aa)
    )
    nt_scFv, nt_scFv_start, nt_scFv_end, \
    nt_linker, nt_linker_start, nt_linker_end, \
    nt_VH, nt_VL, \
    aa_scFv, aa_scFv_start, aa_scFv_end, \
    aa_linker, aa_linker_start, aa_linker_end = get_coords(seq = sequence,
                                                           aa_sequence = aa_read, 
                                                           trans_frame = frame, 
                                                           aa_vh=str(VH_sequence_alignment_aa),
                                                           aa_vl=str(VL_sequence_alignment_aa))
    nt_scFv_id = id_generator.encode(nt_scFv) + "_nt_scFv"
    aa_scFv_id = id_generator.encode(aa_scFv) + "_aa_scFv"
    return (
        nt_scFv_id, nt_scFv, nt_scFv_start, nt_scFv_end, 
        nt_linker, nt_linker_start, nt_linker_end,
        nt_VH, nt_VL,
        aa_read, frame,         
        aa_scFv_id, aa_scFv, aa_scFv_start, aa_scFv_end,
        aa_linker, aa_linker_start, aa_linker_end
    )


def find_aa_boundaries(row_tuple):
    """
    Finds amino acid boundaries for FWR and CDR regions within VH and VL sequences.
    
    Args:
    row_tuple (tuple): A tuple containing sequence ID, VH, and VL sequence alignments and regions.
    
    Returns:
    dict: Dictionary containing start and end positions for FWR and CDR regions.
    """
    (
        seqid,
        VH_sequence_alignment_aa,
        VH_fwr1_aa, 
        VH_cdr1_aa, 
        VH_fwr2_aa, 
        VH_cdr2_aa, 
        VH_fwr3_aa, 
        VH_cdr3_aa, 
        VH_fwr4_aa,
        VL_sequence_alignment_aa,
        VL_fwr1_aa, 
        VL_cdr1_aa, 
        VL_fwr2_aa, 
        VL_cdr2_aa, 
        VL_fwr3_aa, 
        VL_cdr3_aa, 
        VL_fwr4_aa 
    ) = row_tuple
    
    aa_chain_regions = {
        'VH': [VH_fwr1_aa, VH_cdr1_aa, VH_fwr2_aa, VH_cdr2_aa, VH_fwr3_aa, VH_cdr3_aa, VH_fwr4_aa],
        'VL': [VL_fwr1_aa, VL_cdr1_aa, VL_fwr2_aa, VL_cdr2_aa, VL_fwr3_aa, VL_cdr3_aa, VL_fwr4_aa]
    }
    
    tuple_pos = {'sequence_id': seqid}
    
    for chain, aa_seq in [('VH', VH_sequence_alignment_aa), ('VL', VL_sequence_alignment_aa)]:
        regions = ['fwr1_aa', 'cdr1_aa', 'fwr2_aa', 'cdr2_aa', 'fwr3_aa', 'cdr3_aa', 'fwr4_aa']
        for region, region_aa in zip(regions, aa_chain_regions[chain]):
            start_key = f"{chain}_{region}_start"
            end_key = f"{chain}_{region}_end"
            match = re.search(re.escape(region_aa), aa_seq)
            if match:
                tuple_pos[start_key] = match.start() + 1 # to match AIRR format 1-based
                tuple_pos[end_key] = match.end() # closed interval
    return tuple_pos


def set_aa_boundaries(df):
    """
    Sets amino acid boundaries for FWR and CDR regions for all rows in the DataFrame.
    
    Args:
    df (pd.DataFrame): Input DataFrame.
    
    Returns:
    pd.DataFrame: DataFrame with updated boundaries.
    """
    df.fillna('', inplace=True)
    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(
            find_aa_boundaries, 
            df[['sequence_id',
                'VH_sequence_alignment_aa',
                'VH_fwr1_aa', 
                'VH_cdr1_aa', 
                'VH_fwr2_aa', 
                'VH_cdr2_aa', 
                'VH_fwr3_aa', 
                'VH_cdr3_aa', 
                'VH_fwr4_aa',
                'VL_sequence_alignment_aa',
                'VL_fwr1_aa', 
                'VL_cdr1_aa', 
                'VL_fwr2_aa', 
                'VL_cdr2_aa', 
                'VL_fwr3_aa', 
                'VL_cdr3_aa', 
                'VL_fwr4_aa']].itertuples(index=False, name=None))
    
    newcols_df = pd.DataFrame(results)
    return newcols_df

def eval_consistency(df):
    """
    Evaluates the consistency between IgBLAST annotations and Antpack numbering for VH and VL sequences.
    
    Args:
    df (pd.DataFrame): Input DataFrame.
    
    Returns:
    pd.DataFrame: DataFrame with consistency evaluations added.
    """

    logger.info('Evaluating consistency in annotation: IgBLAST vs Antpack numbering')
    # VH and VL length: numbering vs alignment
    df['VH_len_consistency'] = df.apply(lambda row: 'T' 
                                        if len(row['VH_numbering']) == len(row['VH_sequence_alignment_aa']) 
                                        else 'F', axis=1)
    df['VL_len_consistency'] = df.apply(lambda row: 'T' 
                                        if len(row['VL_numbering']) == len(row['VL_sequence_alignment_aa']) 
                                        else 'F', axis=1)
    # VH and VL starting aa: numbering vs alignment
    df['VH_start_consistency'] = df.apply(lambda row: 'T' 
                                          if row['VH_numbering'] and row['VH_numbering'][0] == '1' 
                                          else 'F' if row['VH_numbering'] else 'NA', axis=1)
    df['VL_start_consistency'] = df.apply(lambda row: 'T' 
                                          if row['VL_numbering'] and row['VL_numbering'][0] == '1' 
                                          else 'F' if row['VL_numbering'] else 'NA', axis=1)
    df.loc[(df['VH_start_consistency'] == 'F') & 
           (df['nt_scFv_start'] == 0), 'VH_start_consistency'] = 'Inferred incomplete sequencing of chain'
    return df


def update_coordinates(df): 
    """
    Updates FWR and CDR coordinates to be relative to the start of each VH/VL sequence.
    
    Args:
    df (pd.DataFrame): Input DataFrame.
    
    Returns:
    pd.DataFrame: DataFrame with updated coordinates.
    """
    chains = ['VH', 'VL']
    nt_regions = ['fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4']
    for chain in chains: 
        fwr1_initial_pos = df[f'{chain}_fwr1_start']
        for region in nt_regions: 
            coordinates = ['_'.join([chain, region, 'start'])] + ['_'.join([chain, region, 'end'])]
            for position in coordinates: 
                df[position] = df[position] - fwr1_initial_pos + 1
    return df