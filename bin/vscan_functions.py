#!/usr/bin/env python3

"""
This script performs iterative IgBLAST analysis on DNA sequences.
It splits sequences based on V, D, and J gene alignments to identify potential V-domains.
"""

# Import necessary libraries
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd 
import os
from optparse import OptionParser
import logging
import subprocess


# Parse input arguments
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("--query", dest="query", default=None, type="string", action="store", help="Input file name for DNA sequences in fasta format.")
parser.add_option("--germline_V", dest="dbv", default=None, type="string", action="store",help="Germline database name for V gene.")
parser.add_option("--germline_D", dest="dbd", default=None, type="string", action="store",help="Germline database name for D gene.")
parser.add_option("--germline_J", dest="dbj", default=None, type="string", action="store",help="Germline database name for J gene.")
parser.add_option("--organism", dest="organism", default=None, type="string", action="store", help="Organism of the query sequence: 'human', 'mouse', 'rat', 'rabbit' or 'custom' (see IgBLAST web site).")
parser.add_option("--ig_seqtype", dest="igtype", default="Ig", type="string", action="store", help="Receptor sequence type: 'Ig' or 'TCR'.")
parser.add_option("--domain_system", dest="domainsys", default="imgt", type="string", action="store", help="Domain system to be used for segment annotation: 'imgt' or 'kabat'.")
parser.add_option("--auxiliary_data", dest="auxdata", default=None, type="string", action="store", help="File containing the coding frame start positions for sequences in germline J database.")
parser.add_option("--num_threads", dest="nthreads", default=1, type=int, action="store", help="Number of threads (CPUs) to use in the BLAST search. Default: 1.")
parser.add_option("--minlen", dest="minlen", default=150, type=int, action="store",help="Minimum length (bp) required to continue seraching for a V-domain. Default: 150 bp.")
parser.add_option("--escore", dest="escore", default=1.0e-02, type=float, action="store",help="Minimum V and J alignment support (E, expectation value) to consider a hit fully characterizing a V-(D)-J region. Default: 0.01.")
parser.add_option("--log", dest="loglevel", type="string", default="info")

(options, args) = parser.parse_args()

input = options.query
germlv = options.dbv
germld = options.dbd
germlj = options.dbj
org = options.organism
ig = options.igtype
domain = options.domainsys
auxdata = options.auxdata
nthreads = options.nthreads
minlen = options.minlen
escore = options.escore
loglevel = options.loglevel

# Set logging level
numeric_level = getattr(logging, loglevel.upper(), None)
if not isinstance(numeric_level, int):    
    raise ValueError('Invalid log level: %s' % loglevel)
logging.basicConfig(level=numeric_level)

logging.debug(options)

# Validate essential options
if input is None:
    raise RuntimeError("Error! no query file provided (--query)")
if germlv is None:
    raise RuntimeError("Error! no germline database for V genes provided (--germline_V)")
if germld is None:
    raise RuntimeError("Error! no germline database for D genes provided (--germline_D)")
if germlj is None:
    raise RuntimeError("Error! no germline database for J genes provided (--germline_J)")
if org is None:
    raise RuntimeError("Error! no organism specified (--organism)")
if auxdata is None:
    raise RuntimeError("Error! no auxiliary file provided (--auxiliary_data)")


# List of relative positions in the sequence
relative_pos = ["v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end", 
                "j_sequence_start", "j_sequence_end", "cdr1_start", "cdr1_end", "cdr2_start", 
                "cdr2_end", "cdr3_start", "cdr3_end", "fwr1_start", "fwr1_end", "fwr2_start", 
                "fwr2_end", "fwr3_start", "fwr3_end", "fwr4_start", "fwr4_end"]

def run_igblast(germlineV, germlineD, germlineJ, organism, igseqtype, threads, domain, auxiliary, queryfile, outfile):
    """
    Run IgBLAST to align sequences to germline genes.
    """
    threads = str(threads)
    command = [
        'igblastn', '-germline_db_V', germlineV, '-germline_db_D', germlineD, 
        '-germline_db_J', germlineJ, '-organism', organism, '-ig_seqtype', igseqtype, 
        '-num_threads', threads, '-domain_system', domain, '-num_alignments_V', '1', 
        '-num_alignments_D', '1', '-num_alignments_J', '1', '-outfmt', '19', 
        '-extend_align5end', '-extend_align3end', '-auxiliary_data', auxiliary, 
        '-query', queryfile, '-out', outfile
    ]
    output = subprocess.check_output(command, universal_newlines=True)
    return output

def igtable(input):
    """
    Read IgBLAST output and add coordinates to the sequence ID.

    Args:
        input (str): Path to IgBLAST output file.

    Returns:
        pd.DataFrame: Formatted IgBLAST output.
    """
    blast = pd.read_csv(input, sep='\t', header=0)
    blast['sequence'] = blast['sequence'].astype(str)
    if not blast["sequence_id"].str.contains('\|').any():
        blast["sequence_id"] = blast.apply(lambda row: row["sequence_id"] + "|0:" + str(len(row["sequence"])), axis=1)
        blast["origin"] = "uncut"
    return blast


def blastfilt(blast, escore, output2):
    """
    Filter IgBLAST output based on E score to select well-delimited V-domains.

    Args:
        blast (pd.DataFrame): IgBLAST output.
        escore (float): Minimum E score.
        output2 (str): Path to save filtered output.

    Returns:
        pd.DataFrame: Filtered IgBLAST output.
    """
    filtered_blast = blast[(blast["v_support"] < float(escore)) & (blast["j_support"] < float(escore))]
    filtered_blast["short"] = filtered_blast["sequence_id"].str.split('|').str[0]
    filtered_blast["coord"] = filtered_blast["sequence_id"].str.split('|').str[1].str.split(':').str[0]
    filtered_blast.reset_index(drop=True, inplace=True)
    filtered_blast.to_csv(output2, sep="\t")
    return(filtered_blast) 

def scan5p(filtered_blast, minlen): 
    """
    Scan the 5' end for sufficient sequence length. Returns sequences meeting length requirement on 5' end.
    """
    seq5p = filtered_blast[filtered_blast["v_sequence_start"] > int(minlen)]
    if not seq5p.empty: 
        seq5p["sequence_id"] = seq5p.apply(lambda row: row["short"] + "|" + str(row["coord"]) + ":" + str(round(row["v_sequence_start"])), axis=1)  
        seq5p["sequence"] = seq5p.apply(lambda row: row["sequence"][:int(row["v_sequence_start"])], axis=1)
    else: 
        seq5p = pd.DataFrame()
    return seq5p

def scan3p(filtered_blast, minlen): 
    """
    Scan the 3' end for sufficient sequence length. Returns sequences meeting length requirement on 3' end.
    """
    seq3p = filtered_blast[(filtered_blast["sequence"].str.len() - filtered_blast["j_sequence_end"]) > int(minlen)]
    if not seq3p.empty: 
        seq3p["sequence_id"] = seq3p.apply(lambda row: row["short"] + "|" + str(round(int(row["j_sequence_end"]) + int(row["coord"]))) + ":" + str(len(row["sequence"])+ int(row["coord"])), axis=1)
        seq3p["sequence"] = seq3p.apply(lambda row: row["sequence"][(int(row["j_sequence_end"])-1):], axis=1)
    else: 
        seq3p = pd.DataFrame()
    return seq3p 


def vscan(filtered_blast, minlen):  
    """
    Combine 5' and 3' end scans.

    Args:
        filtered_blast (pd.DataFrame): Filtered IgBLAST output.
        minlen (int): Minimum length required.

    Returns:
        pd.DataFrame: Combined sequences from 5' and 3' scans.
    """
    seq5p = scan5p(filtered_blast, minlen)
    seq3p = scan3p(filtered_blast, minlen)
    if seq5p.empty and seq3p.empty: 
        newset = pd.DataFrame()
    elif seq5p.empty: 
        newset = seq3p[["sequence_id", "sequence"]]
    elif seq3p.empty: 
        newset = seq5p[["sequence_id", "sequence"]]
    else:
        newset = pd.concat([seq5p[["sequence_id", "sequence"]], seq3p[["sequence_id", "sequence"]]], axis=0)
    return newset

def df_to_fasta(newset):
    """
    Convert DataFrame to FASTA format.

    Args:
        newset (pd.DataFrame): DataFrame of sequences.

    Returns:
        list: List of SeqRecord objects.
    """
    if not newset.empty: 
        newset.reset_index(drop=True, inplace=True)
        newsetlist = list(newset.itertuples(index=False, name=None))
        newfasta = []
        for i in newsetlist:
            newfasta.append(SeqRecord(Seq(i[1]), id=i[0], name="", description=""))
    else: 
        newfasta = []
    return newfasta

def splitfasta(input, minlen, escore, output1, output2):
    """
    Process IgBLAST output and prepare new fragments for search.

    Args:
        input (str): Path to IgBLAST output file.
        minlen (int): Minimum length required.
        escore (float): Minimum E score.
        output1 (str): Path to save new fragments in FASTA format.
        output2 (str): Path to save filtered output.

    Returns:
        None
    """   
    df = igtable(input)
    filt_df = blastfilt(df, escore, output2)
    vscan_df = vscan(filt_df, minlen)
    SeqIO.write(df_to_fasta(vscan_df), output1, "fasta")


def rescanner(germlineV, germlineD, germlineJ, organism, igseqtype, threads, domain, auxiliary, queryfile, minlen, escore): 
    """
    Perform iterative IgBLAST analysis and split sequences. Takes as input IgBLAST args and returns list of paths to filtered output files. 
    """    
    counter = 0 
    res = []
    while (os.stat(queryfile).st_size != 0):
        counter +=1
        igblasttsv = "ddDNA_igblast_" + str(counter) + ".tsv" 
        run_igblast(germlineV, germlineD, germlineJ, organism, igseqtype, threads, domain, auxiliary, queryfile, igblasttsv)
        splitfastatsv = "ddDNA_splitfasta_" + str(counter) + ".tsv"
        queryfile = "ddDNA_splitfasta_" + str(counter) + ".fasta"
        splitfasta(igblasttsv, minlen, escore, queryfile, splitfastatsv) 
        res.append(splitfastatsv)
    return res


def gatherdf(result, changepos): 
    """
    Combine all filtered results and adjust coordinates.

    Args:
        result (list): List of paths to filtered output files.
        changepos (list): List of columns to adjust coordinates.

    Returns:
        pd.DataFrame: Combined and adjusted results.
    """    
    dfs = []
    # read and concat all dfs
    for i in result: 
        df = pd.read_csv(i, sep='\t', header = 0, index_col = 0)
        dfs.append(df)
    result_df = pd.concat(dfs, ignore_index=True)
    # make sure all coordinates in the final table are relative to the OG seq
    for j in changepos: 
        result_df.loc[result_df[j].notna(), j] += result_df.loc[result_df[j].notna(), "coord"]
    # sequence field represents the full length sequence
    fullseq = result_df.groupby("short").filter(lambda x: any(x['origin'] == 'uncut'))
    uncut_mapping = (
        fullseq.loc[fullseq['origin'] == 'uncut']
        .set_index('short')
        ['sequence']
        .to_dict()
    )
    result_df.loc[(result_df['origin'] != 'uncut') & (result_df['short'].isin(uncut_mapping)), 'sequence'] = (
        result_df['short'].map(uncut_mapping)
    )
    # remove coordinate data (original position now taken into account in table)
    result_df = result_df.drop(columns = ["sequence_id", "coord", "origin"])
    newcol = result_df.pop("short")
    result_df.insert(0, "sequence_id", newcol)
    result_df.to_csv("igBLAST.tsv", sep="\t")
    return(result_df)

def main(germlineV, germlineD, germlineJ, organism, igseqtype, threads, domain, auxiliary, queryfile, minlen, escore, relative_pos): 
    res = rescanner(germlv, germld, germlj, org, ig, nthreads, domain, auxdata, input, minlen, escore)
    finaldf = gatherdf(res, relative_pos)    


if __name__ == '__main__':
    main(germlv, germld, germlj, org, ig, nthreads, domain, auxdata, input, minlen, escore, relative_pos)