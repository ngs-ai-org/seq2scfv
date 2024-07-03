#!/usr/bin/env python3

import sys
import typing as tp
from Bio import SeqIO
import getID


def main() -> None:
    """
    Generate unique sequence IDs for sequences in a FASTA file.
    
    It reads a FASTA file, generates unique IDs using SEGUID, writes correspondence 
    information to a TSV file, and outputs the modified sequences to a new FASTA file.
    """
    lib_name = sys.argv[1]
    fasta = sys.argv[2]
    corr_file = open(f"{lib_name}.nt_correspondence.tsv", "w")
    out_fasta = open(f"{lib_name}.nt_seguid.fa", "w")
    id_gen = getID.SeqIDGenerator("seguid")
    has_seen = {}
    for rec in SeqIO.parse(fasta, "fasta"):
        seguid =id_gen.encode(rec.seq)
        new_id = f"{seguid}_nt_qual"
        print(new_id, rec.id, lib_name, sep="\t", file=corr_file)
        try:
            _ = has_seen[seguid]
        except KeyError:
            rec.id = new_id
            rec.description = ""
            print(rec.format("fasta"), file=out_fasta, end="")
            has_seen[seguid] = True

if __name__ == "__main__":
    main()
