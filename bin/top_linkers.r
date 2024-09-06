#!/usr/bin/env Rscript

# usage: top_linkers.r  in_frame_igBLAST_paired_delim.tsv aa_reference_linker.fasta

library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("No enough arguments provided. Two are required: in_frame_igBLAST_paired_delim.tsv and aa_reference_linker.fasta", call.=FALSE)
} 

df <- read.table(
                 args[1],
                 fill = TRUE,
                 header = TRUE, 
                 sep = "\t",
                 row.names = NULL)

reflink <- read.table(
                      args[2],
                      header = FALSE, 
                      row.names = NULL)[2, 1]

aa_linkers <- df %>%
  dplyr::count(aa_linker, sort = TRUE) %>%
  mutate(freq = round((n / sum(n)) * 100, 2)) %>%
  head(10) %>%
  mutate(ref_linker_present = str_detect(aa_linker, as.character(reflink)))

nt_linkers <- df %>%
  dplyr::count(nt_linker, sort = TRUE) %>%
  mutate(freq = round((n / sum(n)) * 100, 2)) %>%
  head(10)

write.table(aa_linkers,
            file = "aa_top10_linkers.tsv",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")


write.table(nt_linkers,
            file = "nt_top10_linkers.tsv",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
