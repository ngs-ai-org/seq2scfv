#!/usr/bin/env Rscript

# usage: linklengths.r nt_inferred_linkers_lengths.tsv aa_inferred_linkers_lengths.tsv

library(ggplot2)
library(dplyr)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("No enough arguments provided. Two are required: nt_inferred_linkers_lengths.tsv and aa_inferred_linkers_lengths.tsv", call.=FALSE)
} 


read_data <- function(input_file) {
  df <- read.table(input_file, sep = "\t", header = TRUE)
  return(df)
}

get_quantiles <- function(df) {
    df_5_95 <- df %>%
    filter(length >= quantile(df$length, probs = 0.05)[[1]],
            length <= quantile(df$length, probs = 0.95)[[1]])
    return(df_5_95)

}

plot_linker_lengths <- function(df, df_5_95, seqtype, seqcol) {
    plot1 <- ggplot(df, aes(x = length)) +
        geom_histogram(binwidth = 0.5, fill = seqcol) +
        theme_bw() +
        labs(y = "Count",
            x = paste0("Length (", seqtype, ")"),
            title = paste0("Length distribution of inferred linkers (",
                seqtype, ")")) +
        theme(axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            strip.text = element_text(size = 12),
            title = element_text(size = 12))

    plot2 <- ggplot(df_5_95, aes(x = length)) +
        geom_histogram(binwidth = 0.5, fill = seqcol) +
        theme_bw() +
        labs(y = "Count", 
            x = paste0("Length (", seqtype, ")"),
            title = "Subset from 5th to 95th pctl.") +
        theme(axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            strip.text = element_text(size = 12),
            title = element_text(size = 12))

    grid.arrange(plot1, plot2, nrow = 1, widths = c(1, 0.5))
}

write_q90_data <- function(df_5_95, seqtype) {
    q90_table <- df_5_95 %>% dplyr::count(length)
    colnames(q90_table) <- c(paste0("Length (", seqtype, ")"), "Counts")
    write.table(q90_table,
                paste0(seqtype, "_5_95_percentiles_freq.tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
}

nt <- read_data(args[1])
nt_5_95 <- get_quantiles(nt)
write_q90_data(nt_5_95, "nt")
png(file = "nt_inferred_linkers_length.png", width = 720, height = 480)
plot_linker_lengths(nt, nt_5_95, "nt", "aquamarine4")
dev.off()

aa <- read_data(args[2])
aa_5_95 <- get_quantiles(aa)
write_q90_data(aa_5_95, "aa")
png(file = "aa_inferred_linkers_length.png", width = 720, height = 480)
plot_linker_lengths(aa, aa_5_95, "aa", "steelblue3")
dev.off()