rm(list = ls())

library(dplyr)
library(data.table)
library(tidyr)


args <- commandArgs(trailingOnly = TRUE)

# test if there are two arguments: if not, return an error
if (length(args) != 1) {
  stop("Please supply degradation rates file", call. = FALSE)
}

if (file.exists(args[1])) {
  path <- args[1]
} else {
  stop(paste("file", args[1], "doesn't exist"))
}

infile <- fread(path)

formatted_file <- infile %>%
  select(id, 'half-life') %>%
  rename(gene_name = id, param_val = 'half-life') %>%
  mutate(param_name = 'half_life',
         category = 'all')

fwrite(formatted_file, gsub(".txt", ".kmer.tsv", path), sep = '\t')
  