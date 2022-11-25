#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Making tidyverse functions available. We need this for str_remov and piping.
library(tidyverse)

# Reading in the converted and reformatted EMBL file that contains all novel
# sequences.
full_embl <- read.delim(args[1], header = F)

# Determine where each record in the EMBL file starts and stops, and save that
# as a sort of roadmap for the code downstream
seq_limits <- data.frame(deflines_indices = which(grepl("ID", full_embl$V1, fixed = T)),
                         seq_ends = c(which(grepl("ID", full_embl$V1, fixed = T))[-1]-1, nrow(full_embl)))

# Separate out each record, and make sure the accession line is properly
# formatted with no special characters. We will use that accession line to name
# the EMBL files.
for (i in 1:nrow(seq_limits)){
  
  sub <- data.frame(V1 = full_embl[seq_limits$deflines_indices[i]:seq_limits$seq_ends[i],])
  
  ac_line <- sub$V1[which(grepl("AC", sub$V1))] %>%
    str_remove("AC   ") %>%
    str_remove(";") %>%
    str_replace("-", "_") %>%
    str_replace(fixed("*"), "_")
    
  write.table(sub,
              file = paste(ac_line, ".embl", sep = ""),
              quote = F, col.names = F, row.names = F)
  
}
