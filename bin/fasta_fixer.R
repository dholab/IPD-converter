#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Make sure tidyverse functions are available.
library(tidyverse)

# Read in fasta file, making sure to specify that there are no column names
fasta = read.delim(args[1], header = F)

# Finding which line numbers are deflines, starting with ">"
deflines <- which(grepl(">", fasta$V1))

# Looping through each defline, and removing the pipe symbol and any text after it
for (i in deflines){
  
  defline <- fasta$V1[i]
  fasta$V1[i] <- unlist(str_split(defline, fixed("|")))[1]
  
}

# Exporting the edited FASTA
write.table(fasta, paste(str_remove(args[1], ".fasta"),
                         "_short_names.fasta", sep = ""),
            quote = F, row.names = F, col.names = F)
