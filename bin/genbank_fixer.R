#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## THIS SCRIPT IS MEANT TO REMOVE ALL TRACES OF GENEIOUS FROM THE GENBANK FILE
## BEFORE IT PROCEEDS INTO PARSING WITH BIOPYTHON. GENEIOUS USES CONVENTIONS
## THAT USUALLY LEAD TO ERRORS WITH BIOPYTHON'S STRICT PARSING, SO WE SOLVE THAT
## BY REMOVING DISCREPANCIES FIRST.

# loading tidyverse
library(tidyverse)

# importing genbank file
gbk <- read.delim(args[1], header = F)

# # removing Geneious's odd local version information, which contains characters that confuse parsers
# gbk <- data.frame(V1 = str_remove_all(gbk$V1, fixed("geneious|urn:local:")))
# 
# # recording all the lines where the version, accession, and locus features occur
# version_lines <- which(grepl("VERSION", gbk$V1))
# accession_lines <- which(grepl("ACCESSION", gbk$V1))
locus_lines <- which(grepl("LOCUS", gbk$V1))

# removing Geneious notes
no_geneious_notes <- which( !grepl("Derived using Geneious Prime 2022.0.2", gbk$V1))
gbk <- data.frame(V1 = gbk[no_geneious_notes, ])
no_geneious_notes <- which(!grepl("from Database' based on nucleotide similarity", gbk$V1))
gbk <- data.frame(V1 = gbk[no_geneious_notes, ])
no_geneious_notes <- which( !grepl("/Source=Geneious", gbk$V1))
gbk <- data.frame(V1 = gbk[no_geneious_notes, ])
no_geneious_notes <- which( !grepl("/modified_by=", gbk$V1))
gbk <- data.frame(V1 = gbk[no_geneious_notes, ])

# 
# # renaming version and accession lines
# for (i in 1:length(version_lines)){
#   
#   putative <- unlist(str_split(gbk$V1[locus_lines[i]], "       "))[2]
#   locus <- unlist(str_split(putative, fixed("|")))[1]
#   gbk[version_lines[i], "V1"] <- paste("VERSION     ", locus, sep = "")
#   gbk[accession_lines[i], "V1"] <- paste("ACCESSION   ", locus, sep = "")
#   
# }

# Making sure the locus name is also in the definition
deflines <- which(grepl("DEFINITION", gbk$V1))
for (i in 1:length(deflines)){
  
  definition <- gbk$V1[deflines[i]]
  
  if ( grepl("|", definition, fixed = T) ){
    next
  } else {
    
    locus <- gbk$V1[locus_lines[i]] %>%
      str_remove(pattern = "LOCUS       ") %>%
      str_split(fixed("|")) %>%
      unlist()
    locus <- paste(locus[1], "|", sep = "")
    
    definition <- str_replace(definition, "DEFINITION  ", 
                              paste("DEFINITION  ", locus, sep = ""))
    
    gbk$V1[deflines[i]] <- definition
    
  }
  
}

# Write a new new GenBank file
write.table(gbk, paste(str_remove(args[1], ".gb"),
                       "_corrected.gb", sep = ""),
            quote = F, row.names = F, col.names = F)
