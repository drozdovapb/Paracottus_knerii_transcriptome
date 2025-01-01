## packages
library(ggplot2)
library(openxlsx)
library(dplyr)

## get names
## the whole list
cleaned_by_NCBI <- readLines("./Pkn_rnaspades_ssrf_FCS_GX_cleaned.names.txt")

## eggNOG annotation from Chordata
Chordata_eggNOG_tbl <- read.xlsx("./out.emapper.annotations.fromChordata.xlsx", startRow = 3)
Chordata_eggNOG_tbl$query <- sapply(Chordata_eggNOG_tbl$query, function(x) {substr(x, 1, nchar(x)-3)})
Chordata_eggNOG <- Chordata_eggNOG_tbl$query

## TRAPID
trapid <-read.delim("./transcripts_tax_exp7497.txt")
## Chordata from TRAPID
Chordata_TRAPID <- trapid[grep(pattern = "Chordata", trapid$lineage), ]$transcript_id

## merge (intersect) Chordata annotation, contamination-free
chordata <- c(Chordata_eggNOG, Chordata_TRAPID)
chordata <- intersect(chordata, cleaned_by_NCBI)
#writeLines(chordata, "fish_transcript_names.txt")

## assemble name + annotation table
kv <- data.frame(query = chordata)
kv2 <- left_join(x=kv, y=Chordata_eggNOG_tbl[ , c("query", "Preferred_name", "Description")]) 
kv2$Preferred_name[is.na(kv2$Preferred_name)] <- ""
kv2$Description[is.na(kv2$Description)] <- ""

## make the transcript names shorter (NCBI only allows 50 symbols)
kv2$short_query <- gsub("\\..*", "", kv2$query)
kv2$short_query <- gsub("length", "len", kv2$short_query)

## get rid of non-accepted symbols
kv2$Annotation <- paste(kv2$Preferred_name, kv2$Description, sep="_")
kv2$Annotation <- gsub("\\,", "", kv2$Annotation)
kv2$Annotation <- gsub("\\[", "", kv2$Annotation)
kv2$Annotation <- gsub("\\]", "", kv2$Annotation)
kv2$Annotation <- gsub("\\'", "", kv2$Annotation)

## combine names and annotations 
kv2$Name <- paste(kv2$short_query, kv2$Annotation, sep = "_")
## and shorten the results
kv2$Name <- substr(kv2$Name, 1, 50)

## write out transcripts for filtering by name with seqkit
writeLines(kv2$query, "transcripts_for_TSA_submission_names.txt")

## write out table for seqkit -kv
write.table(kv2[, c("query", "Name")], "annotated_transcripts.tsv",  
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")