library(ggplot2)
library(openxlsx)

## code for Fig. 1
Pkn_filtlenN <- readLines("./Pkn_rnaspades_ssrf_filtlength.names.txt") #120k
Chordata_eggNOG <- read.xlsx("Pkn_rnaspades_ssrf_out.emapper.annotations.fromChordata.xlsx", startRow = 3)
Chordata_eggNOG <- Chordata_eggNOG$query
Chordata_eggNOG <- sapply(Chordata_eggNOG, function(x) {substr(x, 1, nchar(x)-3)})
Chordata_eggNOG <- intersect(Chordata_eggNOG, Pkn_filtlenN)

trapid <-read.delim("./transcripts_tax_exp7497.txt")
Chordata_TRAPID <- trapid[grep(pattern = "Chordata", trapid$lineage), ]$transcript_id

nVennR::plotVenn(list(all=Pkn_filtlenN, fish_eggnog=Chordata_eggNOG, fish_trapid=Chordata_TRAPID), 
                   outFile="fishos.svg", fontScale=1.5)

#library(ggvenn)
#ggvenn(data = list(all=Pkn_filtlenN, eg=Chordata_eggNOG, tra=Chordata_TRAPID))

library(ggVennDiagram)
ggVennDiagram(x=list(all=Pkn_filtlenN, fish_eggnog=Chordata_eggNOG, fish_trapid=Chordata_TRAPID), label_percent_digit = FALSE)

library(BioVenn)
svg("euler_fishes.svg")
BioVenn::draw.venn(list_x=Pkn_filtlenN, list_y=Chordata_eggNOG, list_z=Chordata_TRAPID, 
                   title = "", subtitle = "",
                   xtitle = "All transcripts \n\n\n", xt_c = "black", x_c = "grey90",
                   ytitle = "Chordata / eggNOG", yt_c = "brown4", y_c = "lightgoldenrod",
                   ztitle = "\n Chordata / TRAPID", zt_c="royalblue", z_c="skyblue")
dev.off()

## TRAPID
trapid_groups <- read.xlsx("TRAPID_selected_groups.xlsx")
trapid_groups$Taxon <- factor(trapid_groups$Taxon, 
                              levels = c( "Chordata", "Bacteria", "Platyhelminthes", "Arthropoda", 
                                         "Fungi", "Viridiplantae", "Ciliophora", "Archaea", "Viruses",
                                         "Other", "Unclassified"))

ggplot(trapid_groups, aes(x=Taxon, y=Nseqs)) + 
  geom_bar(stat='identity', col='#0000ff', fill='#85cceb80') + 
  coord_flip() +
  xlab('Taxon') + ylab('Number of annotated transcripts') + 
  theme_bw(base_size = 14) + theme(axis.text = element_text(size=14), axis.title = element_text(size=16)) + 
  scale_x_discrete(limits=rev)
ggsave('kegg_transcripts.svg', device = svg, width = 8, height = 5)
ggsave('trapid_transcripts.png', device = png, width = 8, height = 5)


#### code not used in the manuscript
## full transcript list
Pkn_all <- readLines("./Pkn_rnaspades_ssrf.names.txt") #192k
Pkn_filtlenN <- readLines("./Pkn_rnaspades_ssrf_filtlength.names.txt") #120k

cleaned_by_NCBI <- readLines("./Pkn_rnaspades_ssrf_filtgfx.names.txt")
filtout_by_NCBI <- readLines("./Pkn_rnaspades_ssrf_filteredOutNCBI.names.txt")


filtout_by_NCBI_on_filtlenN <- intersect(filtout_by_NCBI, Pkn_filtlenN) #even 817, it's too few

Chordata_eggNOG <- read.xlsx("./out.emapper.annotations.fromChordata.xlsx", startRow = 3)
Chordata_eggNOG <- Chordata_eggNOG$query
Chordata_eggNOG <- sapply(Chordata_eggNOG, function(x) {substr(x, 1, nchar(x)-3)})
Chordata_eggNOG <- intersect(Chordata_eggNOG, Pkn_filtlenN)

trapid <-read.delim("./transcripts_tax_exp7497.txt")
Chordata_TRAPID <- trapid[grep(pattern = "Chordata", trapid$lineage), ]$transcript_id

Bacteria_TRAPID <- trapid[grep(pattern = "Bacteria", trapid$lineage), ]$transcript_id
Archaea_TRAPID <- trapid[grep(pattern = "Archaea", trapid$lineage), ]$transcript_id
Worms_TRAPID <- trapid[grep(pattern = "Platyhelminthes", trapid$lineage), ]$transcript_id
Fungi_TRAPID <- trapid[grep(pattern = "Fungi", trapid$lineage), ]$transcript_id
Viruses_TRAPID <- trapid[grep(pattern = "Viruses", trapid$lineage), ]$transcript_id
Arthropoda_TRAPID <- trapid[grep(pattern = "Arthropoda", trapid$lineage), ]$transcript_id
Ciliata_TRAPID <- trapid[grep(pattern = "Ciliophora", trapid$lineage), ]$transcript_id
Plantae_TRAPID <- trapid[grep(pattern = "Viridiplantae", trapid$lineage), ]$transcript_id

library(UpSetR)
lst <- list(all=Pkn_filtlenN, NCBI_out=filtout_by_NCBI_on_filtlenN, 
            chordata_ann=Chordata_eggNOG, chordata_tax=Chordata_TRAPID,
            bact = Bacteria_TRAPID, vir = Viruses_TRAPID)
upset(fromList(lst))


venneuler::venneuler()

library(ggvenn)
ggvenn(data=lst)

#devtools::install_github("vqf/nVennR")
library(nVennR)
myV <- plotVenn(lst, outFile="1.png")

myV <- plotVenn(list(all=Pkn_filtlenN, 
                     bact=Bacteria_TRAPID, vir=Viruses_TRAPID, worms=Worms_TRAPID, arch=Archaea_TRAPID, fungi=Fungi_TRAPID), 
                outFile="ex.pdf")


myV <- plotVenn(list(all=Pkn_filtlenN, 
                     bact=Bacteria_TRAPID, vir=Viruses_TRAPID, worms=Worms_TRAPID, crust = Arthropoda_TRAPID, exNCBI = filtout_by_NCBI_on_filtlenN), 
                outFile="ex.pdf")

excludeVenn <- plotVenn(list(bact=Bacteria_TRAPID, vir=Viruses_TRAPID, worms=Worms_TRAPID, crust = Arthropoda_TRAPID, exNCBI = filtout_by_NCBI_on_filtlenN), 
                outFile="ex.svg", fontScale=1.5)


excludeVenn <- plotVenn(list(Bacteria=Bacteria_TRAPID, Viruses=Viruses_TRAPID, Ciliata=Ciliata_TRAPID, Plants = Plantae_TRAPID,
                             Worms=Worms_TRAPID, Crustacea = Arthropoda_TRAPID, exNCBI = filtout_by_NCBI_on_filtlenN), 
                        outFile="ex.svg", fontScale=1.5, nCycles = 10^5)

excludeVenn <- plotVenn(list(bact=Bacteria_TRAPID, vir=Viruses_TRAPID, worms=Worms_TRAPID, 
                             crust = Arthropoda_TRAPID, exNCBI = filtout_by_NCBI_on_filtlenN, fish=Chordata_eggNOG), 
                        outFile="ex2.svg", fontScale=1.5)


excludeTranscripts <- c(Bacteria_TRAPID, Viruses_TRAPID, Worms_TRAPID, Arthropoda_TRAPID, filtout_by_NCBI_on_filtlenN)
parasites <- c(Bacteria_TRAPID, Viruses_TRAPID, Worms_TRAPID, Arthropoda_TRAPID)

clearedVenn <- plotVenn(list(all=Pkn_filtlenN, toExclude = excludeTranscripts), 
                        outFile="allminusEx.svg", fontScale=1.5)

alloperationsVenn <- plotVenn(list(all=Pkn_filtlenN, toExclude = excludeTranscripts, fish=Chordata_eggNOG, fish2=Chordata_TRAPID), 
                        outFile="allPlusMinus.svg", fontScale=1.5)

alloperationsVenn <- plotVenn(list(all=Pkn_filtlenN, toExclude = excludeTranscripts, fish2=Chordata_TRAPID), 
                              outFile="allPlusMinus.svg", fontScale=1.5)

fishos <- plotVenn(list(all=Pkn_filtlenN, fish_eggnog=Chordata_eggNOG, fish_trapid=Chordata_TRAPID), 
                   outFile="fishos.svg", fontScale=1.5)

excludeVenn$set


setdiff(filtout_by_NCBI_on_filtlenN, parasites)

