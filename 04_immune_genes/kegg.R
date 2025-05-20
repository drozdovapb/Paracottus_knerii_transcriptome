library(ggplot2)
library(openxlsx)

kegg <- read.xlsx("KEGG 5-1 immune system.xlsx")
kegg <- kegg[!(kegg$ID %in% c("map04625", "map04613")), ] #these did not exist at the moment eggNOG5.0 was created
kegg$ID

eggnog <- read.xlsx("../03_annotation/Pkn_rnaspades_ssrf_out.emapper.annotations.fromChordata.xlsx", startRow = 3)

immune_rows <- unlist(sapply(kegg$ID, function(k){c(grep(pattern = k, x = eggnog$KEGG_Pathway))}))
eggnog_immune <- eggnog[immune_rows, ]
write.xlsx(eggnog_immune, "../04_immune_genes/eggnog_immune_genes.xlsx")

kegg$ntranscripts <- sapply(kegg$ID, function(k){
  sum(grepl(pattern = k, x = eggnog$KEGG_Pathway))})

ggplot(kegg, aes(x=Name, y=ntranscripts)) + 
  geom_bar(stat = 'identity', fill='#00ff0080', col='#003300ff') + 
  coord_flip() + 
  xlab('KEGG pathway') + ylab('Number of annotated transcripts') + 
  theme_bw() + 
  theme(axis.text = element_text(size=10)) + #, axis.title = element_text(size=16)) + 
  scale_x_discrete(limits=rev)

ggsave('kegg_transcripts.png', device = png, width = 8, height = 5)
