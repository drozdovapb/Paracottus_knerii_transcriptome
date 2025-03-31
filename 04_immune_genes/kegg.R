library(ggplot2)
library(openxlsx)

kegg <- read.xlsx("KEGG 5-1 immune system.xlsx")
kegg$ID

eggnog <- read.xlsx("../03_annotation/Pkn_rnaspades_ssrf_out.emapper.annotations.fromChordata.xlsx", startRow = 3)


kegg$ntranscripts <- sapply(kegg$ID, function(k){
  sum(grepl(pattern = k, x = eggnog$KEGG_Pathway))})

ggplot(kegg, aes(x=Name, y=ntranscripts)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  xlab('KEGG pathway') + ylab('Number of annotated transcripts') + 
  theme_bw() + 
  scale_x_discrete(limits=rev)

ggsave('kegg_transcripts.png', device = png, width = 8, height = 5)
