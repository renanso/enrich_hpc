rm(list=ls())

##################
### GC filter
##################

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("stringr", quietly = TRUE))
  install.packages("stringr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")

dna1<- readDNAStringSet("candidates.fasta", format="fasta",
                        nrec=-1L, skip=0L, seek.first.rec=FALSE,
                        use.names=TRUE, with.qualities=FALSE)
dna1[1]

seq_name = names(dna1)
sequence = paste(dna1)
df <- data.frame(seq_name, sequence)
colnames(df)<-c('qseqid','sequence')

# GC content 

gc_fun <- function(x){ 
  num_g <- str_count(x, "G")
  num_c <- str_count(x, "C")
  return ((num_g + num_c) / str_length(x) * 100 )} 

df$gc<-round(as.numeric(sapply(df[,2],gc_fun)),1)

nrow(df)

##Apply filters
str(df)
df2<- df %>% filter(gc >= 40, gc <= 46)

#gc after
jpeg(file="gc_content.jpeg", width = 1280, height = 1280, units = "px",  quality = 100)
par(mfrow=c(1,2))
hist(df$gc, main = "Original GC %", breaks = 5)
hist(df2$gc, main = "Filtered GC %", breaks = 5)
dev.off()

nrow(df2)
##changing table to fasta format
#add ">" to headers
df2[,1] <- paste0(">",df2[,1])

#bind rows of headers and seqs
probes_fasta <- c(rbind(df2[,1], df2[,2]))
probes_fasta[1:10]

write.table(probes_fasta, "candidates_gc_filtered.fasta", row.names=FALSE,sep="\t", quote = FALSE, col.names = FALSE)

