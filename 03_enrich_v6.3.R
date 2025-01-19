rm(list=ls())

if (!require("data.table", quietly = TRUE))
  install.packages("data.table")
#library(data.table)

#blast_results1<-read.table("results_high_v2.txt")
blast_results1 <- fread("blast_results.txt", showProgress = FALSE)

names<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

colnames(blast_results1)<- names

blast_results1[1:10,1:12]

## Applying filters

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
#library(dplyr)

## Filter 1: % of identity and fragment length
#probes binding in the correct place and with high indentity to off targets will pass
#everything else will be discarded
probes1<- blast_results1 %>% filter(length > 110) %>% filter(pident > 90)

#probes1a<- blast_results1 %>% filter(pident > 80) %>% filter(length > 110)

## Filter 2: Filter out probes that binds to other parts of the genome with more than 90% match

# Return names which have only a single row of data
probes2<-probes1 %>% 
  group_by(qseqid) %>% 
  filter(n()==1)

#Filter 3: Matching pattern to filter probes with specific match to target region

if (!require("stringr", quietly = TRUE))
  install.packages("stringr")
#library(stringr)

##chromosome match test

###stopped here match strings
#original_chr<- str_match(probes2$qseqid, "_S\\s*(.*?)\\_")

original_chr<-(str_split_fixed(probes2$qseqid, "_", 2))[,1]

probes2$original_chr<- original_chr

probes2$target_chr<-probes2$sseqid

probes2$chr_match<- probes2$original_chr==probes2$target_chr

## filter matching results

probes3<-filter(probes2, chr_match == "TRUE")

##position test
##extract original SNP position
#position<-str_match(probes3$qseqid, "S[0-7]{2}_\\s*(.*?)\\|")
position<-(str_split_fixed(probes3$qseqid, '\\|', 2))[,1]
position2<-(str_split_fixed(position, '_',2))[,2]
position2<-as.numeric(position2)
#position2<-as.numeric(position[,2])

#logical test
probes3$position<- position2
probes3$test<-position2 - probes2$sstart
probes3$pos_match<- probes3$test < 300 ##assuming the maximum distance base on the design
##filter position

probes4<-probes3

#drop the initial dataframe to free up memory

rm(blast_results1, probes1, probes2, probes3)
gc()

##################
### Thermodynamics
##################

#prepare data set adding sequence to blast results filtered

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#library(BiocManager)

if (!require("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")

#library("Biostrings")

dna1<- readDNAStringSet("candidates2.fasta", format="fasta",
                        nrec=-1L, skip=0L, seek.first.rec=FALSE,
                        use.names=TRUE, with.qualities=FALSE)
dna1[1]

seq_name = names(dna1)
sequence = paste(dna1)
df <- data.frame(seq_name, sequence)
colnames(df)<-c('qseqid','sequence')

probes5 <- (merge(df, probes4, by = 'qseqid'))

# Install the latest version directly from GitHub

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

if (!require("primer3", quietly = TRUE))
  devtools::install_github("jensenlab/primer3")

##Tm
#iterate on the column (apply function to each row in a column)

probes5$tm<-round(as.numeric(sapply(probes5[,2],calculate_tm)),1)

# GC content 

gc_fun <- function(x){ 
  num_g <- str_count(x, "G")
  num_c <- str_count(x, "C")
  return ((num_g + num_c) / str_length(x) * 100 )} 

probes5$gc<-round(as.numeric(sapply(probes5[,2],gc_fun)),1)

##checking distributions
jpeg(file="TM_gc_1.jpeg")
par(mfrow=c(1,2))
hist(probes5$gc, main = "GC content 1")
hist(probes5$tm, main = "Tm 1")
dev.off()

##Apply filters
str(probes5)
probes6<- probes5 %>% filter(gc >= 40, gc <= 46) %>% filter(tm >= 70, tm <= 76)

#probes6<- probes5 %>% filter(gc >= 42, gc <= 45) %>% filter(tm >= 72, tm <= 75)

##Applying a more stringent filter for GC and Tm
##it could be 71 < tm < 75

jpeg(file="TM_gc_after_filter.jpeg")
par(mfrow=c(1,2))
hist(probes6$gc, main = "GC content - filtered")
hist(probes6$tm, main = "Tm -filtered")
dev.off()

#####################
#secondary structure

##Note that the maximum length of ``seq`` is 60 bp. This is a cap suggested 
#by the Primer3 team as the longest reasonable sequence length for which a 
#two-state NN model produces reliable results.

# We can still analyze the probes by creating fragments up to 60bp.

#frag1 = first 60bp (1 to 60)
#frag2 = last 60bp (61 to 120)
#frag3 = first 30bp + last 30 bp (1 to 30 + 91 to 120)
#frag4 = middle 60bp (31 to 90)
#frag5 = first 30bp + 60-90bp (1 to 30 + 61 to 90)
#frag6 = 30-60 + last 30 bp (31 to 60 + 91 to 120)

probes6$frag1<-substr(probes6$sequence, 1, 60)
probes6$frag2<-substr(probes6$sequence, 61, 120)
probes6$frag3<-paste0(substr(probes6$sequence, 1, 30), substr(probes6$sequence, 91, 120))
probes6$frag4<-substr(probes6$sequence, 31, 90)
probes6$frag5<-paste0(substr(probes6$sequence, 1, 30), substr(probes6$sequence, 61, 90))
probes6$frag6<-paste0(substr(probes6$sequence, 31, 60), substr(probes6$sequence, 91, 120))

##hairpin Tm calculation
probes6_hairpin_frag1<-lapply(probes6$frag1,calculate_hairpin)
probes6_hairpin_frag2<-lapply(probes6$frag2,calculate_hairpin)
probes6_hairpin_frag3<-lapply(probes6$frag3,calculate_hairpin)
probes6_hairpin_frag4<-lapply(probes6$frag4,calculate_hairpin)
probes6_hairpin_frag5<-lapply(probes6$frag5,calculate_hairpin)
probes6_hairpin_frag6<-lapply(probes6$frag6,calculate_hairpin)

probes6$hairpin_temp_frag1<-sapply(probes6_hairpin_frag1,"[[",2)
probes6$hairpin_temp_frag2<-sapply(probes6_hairpin_frag2,"[[",2)
probes6$hairpin_temp_frag3<-sapply(probes6_hairpin_frag3,"[[",2)
probes6$hairpin_temp_frag4<-sapply(probes6_hairpin_frag4,"[[",2)
probes6$hairpin_temp_frag5<-sapply(probes6_hairpin_frag5,"[[",2)
probes6$hairpin_temp_frag6<-sapply(probes6_hairpin_frag6,"[[",2)

#probes6_hairpin[1]
#structure_found<-sapply(probes6_hairpin,"[[",1)
#ds<-sapply(probes6_hairpin,"[[",3)
#dh<-sapply(probes6_hairpin,"[[",4)
#dg<-sapply(probes6_hairpin,"[[",5)
#align_end_1<-sapply(probes6_hairpin,"[[",6)
#align_end_2<-sapply(probes6_hairpin,"[[",7)

probes6_homodimer_frag1<-lapply(probes6$frag1,calculate_homodimer)
probes6_homodimer_frag2<-lapply(probes6$frag2,calculate_homodimer)
probes6_homodimer_frag3<-lapply(probes6$frag3,calculate_homodimer)
probes6_homodimer_frag4<-lapply(probes6$frag4,calculate_homodimer)
probes6_homodimer_frag5<-lapply(probes6$frag5,calculate_homodimer)
probes6_homodimer_frag6<-lapply(probes6$frag6,calculate_homodimer)

probes6$homodimer_temp_frag1<-sapply(probes6_homodimer_frag1,"[[",2)
probes6$homodimer_temp_frag2<-sapply(probes6_homodimer_frag2,"[[",2)
probes6$homodimer_temp_frag3<-sapply(probes6_homodimer_frag3,"[[",2)
probes6$homodimer_temp_frag4<-sapply(probes6_homodimer_frag4,"[[",2)
probes6$homodimer_temp_frag5<-sapply(probes6_homodimer_frag5,"[[",2)
probes6$homodimer_temp_frag6<-sapply(probes6_homodimer_frag6,"[[",2)

## filter for hairpin Tm

probes7<- probes6 %>% filter(hairpin_temp_frag1 < 50) %>% 
  filter(hairpin_temp_frag2 < 50) %>%
  filter(hairpin_temp_frag3 < 50) %>%
  filter(hairpin_temp_frag4 < 50) %>%
  filter(hairpin_temp_frag5 < 50) %>%
  filter(hairpin_temp_frag6 < 50)

## filter for homodimer Tm

probes8<- probes7 %>% filter(homodimer_temp_frag1 < 50) %>% 
  filter(homodimer_temp_frag2 < 50) %>%
  filter(homodimer_temp_frag3 < 50) %>%
  filter(homodimer_temp_frag4 < 50) %>%
  filter(homodimer_temp_frag5 < 50) %>%
  filter(homodimer_temp_frag6 < 50)
  

# Thermodynamics calculations. The most important is the Structure Tm. If the Tm is lower than the reaction temperature,
# there will be no issues with the probe. 

# Other parameters: structure found, ds - change in entropy, dh - change in entalpy, dg - Gibbs free energy


##Filter to keep only one probe per site

#probes8$probe_sites<-data.frame(str_match(probes8[,1], "S\\s*(.*?)\\|"))[,2]
#probes8$probe_type<-str_sub(data.frame(str_match(probes8[,1], "\\|[A-Z]+\\|"))[,1],2,-2)


## sort to benefit highest temperature probes

probes9 <- probes8[order(probes8$tm, decreasing = TRUE),]

## Filter for one probe per gene
#probes9$gene<-str_match(probes9[,1], "[a-z]{6}[0-9]{2}[a-z][0-9]{5}")[,1]
#probes10<- probes9 %>% distinct(gene, .keep_all = TRUE)

probes10<-probes9

saveRDS(probes10, "probes10")
write.table(probes10, "probes10.txt")

###############################
######Visualize probe positions
##plot Positions
str(probes10)

## making the plot

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

##bins every 5M
bins<-seq(from=1, to=0.8e10, by=0.5e8)

##adding the bins to dataframe
probes10$bins<-findInterval(probes10$position, bins)

# extracting the info needed for plots
df2<-probes10[,c(14,17,40)]
df2$bins2<-paste0(df2$original_chr,"-",df2$bins)
df3<-data.frame(table(df2$bins2))
colnames(df3)<-c("bins2","freq")
df4 <- merge(df2,df3,by="bins2")
##keep unique bins
df5<- df4 %>% distinct(bins2, .keep_all = TRUE)

df5$original_chr <- as.numeric(gsub('Chr', '', df5$original_chr))

#plot1<- ggplot(probes10,aes(x=position, y= original_chr, colour = probe_type)) +
#geom_point() +
#scale_y_continuous(breaks = seq(min(probes10$original_chr), max(probes10$original_chr), by = 1)) +
#theme_bw() +
#theme(axis.line = element_line(colour = "black"),
#        panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#      panel.border = element_blank(),
#     panel.background = element_blank()) +
#annotate("segment", x = 1, xend = c(1.84e9,1.39e9,1.1e9,1.1e9,0.9e9,0.75e9,0.76e9), y = 1:7, yend = 1:7,
#        colour = "black") +
#  theme(
#    legend.position = c(0.92, 0.92),
#    legend.justification = c("right", "top"),
#    legend.box.just = "right",
#    legend.margin = margin(6, 6, 6, 6)
#  ) + ggtitle("All probes")

plot2<- ggplot(df5,aes(x=position, y= original_chr)) + 
  geom_count(aes(size=freq)) + 
  scale_y_continuous(breaks = seq(min(df5$original_chr), max(df5$original_chr), by = 1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotate("segment", x = 1, xend = c(1.84e9,1.39e9,1.1e9,1.1e9,0.9e9,0.75e9,0.76e9), y = 1:7, yend = 1:7,
           colour = "black") +
  annotate("text", x = 1.7e9, y = 4.5, label = "Bin= 5M") +
 theme(
  legend.position = c(0.92, 0.92),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
) + ggtitle("All probes")

pdf("plot_all_probes.pdf")
#print(plot1)     # Plot 1 --> in the first page of PDF
print(plot2)     # Plot 2 ---> in the second page of the PDF
dev.off() 


###########################
##Formating the final files
###########################

##ordering by probe type and selecting the high and moderate impact
#probes10_srt<-probes10[order(probes10$probe_type),,drop=FALSE] %>% filter(probe_type != "LOW")

#only reorder
probes10_srt<-probes10[order(probes10$original_chr),]

##Keep only one probe per target
probes10_srt_fltr<-  probes10_srt %>% distinct(position, .keep_all = TRUE)

##probes summary

probes10_sum<- probes10_srt_fltr[,c(1,2,20,21,28:39)]

write.csv(probes10_sum, "final_probes_summary.csv")

##visualizing the probes again

#tm and gc
jpeg(file="final_TM_gc.jpeg", width = 1280, height = 1280, units = "px",  quality = 100)
par(mfrow=c(1,2))
hist(probes10_srt_fltr$gc, main = "GC content", breaks = 5)
hist(probes10_srt_fltr$tm, main = "Tm", breaks = 5)
dev.off()

##distribution of probes across the chromosomes

#plot3<-ggplot(probes10_1k,aes(x=position, y= original_chr, colour = probe_type)) +
#  geom_point() +
#  scale_y_continuous(breaks = seq(min(probes10$original_chr), max(probes10$original_chr), by = 1)) +
#  theme_bw() +
#  theme(axis.line = element_line(colour = "black"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank()) +
#  annotate("segment", x = 1, xend = c(1.84e9,1.39e9,1.1e9,1.1e9,0.9e9,0.75e9,0.76e9), y = 1:7, yend = 1:7,
#           colour = "black") +
#  theme(
#    legend.position = c(0.92, 0.92),
#    legend.justification = c("right", "top"),
#    legend.box.just = "right",
#    legend.margin = margin(6, 6, 6, 6)
#  ) + ggtitle("1k probes")

# extracting the info needed for plot 2
#df2<-probes10_1k[,c(14,17,37)]
#df2$bins2<-paste0(df2$original_chr,"-",df2$bins)
#df3<-data.frame(table(df2$bins2))
#colnames(df3)<-c("bins2","freq")
#df4 <- merge(df2,df3,by="bins2")
##keep unique bins
#df5<- df4 %>% distinct(bins2, .keep_all = TRUE)

#plot4<-ggplot(df5,aes(x=position, y= original_chr)) + 
#  geom_count(aes(size=freq)) + 
#  scale_y_continuous(breaks = seq(min(probes10$original_chr), max(probes10$original_chr), by = 1)) +
#  theme_bw() +
#  theme(axis.line = element_line(colour = "black"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank()) +
#  annotate("segment", x = 1, xend = c(1.84e9,1.39e9,1.1e9,1.1e9,0.9e9,0.75e9,0.76e9), y = 1:7, yend = 1:7,
#           colour = "black") +
#  annotate("text", x = 1.7e9, y = 4.5, label = "Bin= 5M") +
#  theme(
#    legend.position = c(0.92, 0.92),
#    legend.justification = c("right", "top"),
#    legend.box.just = "right",
#    legend.margin = margin(6, 6, 6, 6)
#  )+ ggtitle("1k probes")

#pdf("plot_1k_probes.pdf")
#print(plot3)     # Plot 3 --> in the first page of PDF
#print(plot4)     # Plot 4 ---> in the second page of the PDF
#dev.off() 

names<-(str_split_fixed(probes10_srt_fltr$qseqid, '\\|',6))[,c(1:5)]
final_probes<- cbind(paste0(names[,1],"-",names[,4],"-",names[,5],"-Tm=", probes10_srt_fltr$tm,"-GC=", probes10_srt_fltr$gc), probes10_srt_fltr$sequence)

#final_probes<- cbind(paste0(probes10_srt_fltr$qseqid,"|position=",probes10_srt_fltr$sstart,"-", probes10_srt_fltr$send,"|Tm=", probes10_srt_fltr$tm,"|GC=", probes10_srt_fltr$gc), probes10_srt_fltr$sequence)

##changing table to fasta format
#add ">" to headers
final_probes[,1] <- paste0(">",final_probes[,1])

#bind rows of headers ans seqs
probes_fasta <- c(rbind(final_probes[,1], final_probes[,2]))
probes_fasta[1:10]

write.table(probes_fasta, "final_probes.fasta", row.names=FALSE,sep="\t", quote = FALSE, col.names = FALSE)

