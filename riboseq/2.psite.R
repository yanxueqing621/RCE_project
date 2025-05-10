library(riboWaltz)
library(ggplot2)
library(cowplot)
library(data.table)
library(magrittr)
library(patchwork)
library(ggsignif)
library(reshape2)
library("ggsci")
library(RColorBrewer)
library(eulerr)

df.anno = fread("gencode.v44.annotation_longest_cds.gtf.anno.first_inframe_3TC.anno", header = TRUE, sep="\t", stringsAsFactors = TRUE)
name_of_bams <- c('test')   #sample names used in variable
names(name_of_bams) <- name_of_bams  # bam file prefix

reads_list <- bamtolist(bamfolder = "./", annotation = df.anno, name_samples = name_of_bams)
reads_list_filter <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 25:32)
psite_offset=psite(reads_list_filter, flanking = 6, extremity = "5end")
reads_psite_list <- psite_info(reads_list_filter, psite_offset)
write.table(reads_psite_list[['test']][,c('transcript', 'psite')], "test_psite_25_32.xls",  row.names = FALSE,col.names = TRUE, quote = FALSE)


