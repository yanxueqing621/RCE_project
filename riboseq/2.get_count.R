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
#df.anno = fread("gencode.v44.annotation_longest_cds.gtf.anno.TGA", header = TRUE, sep="\t", stringsAsFactors = TRUE) 


############################# quality control  QC ################################
##### 1. 读取bam文件 

GFP = 'GFP12'
RCE = 'RCE12'
WT = 'WT12'

df.anno = fread("./gencode.v44.annotation_longest_cds.gtf.anno", header = TRUE, sep="\t", stringsAsFactors = TRUE) 
name_of_bams <- c(GFP, RCE, WT)   #sample names used in variable
names(name_of_bams) <- c(GFP, RCE, WT)  # bam file prefix
reads_list <- bamtolist(bamfolder = "./", annotation = df.anno, name_samples = name_of_bams)


##### 2. ribo-seq mapping reads length distribution
plot_len_dis  <- function(reads_list_tmp, sample_name, min_len= 20, max_len = 40) {
  rlen <- rlength_distr(reads_list_tmp, sample = sample_name, cl = 100)
  df_len = rlen$dt  
  print(names(df_len))
  sample = unlist(strsplit(names(df_len)[2], "_"))[1]
  names(df_len) = c("len", "count", 'percent')
  df_len = df_len[df_len$len>=min_len & df_len$len<=max_len, ]
  p = ggplot(data = df_len, aes(x = len, y = percent)) + geom_bar(stat = 'identity')+ theme_bw()
  p = p + labs(title = sample, x = "Read length", y = "Count (%)")
  p = p + theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
  #p
  return(p)
}
mgta_len_plot = plot_len_dis(reads_list, GFP)
rce_len_plot = plot_len_dis(reads_list, RCE)
wt_len_plot = plot_len_dis(reads_list, WT)
p.all = mgta_len_plot + rce_len_plot + wt_len_plot
p.all

##### 3. psite analyses  选取长度为25:32的进行分析
reads_list_filter <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 25:32)
psite_offset=psite(reads_list_filter, flanking = 6, extremity = "5end")
reads_psite_list <- psite_info(reads_list_filter, psite_offset)


######## 4. plot  psite in three annotated transcript region(5'UTR, CDS, 3'UTR)
wt_psite_region <- region_psite(reads_psite_list, df.anno, sample = WT)
rce_psite_region <- region_psite(reads_psite_list, df.anno, sample = RCE)
mgta_psite_region <- region_psite(reads_psite_list, df.anno, sample = GFP)
plot_grid(wt_psite_region$plot, rce_psite_region$plot,mgta_psite_region$plot, nrow = 1, align = "h")


######## 5. plot  psite  trinucleotide periodicity of all read length  
TCAWT_periodicity <- frame_psite_length(reads_psite_list, sample = WT, region = "all", cl = 90)
RCETGA_periodicity <- frame_psite_length(reads_psite_list, sample = RCE, region = "all", cl = 90)
GFP_periodicity <- frame_psite_length(reads_psite_list, sample = GFP, region = "all", cl = 90)
plot_grid(TCAWT_periodicity$plot, RCETGA_periodicity$plot,GFP_periodicity$plot, nrow = 1, align = "h")


######## 6. plot  Psite in frame of CDS
TCAWT_frame <- frame_psite(reads_psite_list, sample = WT, region = "all")
RCETGA_frame <- frame_psite(reads_psite_list, sample = RCE, region = "all")
GFP_frame <- frame_psite(reads_psite_list, sample = GFP, region = "all")
plot_grid(TCAWT_frame$plot, RCETGA_frame$plot,GFP_frame$plot, nrow = 1, align = "h")



######## 7. Metaplot

TCAWT_metaprofile1 <- metaprofile_psite(reads_psite_list, df.anno, sample = WT,
                                          utr5l = 100, cdsl = 100, utr3l =100,
                                          plot_title = "sample.transcript")
TCAWT_metaprofile1
RCETGA_metaprofile2 <- metaprofile_psite(reads_psite_list, df.anno, sample = RCE,
                                          utr5l = 100, cdsl = 100, utr3l =100,
                                          plot_title = "sample.transcript")
RCETGA_metaprofile2
MTGA_metaprofile3 <- metaprofile_psite(reads_psite_list, df.anno, sample =  GFP ,
                                          utr5l = 100, cdsl = 100, utr3l =100,
                                          plot_title = "sample.transcript")
MTGA_metaprofile3
sample_list <- list( WT  = c( WT ),  RCE  = c(RCE), GFP = c(GFP))
metaprofile_overlaid <- metaprofile_psite(reads_psite_list, df.anno, sample = sample_list,
                                                  multisamples = "average", plot_style = "overlaid",
                                                  utr5l = 100, cdsl = 3, utr3l = 100,length_range = 29,
                                                  frequency = TRUE,  plot_title = "transcript.length_range",
                                                  colour = c("gray70","aquamarine4", 'red'))
metaprofile_overlaid[["plot"]]







#################################  定制化分析  #########################################

########### 1. 样本整体通读率  3'UTR / CDS   只对TGA基因进行分析  
df.anno = fread("gencode.v44.annotation_longest_cds.gtf.anno.TGA", header = TRUE, sep="\t", stringsAsFactors = TRUE) 
name_of_bams <- c(GFP, RCE, WT)   #sample names used in variable
names(name_of_bams) <- c(GFP, RCE, WT)  # bam file prefix
reads_list <- bamtolist(bamfolder = "./", annotation = df.anno, name_samples = name_of_bams)
reads_list_filter <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 25:32)
psite_offset=psite(reads_list_filter, flanking = 6, extremity = "5end")
reads_psite_list <- psite_info(reads_list_filter, psite_offset)

pos_count = function(sample) {
  mgta = reads_psite_list[[sample]]
  print(sample)
  cds.total = sum(mgta$psite_region == 'cds')
  total = nrow(mgta)
  
  mgta = mgta[mgta$psite_from_stop <= 100 & mgta$psite_from_stop >= -100, ]
  total_cds = sum(mgta$psite_region == 'cds')
  total_utr3 = sum(mgta$psite_region == '3utr'& mgta$psite_from_stop >= 5)
  
  print(paste(total_utr3, total_cds,  total_utr3/total_cds, total_utr3/total, sep= "\t"))
  mgta$transcript = as.character(mgta$transcript)
  mgta_stop = mgta[mgta$psite_from_stop <= 100 & mgta$psite_from_stop >=-100, ]
  mgta_stop_count = as.data.frame(table(mgta_stop$psite_from_stop),stringsAsFactors = FALSE)
  names(mgta_stop_count) = c("pos", "count")
  mgta_stop_count$pos = as.integer(mgta_stop_count$pos)
  mgta_stop_count$samples = sample
  row.names(mgta_stop_count) = mgta_stop_count$pos
  #mgta_stop_count$sum = sum(mgta_stop_count[-100:-1, 'count'])
  mgta_stop_count$norm = mgta_stop_count$count / total * 100
  return(list(count = mgta_stop_count, ratio = total_utr3/total_cds))
}

tcawt.result = pos_count(WT)
mtga.result = pos_count( GFP )
rectga.result = pos_count(RCE)


tcawt.count = tcawt.result$count
mtga.count = mtga.result$count
rectga.count = rectga.result$count
tcawt.ratio = tcawt.result$ratio
mtga.ratio = mtga.result$ratio
rectga.ratio = rectga.result$ratio


############ 2. metaplot ############################# 
tcawt_mtga_rectga = rbind(tcawt.count, mtga.count, rectga.count)
tcawt_mtga_rectga$samples = factor(tcawt_mtga_rectga$samples, levels = c(GFP, RCE, WT) )
tcawt_mtga_rectga2 = tcawt_mtga_rectga[tcawt_mtga_rectga$pos >0, ]


colors = c(rgb(223,111, 106, maxColorValue = 255), rgb(232, 199, 100, maxColorValue = 255), "#a7c957")
p1 = ggplot(tcawt_mtga_rectga, aes(x = pos, y = norm))  
p1 = p1 + geom_line(aes(color = samples), linewidth = 0.5)
p1 = p1 + theme_classic() + xlab("Position relative to stop codon")+ theme(legend.position = 'top') + ylab('Percent')
p1 = p1 + scale_color_manual(values = colors)+theme(text = element_text(size = 8))
p1

p2 = ggplot(tcawt_mtga_rectga2, aes(x = pos, y = norm))  
p2 = p2 + geom_line(aes(color = samples), linewidth = 0.6)
p2 = p2 + theme_classic() + xlab("Position relative to stop codon")+ theme(legend.position = 'top')+ ylab('Percent')
#p2 = p2 +  scale_color_npg(alpha = 0.9) 
p2 = p2 + scale_color_manual(values = colors) +theme(text = element_text(size = 8)) + theme(text = element_text(size = 15))

p2
 
ggsave("cds100_3utr100.TGA.png",p1, width = 12, height = 10, units = 'cm', dpi = 300)
ggsave("cds0_3utr100.TGA.png",p2, width = 12, height = 10, units = 'cm',dpi = 300)
 
#colors2 = rev(colors)
df = data.frame(type = factor(c(GFP, RCE, WT),levels =  c(GFP, RCE, WT)) ,  value =c(mtga.ratio, rectga.ratio, tcawt.ratio) ) 
p = ggplot(df, aes(x = type, y = value, fill = type))+geom_bar(stat ='identity', width = 0.72) +theme_classic() + ylab("3'UTR/CDS") +xlab("") +  scale_fill_manual(values = colors)
p = p + theme(legend.position = "none") +theme(text = element_text(size = 7))  
p
p23 = ggdraw(p2 + draw_plot(p, 50, 0.0036, 50, 0.01))
p23
#ggsave("cds0_3utr100.TGA_report_readthough.TGA.pdf", dpi = 300, scale = 1 )
ggsave("cds0_3utr100.TGA_report_readthough.TGA.png", width = 10, heigh = 8, unit = 'cm')


####################### 4. RRTS得分计算 ########################

####### 4.1 make data.frame only include transcript id and psite position
df.anno = fread("gencode.v44.annotation_longest_cds.gtf.anno", header = TRUE, sep="\t", stringsAsFactors = TRUE) 
name_of_bams <- c(GFP, RCE, WT)   #sample names used in variable
names(name_of_bams) <- c(GFP, RCE, WT)  # bam file prefix
reads_list <- bamtolist(bamfolder = "./", annotation = df.anno, name_samples = name_of_bams)
reads_list_filter <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 25:32)
psite_offset=psite(reads_list_filter, flanking = 6, extremity = "5end")
reads_psite_list <- psite_info(reads_list_filter, psite_offset)
write.table(reads_psite_list[[GFP]][,c('transcript', 'psite')], "GFP12_psite_25_32.xls",  row.names = FALSE,col.names = TRUE, quote = FALSE)
write.table(reads_psite_list[[RCE]][,c('transcript', 'psite')], "RCE12_psite_25_32.xls",  row.names = FALSE,col.names = TRUE, quote = FALSE)
write.table(reads_psite_list[[WT]][,c('transcript', 'psite')], "WT12_psite_25_32.xls",  row.names = FALSE,col.names = TRUE, quote = FALSE)

