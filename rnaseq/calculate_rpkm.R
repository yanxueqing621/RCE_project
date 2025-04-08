
library(getopt)
library(DESeq2)

spec = matrix(c(
                        'infile', 'i', 1, 'character',
            'control','c', 1, 'character',
            'treat','t', 1, 'character',
                        'help',  'h', 0, 'logical'
                        ), byrow = TRUE, ncol = 4)

args = getopt(spec)
if (!is.null(args$help) || is.null(args$infile)) {
    #gid length dmso1 dmso2  stm1 stm2
        cat("Rscript rna_feature_merge_rpkm_r -i all_sample_count.xls -c control_sample_name -t treat_sample_name \n")
        cat(paste(getopt(spec=spec, usage = T), "\n"))
        quit()
}

fpkm_formula <- function(df, len) {
  head(len)
  fpkm = function(x,len) {
    N = sum(x)
    exp(log(x) + log(1e9) - log(len) - log(N))
  }
  df.fpkm = apply(df, 2, fpkm, len)
  return(df.fpkm)
}


calc_fpkm_tpm <- function(df){
  df.fpkm = fpkm_formula(df[,2:ncol(df)], df$Length)
  df.fpkm = round(df.fpkm, 2)
  colnames(df.fpkm) = paste(colnames(df.fpkm), "_FPKM", sep = "")
  df.fpkm2 = data.frame(length = df$Length, df.fpkm)
  fpkm_tpm =  df.fpkm2 
  if (!file.exists("all_fpkm.txt")) {
    write.table(fpkm_tpm, file = "all_fpkm_tpm.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  }
  return(df.fpkm)
}


df = read.table(args$infile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
print(head(df))
fpkm = calc_fpkm_tpm(df)

