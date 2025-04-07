library("methylSig")
library(getopt)
library(ggplot2)
library(dplyr)
library(bsseq)
spec = matrix(c(
                        'c1',    '1', 1, 'character',
                        'c2',    '2', 1, 'character',
                        't1',      '3', 1, 'character',
                        't2',      '4', 1, 'character',
                        'outfile','o', 1, 'character',
                        'help',  'h', 0, 'logical'
                        ), byrow = TRUE, ncol = 4)
args = getopt(spec)
if (!is.null(args$help) || is.null(args$n1)) {
        cat("Rscript identify_offtarget_site.R --c1 control1_for_statistic.xls --c2 control2_for_statistic.xls  -o offtarget.xls\n")
        cat(paste(getopt(spec=spec, usage = T), "\n"))
        quit()
}

files1 = c(
  args$c1,
  args$c2,
  args$t1,
  args$t2
)
bsseq1 = bsseq::read.bismark(
  files = files1,
  colData = data.frame(type=c('hela','hela','hypoxia', 'hypoxia')),
  rmZeroCov = FALSE,
  strandCollapse = FALSE
)
bsseq2 = filter_loci_by_group_coverage(
  bs = bsseq1,
  group_column = 'type',
  c('hela' = 2, 'hypoxia' = 2))

diff_gra = diff_methylsig(
  bs = bsseq2,
  group_column = 'type',
  comparison_groups = c('control' = 'hela', 'case' = 'hypoxia'),
  disp_groups = c('case' = TRUE, 'control' = TRUE),
  local_window_size = 0,
  t_approx = TRUE,
  n_cores = 1)


df.offtarget = diff_gra %>% filter(
  fdr < 0.05,
  meth_control < 5,
  meth_case > 10,
  meth_case - meth_control >10,
)

write.table(df.offtarget, file = "RCE_offtarget_dif10.xls", sep = "\t", quote = FALSE, row.names = FALSE)

