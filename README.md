# RCE_project
## Table of content
* Background
* System Requirement
* Example and Usage
* Maintainers and Contributing
* License

## Background

RNA Codon Expansion (RCE) strategy that introduces and decodes bioorthogonally assignable ΨCodons (ΨGA, ΨAA or ΨAG) at specified mRNA transcripts for ncAA incorporation in mammalian cells. This project provided the customized code for assessing sudoU off-target of RCE tools (RNA Codon Expansion)



## System Requirement



### Hardware requirements
This pipeline requires only a standard computer with enough RAM to support the in-memory operations.

### Software requirements
- Linux: Rocky Linux 8.8

### Software 
- PRAISE pipeline (version 1.0)
- Methylsig R package (version 1.17.0)
- Perl (version 5.26)
- Perl packages:
  - Modern::Perl
  - IO::All
- R (version 4.3.1)
- R packages
  - methylSig
  - getopt
  - dplyr
  - bsseq


## Example and Usage

### 0. install related R package and Perl module
install  R package "methylSig":
```bash

# install BiocManager if not exist
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methylSig")
```

install R package getopt

```bash
install.packages("getopt")
```
install R package getopt

```bash
# install BiocManager if not exist
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("bsseq")
```

install Perl module Modern::Perl
```bash
cpanm Modern::Perl
```

install Perl module IO::All
```bash
cpanm IO::All
```

### 1. call transcriptomic-wide Ψ site from PRAISE sequecing method
After obtained raw data of the PRAISE sequencing, we firstly need to call transcriptomic-wide Ψ sites according published PRAISE-pipeline (https://github.com/Zhe-jiang/PRAISE). This step is time-consuming, which usually need several days. Here, we have provided a example file in data dir ("gce_sudoU_sites.xls"). The content of each column in the outfile is as follows:

| Columns | Interpretation |
| :---: | :---: |
| chr_name | transcript name |
| site | sudoU positions in transcripts |
| treated_total_counts | reads coverage at sudoU site in the treated sample |
| treated_deletion_counts | deletion count at sudoU site in the treated sample  |
| treated_deletion_ratio | deletion ratio at sudoU site in the treated sample |
| ctrl_total_counts | reads coverage at sudoU site in the control sample |
| ctrl_deletion_counts |  deletion count at sudoU site in the control sample |
| ctrl_deletion_ratio | deletion ratio at sudoU site in the control sample |
| p_value |  statistically significant whether a true sudoU site |
| chr_site | position in chromosome |

### 2. Transform the format of sudoU_sites file
change the format of sudoU_sites file to conduct statistic analysis, the output file gce_sudoU_sites_for_statistic.xls is used to conduct methylsig analysis, this step usually need less than one minute.
* ``` perl convert_to_methylsig_format.pl  gce_sudoU_sites.xls  gce_sudoU_sites_for_statistic.xls```

### 3. identified offtarget sudoU sites
* ``` Rscript identify_offtarget_site.R --c1 control1_for_statistic.xls --c2 control2_for_statistic.xls --t1 treat1_for_statistic.xls --t2 treat2_for_statistic.xls  -o offtarget.xls ```

The content of each column in the outfile is as follows:
| Columns | Interpretation |
| :---: | :---: |
| gene id | transcript id in Refseq |
| gene symbol | gene symbol |
| gene type | mRNA, lncRNA or other type |
| position | sudoU site position in transcript  |
| control_ratio | sudoU ratio in control sample |
| treat_ratio | sudoU ratio in treatment sample |
| meth_diff |  ratio differenct between control and treatment |
| p_value |  statistically significant of the off-target sudoU site |
| fdr | corrected pvalue |




 
