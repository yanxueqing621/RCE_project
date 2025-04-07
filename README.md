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
- 



## Example and Usage
### 一. call transcriptomic-wide Ψ site from PRAISE sequecing method
After obtained raw data of the PRAISE sequencing, we firstly need to call transcriptomic-wide Ψ sites according published PRAISE-pipeline (https://github.com/Zhe-jiang/PRAISE). This step is time-consuming, which usually need several days. Here, we have provided a example file in data dir ("sudoU_sites.xls"). The content of each column in the file is as follows:

| Columns | Interpretation |
| :---: | :---: |
| Chr | chromosome |
| Sites | genomic loci |
| Strand | strand |
| Gene | annotated gene |
| CR | conversion rate for genes |
| AGcov | reads coverage with A and G |
| Acov | reads coverage with A |
| Genecov | mean coverage for the whole gene |
| Ratio | A rate for the sites/ or methylation level for the sites |
| Pvalue | test for A rate based on the background |
| P_adjust | FDR ajusted P value |



