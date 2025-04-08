#!/usr/bin/perl
use IO::All;
use Modern::Perl;
use Data::Dumper;

my $infile = shift;

   my $basename = io($infile)->filename;
   #say "seqkit rmdup -s -j 20 $infile -o $basename";
   #`seqkit rmdup -s -j 20 $infile -o $basename `;
   say "trim_galore -j 7 -q 20 --length 20  -a AACTGTAGGCACCATCAAT  $basename";
   `/lustre1/chengqiyi_pkuhpc/softwares/anaconda3/bin/trim_galore -j 7 -q 20 --length 20 -a AACTGTAGGCACCATCAAT  $basename`;
   my $index='/lustre2/chengqiyi_pkuhpc/yanxueqing/public_data/human/rRNA/rRNA_uniq';
   my $index2='/lustre2/chengqiyi_pkuhpc/yanxueqing/public_data/human/tRNA-GtRNAdb/hg38-tRNAs_uniq.fa';
   my $index3='/lustre2/chengqiyi_pkuhpc/yanxueqing/public_data/human/mRNA_gencode/gencode.v44.longest.transcripts.fa';
   my $cleanfq = $basename;
   $cleanfq =~s/.fq.gz/_trimmed.fq.gz/;

   my $fq_no_rRNA = $cleanfq;
   my $fq_no_rRNA_tRNA = $cleanfq;
   $fq_no_rRNA =~s/.fq.gz/_no_rRNA.fq/g;
   $fq_no_rRNA_tRNA =~s/.fq.gz/_no_rRNA_tRNA.fq/g;

   say "bowtie -p 20 -S --no-unal --un $fq_no_rRNA  -x $index $cleanfq $cleanfq.sam  2>$cleanfq.log && rm $cleanfq.sam";
   `bowtie -p 20 -S --no-unal --un $fq_no_rRNA  -x $index $cleanfq $cleanfq.sam 2>$cleanfq.log && rm $cleanfq.sam`;

   say "bowtie -p 20 -S --no-unal --un $fq_no_rRNA_tRNA  -x $index2 $fq_no_rRNA $fq_no_rRNA.sam 2>$fq_no_rRNA.log && rm $fq_no_rRNA.sam";
   `bowtie -p 20 -S --no-unal --un $fq_no_rRNA_tRNA  -x $index2 $fq_no_rRNA $fq_no_rRNA.sam 2>$fq_no_rRNA.log && rm $fq_no_rRNA.sam`;

   say "bowtie -p 20 -S --no-unal -x $index3 $fq_no_rRNA_tRNA $fq_no_rRNA_tRNA.sam 2>$fq_no_rRNA_tRNA.log";
   `bowtie -p 20 -S --no-unal -m 1 -v 1 -a --best --strata  -x $index3 $fq_no_rRNA_tRNA $fq_no_rRNA_tRNA.sam 2>$fq_no_rRNA_tRNA.log`;

   #bowtie -p 20 -S --no-unal -x /lustre2/chengqiyi_pkuhpc/yanxueqing/public_data/human/mRNA_gencode/gencode.v44.longest.transcripts.fa -m 3 -a --best --strata test_trimmed_no_rRNA_tRNA.fq test_trimmed_no_rRNA_tRNA.fq2.sam 2>test_trimmed_no_rRNA_tRNA.fq2.log
   say "samtools sort -@ 20 -o $fq_no_rRNA_tRNA.sort.bam $fq_no_rRNA_tRNA.sam";
   `samtools sort -@ 20 -o $fq_no_rRNA_tRNA.sort.bam $fq_no_rRNA_tRNA.sam`;
   `rm $fq_no_rRNA_tRNA.sam`;

