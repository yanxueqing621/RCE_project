
sample=$1
outprefix=$1
echo "trim_galore --paired -j 7 -q 20  --length 35 --basename $outprefix  ${sample}_1.fq.gz ${sample}_2.fq.gz"
trim_galore --paired -j 7 -q 20  --length 35 --basename $outprefix  ${sample}_1.fq.gz ${sample}_2.fq.gz
echo "bowtie2 -p 20 -x rRNA.fas --no-unal -1 ${outprefix}_val_1.fq.gz -2 ${outprefix}_val_2.fq.gz -S $outprefix.sam --un-conc ${outprefix}.fq 2>$outprefix.bt2.log && rm $outprefix.sam"
bowtie2 -p 20 -x rRNA.fas --no-unal -1 ${outprefix}_val_1.fq.gz -2 ${outprefix}_val_2.fq.gz -S $outprefix.sam --un-conc ${outprefix}.fq 2>$outprefix.bt2.log && rm $outprefix.sam
hisat2 -p 20 --no-unal  --dta -x hg38_trans_hisat  -1 $outprefix.1.fq -2 $outprefix.2.fq 2>$outprefix.his.log |samtools sort -@ 10 -o $outprefix.bam
featureCounts -T 20 -t exon -g gene_id -a GCF_000001405.39_GRCh38.p13_genomic.gtf -o $outprefix  $outprefix.bam
