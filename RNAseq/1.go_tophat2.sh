mm10_index=/home/genomewide/refgenome/mm10/mm10
mm10_seq=/home/genomewide/refgenome/mm10/mm10.fa
mm10_gtf=/home/genomewide/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf
fastqc=/home/zhangfeng/tools/FastQC/fastqc
tophat2=/home/zhangfeng/disk/project/Cuffdiff/tools/tophat-2.1.0.Linux_x86_64/tophat2
cpu=30

wp=/home/yzj/ZF/0319/YF-T1/
line=L1
read1=$wp'V300012005_L1_HK500MUSwbrEAAJRAAPEI-591_1.fq'
read2=$wp'V300012005_L1_HK500MUSwbrEAAJRAAPEI-591_2.fq'
$fastqc $read1 $read2 -o $wp'/fastqc'
$tophat2 -p $cpu -o $wp$line\.tophatdir $mm10_index $read1 $read2

line=L2
read1=$wp'V300012005_L2_HK500MUSwbrEAAJRAAPEI-591_1.fq'
read2=$wp'V300012005_L2_HK500MUSwbrEAAJRAAPEI-591_2.fq'
$fastqc $read1 $read2 -o $wp'/fastqc'
$tophat2 -p $cpu -o $wp$line\.tophatdir $mm10_index $read1 $read2

samtools merge ${wp}YF-T1.bam ${wp}L1.tophatdir/accepted_hits.bam ${wp}L2.tophatdir/accepted_hits.bam
cufflinks -p $cpu -u -G $mm10_gtf -o ${wp}cufflinks ${wp}YF-T1.bam
