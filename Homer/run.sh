
1.安装
mkdir homer && cd homer
wget http://homer.salk.edu/homer/configureHomer.pl
perl configureHomer.pl -install
perl configureHomer.pl -install hg19
perl configureHomer.pl -install hg38

2.DNA motif的查找使用方法
(1) Next-Gen Sequencing/Genomic Position Analysis:
findMotifsGenome.pl: Program will find de novo and known motifs in regions in the genome

Usage: findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]

Example:
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' input/H3K27me3_peaks.broadPeak > input/H3K27me3_peaks_homer.bed
findMotifsGenome.pl input/H3K27me3_peaks_homer.bed hg19 output/motifDir -len 8,10,12

(2) Gene/Promoter-based Analysis:
findMotifs.pl: Program will find de novo and known motifs in a gene list  

Usage:  findMotifs.pl <input list> <promoter set> <output directory> [additoinal options]

Example: findMotifs.pl input/genelist.txt human output/motifResults/ -len 8,10,12
