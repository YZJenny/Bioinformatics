


1. get fasta
bedtools getfasta -fi /home/genomewide/refgenome/hg19/hg19.fa -bed pipeline/MEME/data/H3K27me3_peaks.broadPeak -fo pipeline/MEME/data/H3K27me3_peaks.fa &

2. run meme
meme pipeline/MEME/data/H3K27me3_peaks.fa -o pipeline/MEME/output -dna &
