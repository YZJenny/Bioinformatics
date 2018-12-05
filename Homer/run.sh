
1.安装
mkdir homer && cd homer
wget http://homer.salk.edu/homer/configureHomer.pl
perl configureHomer.pl -install
perl configureHomer.pl -install hg19
perl configureHomer.pl -install hg38
perl configureHomer.pl -install human-p
perl configureHomer.pl -install human-o


2.DNA motif的查找使用方法
2.1 FASTA file Motif Discovery
(1) findMotifs.pl - performs motif analysis with lists of Gene Identifiers or FASTA files (See FASTA file analysis)
Usage:  
      findMotifs.pl targets.fa fasta motifResults/ -fasta background.fa
Example:
      findMotifs.pl input/H3K27me3_peaks.fa fasta output/fa_output/ -len 8,10,12

2.1 Gene/Promoter-based Analysis
(1) findMotifs.pl - performs motif and gene ontology analysis with lists of Gene Identifiers, both promoter and mRNA motifs (See Gene ID Analysis Tutorial)
Usage:
      findMotifs.pl <input list> <promoter set> <output directory> [additoinal options]
Example:
      findMotifs.pl input/genelist.txt human output/genelst_output/ -len 8,10,12


(2) findGO.pl - performs only gene ontology analysis with lists of Gene Identifiers
Usage:
      findGO.pl <target ids file> <organism> <output directory> [options]
Example:


2.3 Next-Gen Sequencing/Genomic Position Analysis
(1) findMotifsGenome.pl: Program will find de novo and known motifs in regions in the genome
Usage: 
      findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]

Example:
      awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' input/H3K27me3_peaks.broadPeak > input/H3K27me3_peaks_homer.bed
      findMotifsGenome.pl input/H3K27me3_peaks_homer.bed hg19 output/findMotGen -len 8,10,12
