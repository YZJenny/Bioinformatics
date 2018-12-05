# 1. Installation
mkdir homer && cd homer
wget http://homer.salk.edu/homer/configureHomer.pl
perl configureHomer.pl -install

# 2. Download related data
perl configureHomer.pl -install hg19      #GENOMES: human genome and annotation for UCSC hg19
perl configureHomer.pl -install hg38      #GENOMES: human genome and annotation for UCSC hg38
perl configureHomer.pl -install human-p   #PROMOTERS: promoter set(-2000,2000)bp
perl configureHomer.pl -install human-o   #ORGANISMS: Homo sapiens (human) accession and ontology information, NCBI Gene

# 3. finding motif
# 3.1 FASTA file Motif Discovery
findMotifs.pl - performs motif analysis with lists of Gene Identifiers or FASTA files (See FASTA file analysis)
Usage:  
      findMotifs.pl targets.fa fasta output/ -fasta background.fa
Example:
      findMotifs.pl pipeline/Homer/input/H3K27me3_peaks.fa fasta pipeline/Homer/output/fa_output/ -len 8,10,12

# 3.1 Gene/Promoter-based Analysis
findMotifs.pl - performs motif and gene ontology analysis with lists of Gene Identifiers, both promoter and mRNA motifs (See Gene ID Analysis Tutorial)
Usage:
      findMotifs.pl <input list> <promoter set> <output directory> [additoinal options]
Example:
      findMotifs.pl pipeline/Homer/input/genelist.txt human pipeline/Homer/output/genelst_output/ -len 8,10,12

#3.3 Next-Gen Sequencing/Genomic Position Analysis
findMotifsGenome.pl: Program will find de novo and known motifs in regions in the genome
Usage: 
      findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]

Example:
      awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' input/H3K27me3_peaks.broadPeak > input/H3K27me3_peaks_homer.bed
      findMotifsGenome.pl pipeline/Homer/input/H3K27me3_peaks_homer.bed hg19 pipeline/Homer/output/findMotGen -len 8,10,12

# 4. Another analysis
# 4.1 Functional Enrichment Analysis
findGO.pl pipeline/Homer/input/genelist.txt human pipeline/Homer/output/GO_output/
