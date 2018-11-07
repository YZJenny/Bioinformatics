# 1.Installation
#1.1 Download the software from http://meme-suite.org/doc/download.html/
#1.2 Type the following commands
tar zxf meme_5.0.2.tar.gz
cd meme_5.0.2
./configure --prefix=$HOME/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make test
make install 
export PATH=$HOME/meme/bin:$PATH 
# 1.3 Prerequisite software
#(1) Perl
    cd scripts
    perl dependencies.pl
#(2) Python
#(3) zlib
#(4) Ghostscript - for creating PNG files.
#(5) Assorted common utilities

# 2. Download Databasea
wget http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz
wget http://meme-suite.org/meme-software/Databases/gomo/gomo_databases.3.2.tgz

# 3. get fasta
bedtools getfasta -fi /home/genomewide/refgenome/hg19/hg19.fa -bed pipeline/MEME/data/H3K27me3_peaks.broadPeak -fo pipeline/MEME/data/H3K27me3_peaks.fa &

# 4. run MEME
# 4.1 meme
meme pipeline/MEME/input/H3K27me3_peaks.fa -o pipeline/MEME/output/meme -dna &

# 4.2 dreme
dreme -p pipeline/MEME/input/H3K27me3_peaks.fa -o pipeline/MEME/output/dreme -dna &

# 4.3 meme-chip
dreme -p pipeline/MEME/input/H3K27me3_peaks.fa -o pipeline/MEME/output_dreme -dna -db meme/db/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme &


