mm10_index=/home/genomewide/refgenome/mm10/mm10
mm10_seq=/home/genomewide/refgenome/mm10/mm10.fa
mm10_gtf=/home/genomewide/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf

mkdir GTF/ DIFF_N_T/ Figure/
realpath */cufflinks/transcripts.gtf > GTF/assemblies.txt

wp=/home/yzj/ZF/0319/
cuffmerge -g $mm10_gtf -s $mm10_seq -p 18 -o ${wp}GTF ${wp}GTF/assemblies.txt

N1=${wp}N1/N1.bam
N2=${wp}N2/N2.bam
N3=${wp}N3/N3.bam

T1=${wp}YF-T1/YF-T1.bam
T2=${wp}YF-T2/YF-T2.bam
T3=${wp}YF-T3/YF-T3.bam
T4=${wp}YF-T4/YF-T4.bam

merged_gtf=${wp}GTF/merged.gtf
cuffdiff -o ${wp}DIFF_N_T -p 18 -L N,T -u $merged_gtf $N1,$N2,$N3 $T1,$T2,$T3,$T4 
