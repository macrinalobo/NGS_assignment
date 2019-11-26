# NGS assignment

instructions to run

make sure gatk, samtools are installed and paths are correct

Run:

python q1.py

python q2.py

solutions to assignment:

Files for part 1:

q1.py - code for question 1

read_length_histogram.png - histogram for q1

cumulative_distribution_of_2mer_frequencies.png - distribution for q2

terminal_output_for_q1.txt - terminal output for q1

Part 2

Files for part 2:

q2.py - python code for part 2

INPUT: downloaded bam and bai file

HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam //excluded due to size

HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai

OUTPUT:

HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.sorted.bam - sorted bam //excluded due to size

insert_size_distribution.pdf - distribution of insert sizes

insert_size_information.txt - statistics + histogram bin values of insert sizes

terminal_output_for_q2.txt - terminal output for q2


STEPS:
1) ftp file download
2) Header view: samtools view -H <bamfile>, Fileview: samtools view <bamfile> | less - just to look inside the file
3) sort: samtools sort ( in this case it appears its already sorted )
4) get summary stats - samtools flagstat - for my understanding, to look at duplicate, unmapped reads
5) Picard CollectInsertMetrics to get and plot insert size distribution for paired end mapped reads - unpaired and duplicate reads are excluded

ASSUMPTIONS - documented in code as well
    
    - ASSUMPTION 1: fragment size approximately equals insert size
    Here I ignore adapter length ( usually trimmed ) since actually fragment size = insert size + length of adapters at both ends
    
    - ASSUMPTION 2: insert size approximately equals TLEN ( field 9 of sam ) which I believe Picard CollectInsertMetrics uses for its plots and stats
    This assumption would be violated in case of RNAseq data due to splicing and hence introns not being sequenced but that is not relevant here
    
    - ASSUMPTION 3: In some places it seems "insert size" is used interchangeably with "inner mate size" but not here.
    I assume inner mate size to be TLEN - len(R1) - len(R2) and do not use it anywhere.

NOTE:

CollectInsertMetrics has some filters which are implicitly set and I chose not to change:
1) duplicates are excluded, as they should be
2) extreme insert sizes ( 10 deviations away from median ) are excluded 
3) 1431742 read pairs out of 2931500 reads meet the criteria here

Insert size and somatic alterations:

Abnormal insert sizes can be a sign of somatic structural variation in those regions.
For a typical standardized target size, the insert size distribution should be unimodal with mode around the target size.
tails in the insert size distribution could be due to alignment errors / artifacts. However, multimodal insert size distributions
are indicative of somatic structural variation ( eg deletion, duplication, insertion, inversion, translocation)




