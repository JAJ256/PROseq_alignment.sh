# PROseq_alignment.sh

This is a pipeline script for aligning paired-end PRO-seq data that has cells of a different species spiked in for normalization, and uses some combination of random UMI sequences on the ligation end of either the 5' or 3' adapter, or both.

Run this script in a directory that has one folder named "fastq" which contains the data. Fastq files must have identical names other than ending in _R1.fastq and _R2.fastq.

# Parameters

These parameters are at the top of the file and can be changed depending on your data.

THREADS: (Integer) Number of threads to spawn for each step in the process.

UMI_LEN: (Integer) length of the UMI in basepairs. The if both 5' and 3' UMIs are used, this is the lenght for both. Fastp cannot handle multiple UMIs of different length.


FIVEP_UMI: (String) Set to "Y" if the UMI is on the 5' adapter

THREEP_UMI: (String) Set to "Y" if the UMI is on the 3' adapter

Note: both FIVEP_UMI and THREEP_UMI can be "Y" if there are UMIs on both sides of the insert

ADAPTOR_1 and ADAPTOR_2: (String) adapter sequences to trim. Default is TruSeq Small RNA sequences. These sequences are only here for backup in the case that fastp cannot automatically determine the adapter sequences by overlap analysis.


GENOME_EXP: (String) Path to the bowtie2 index for your experimental genome

GENOME_SPIKE: (String) Path to the bowtie2 index for your spike-in genome. To prepare a spike-in genome, combine your experimental genome with a repeat-masked version of your spike-in organism's genome. You must first modify the chromosome labels of the spike-in genome so alignments can be sorted later. 

SPIKE_PREFIX: (String) This is the prefix you've used on your spike in chromosomes, ie >spikechr1

RDNA: (String) Path to the bowtie2 index for the rDNA repeat for your organism(s)

MAPQ: (Integer) Mapq score cutoff for filtering multimappers




