#!/bin/bash

###########################################################
# Tue Mar  5 14:04:36 EST 2019                            #
# This is a pipeline script for handling paired end       #
# PRO-seq data with UMIs on both ends of the read.        #
# Run this script in a directory that has one folder      #
# named "fastq" which contains the data.                  #
# Fastq files must have identical names other than        #
# ending in _R1.fastq and _R2.fastq.                      #
###########################################################

## Parameters
THREADS=50 # Threads to use for multithreaded applications
UMI_LEN=6  # Length of UMI in basepairs

## UMI Flags (set to Y or N as appropriate)
FIVEP_UMI="Y" # Is there a UMI on the 5' end of the read?
THREEP_UMI="Y" # Is there a UMI on the 3' end of the read?

## Adaptor sequences to clip. Default = Tru-Seq small RNA
ADAPTOR_1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" 
ADAPTOR_2="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"

## Genomes. Fill in paths.
GENOME_EXP="/home/jaj256/genome/dm6/dm6Hsp70AaOnly"
GENOME_SPIKE="/home/jaj256/genome/dm6hg38/dm6hg38" ## USE REPEAT MASKED VERSION!!
SPIKE_PREFIX="hg38" ## This is the prefix you've used on your spike in chromosomes, ie >hg38chr1
RDNA="/home/jaj256/genome/dm3hg38/dm3hg38rDNA"

## Mapq value for filtering multimappers
MAPQ=10

###############################################################
#                           PIPELINE                          #
###############################################################


# Unzipping if needed
echo "unzipping..."
for FILE in fastq/*
    do
        if [[ "$FILE" == *.gz ]]
        then 
            gunzip $FILE & 
        fi
    done
wait

# Removing extra info from filenames.
# This is general and works with files from Cornell BRC.
# If filenames are formatted differently, does nothing.
echo "renaming if needed..."
for FILE in $(ls fastq/)
    do 
        NEW=fastq/"$(echo "$FILE" | 
            sed 's/^[0-9]\+_[0-9]\+_[0-9]\+_[0-9A-Z]\+_//' | 
            sed 's/_[ATCG]\{6,8\}_/_/')"
        if [ ! -s "$NEW" ]
        then
            mv fastq/"$FILE" "$NEW"
        fi
    done

mkdir -p logs
mkdir -p logs/fastqc

### Running fastqc on files
echo "running fastqc if needed..."
for FILE in fastq/*.fastq
    do 
        if [ ! -s logs/fastqc/"$(basename ${FILE/.fastq/_fastqc.zip})" ]
        then
            fastqc "$FILE" -o logs/fastqc --quiet  &
        fi
    done
wait


mkdir -p trimmedFastq

### Autodetecting paired end files
echo "detecting paired end files..."
NUM=$(ls fastq | wc -l) 
NUM_REDUCED=$(ls fastq | sed 's/_R.*//' | uniq | wc -l)
if [[ $NUM == $NUM_REDUCED ]]
    then
        PAIRED="N"
        echo "detected ""$NUM"" single end fastq files. exiting..."
        exit
    else
        PAIRED="Y"
        echo "detected ""$NUM_REDUCED"" paired end fastq files"
fi

### Trimming adapters and filtering rRNA reads
echo "trimming adapters and filtering rDNA reads..."
mkdir -p logs/fastp
mkdir -p logs/rRNA
mkdir -p trimmedFastq
if [[ $PAIRED == "Y" ]]
    then
    ## Branches for either 3' UMI or both UMIs
    if [[ $THREEP_UMI == "Y" ]]
        then
        # Branch for both UMIs
        if [[ $FIVEP_UMI == "Y" ]]
            then
            for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
                do
                    if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]
                        then
                        echo "trimming adapters and filtering rRNA reads for "${PAIR}
                        (fastp \
                        -i fastq/${PAIR}_R1.fastq \
                        -I fastq/${PAIR}_R2.fastq \
                        --adapter_sequence $ADAPTOR_1 \
                        --adapter_sequence_r2 $ADAPTOR_2 \
                        --umi \
                        --stdout \
                        --umi_loc=per_read \
                        --umi_len=${UMI_LEN} \
                        --html logs/fastp/${PAIR}_fastp.html \
                        -w $(echo ${THREADS}/3 | bc)\
                        -c \
                        --overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
                        (bowtie2 \
                        --fast-local \
                        --un-conc trimmedFastq/${PAIR}.fastq \
                        --interleaved - \
                        -x ${RDNA} \
                        --threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
                    fi
                done
        # Branch for just 3' UMI
            else
            for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
                do
                    if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]
                    then
                    echo "trimming adapters and filtering rRNA reads for "${PAIR}
                    (fastp \
                    -i fastq/${PAIR}_R1.fastq \
                    -I fastq/${PAIR}_R2.fastq \
                    --adapter_sequence $ADAPTOR_1 \
                    --adapter_sequence_r2 $ADAPTOR_2 \
                    --umi \
                    --stdout \
                    --umi_loc=read1 \
                    --umi_len=${UMI_LEN} \
                    --html logs/fastp/${PAIR}_fastp.html \
                    -w $(echo ${THREADS}/3 | bc) \
                    -c \
                    --overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
                    (bowtie2 \
                    --fast-local \
                    --un-conc trimmedFastq/${PAIR}.fastq \
                    --interleaved - \
                    -x ${RDNA} \
                    --threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
                    fi
                done
            fi
        # Branch for only 5' UMI or no UMIs
        else
        # Branch for only 5' UMI
        if [[ $FIVEP_UMI == "Y" ]]
            then
            for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
                do
                    if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]
                        then
                        echo "trimming adapters and filtering rRNA reads for "${PAIR}
                        (fastp \
                        -i fastq/${PAIR}_R1.fastq \
                        -I fastq/${PAIR}_R2.fastq \
                        --adapter_sequence $ADAPTOR_1 \
                        --adapter_sequence_r2 $ADAPTOR_2 \
                        --umi \
                        --stdout \
                        --umi_loc=read2 \
                        --umi_len=${UMI_LEN} \
                        --html logs/fastp/${PAIR}_fastp.html \
                        -w $(echo ${THREADS}/3 | bc) \
                        -c \
                        --overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
                        (bowtie2 \
                        --fast-local \
                        --un-conc trimmedFastq/${PAIR}.fastq \
                        --interleaved - \
                        -x ${RDNA} \
                        --threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
                    fi
                done
                # Branch for no UMI
                else
                    for PAIR in $(ls fastq | sed 's/_R[1-2].*//' | uniq )
                    do
                        if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]
                            then
                            echo "trimming adapters and filtering rRNA reads for "${PAIR}
                            (fastp \
                            -i fastq/${PAIR}_R1.fastq \
                            -I fastq/${PAIR}_R2.fastq \
                            --adapter_sequence $ADAPTOR_1 \
                            --adapter_sequence_r2 $ADAPTOR_2 \
                            --stdout \
                            --html logs/fastp/${PAIR}_fastp.html \
                            -w $(echo ${THREADS}/3 | bc) \
                            -c \
                            --overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
                            (bowtie2 \
                            --fast-local \
                            --un-conc trimmedFastq/${PAIR}.fastq \
                            --interleaved - \
                            -x ${RDNA} \
                            --threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null
                        fi
                    done
                fi
            fi
        fi

### Cleaning up filenames in trimmedFastq (bowtie automatically names PE --un output)
for FILE in trimmedFastq/*1.fastq
do
    if [ ! -s ${FILE/.1.fastq/_R1.fastq} ]
    then
        mv "$FILE" ${FILE/.1.fastq/_R1.fastq}
    fi
done

for FILE in trimmedFastq/*2.fastq
do 
    if [ ! -s ${FILE/.2.fastq/_R2.fastq} ]
    then
        mv "$FILE" ${FILE/.2.fastq/_R2.fastq}
    fi
done


### Aligning to spike in genome to get normalization factors
mkdir -p spikeBAM
mkdir -p logs/spikeAlign


if [[ "$PAIRED" == "Y" ]]
then
    for PAIR in $(ls trimmedFastq | sed 's/_R[1-2].*//' | uniq )
    do
        if [ ! -s "spikeBAM/${PAIR}_hg38.BAM" ]
        then
            echo "aligning ${PAIR} to spike in genome"
            (bowtie2 \
            --local \
            --very-sensitive-local \
            --threads $(echo ${THREADS}/3*2 | bc) \
            --no-unal \
            --no-mixed \
            --no-discordant \
            -x "$GENOME_SPIKE" \
            -1 "trimmedFastq/${PAIR}_R1.fastq" \
            -2 "trimmedFastq/${PAIR}_R2.fastq" \
            2> logs/spikeAlign/${PAIR}_spikeAlign.log) |
            samtools view -hS -f 2 -q ${MAPQ} |
            perl -n -e 'print $_ if (/^\@/ || /'${SPIKE_PREFIX}'/ ) ' |
            samtools view -b | 
            samtools sort -@ $(echo ${THREADS}/3 | bc) -o spikeBAM/${PAIR}.BAM
            samtools index spikeBAM/${PAIR}.BAM
        fi
    done
fi



### Aligning to experimental genome

mkdir -p BAM
mkdir -p logs/align



if [[ "$PAIRED" == "Y" ]]
then
    for PAIR in $(ls trimmedFastq | sed 's/_R[1-2].*//' | uniq )
    do
        if [ ! -s "BAM/${PAIR}.BAM" ]
        then
            echo "aligning ${PAIR} to experimental genome"
            (bowtie2 \
            --local \
            --sensitive-local \
            --threads $(echo ${THREADS}/3*2 | bc) \
            -x "$GENOME_EXP" \
            -1 "trimmedFastq/${PAIR}_R1.fastq" \
            -2 "trimmedFastq/${PAIR}_R2.fastq" \
            2> logs/align/${PAIR}_align.log) |
            samtools view -bS -f 2 -q ${MAPQ} |
            samtools sort -@ $(echo ${THREADS}/3 | bc) -o BAM/${PAIR}.BAM 
            samtools index BAM/${PAIR}.BAM
        fi
    done
fi


### deduplicating with UMIs
mkdir -p BAMdeDuped
mkdir -p logs/deDup

for FILE in BAM/*.BAM
do
    if [ ! -s "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" ]
    then
        (umi_tools dedup \
        -I "$FILE" \
	    --umi-separator=":" \
        --paired \
        -S "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" \
        )> "logs/deDup/$(basename ${FILE%.BAM}_deDup.log)" &&
        samtools index "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" 
    fi
done


### deduplicating with UMIs
mkdir -p spikeBAMdeDuped
mkdir -p logs/spikedeDup

for FILE in spikeBAM/*.BAM
do
    if [ ! -s "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" ]
    then
        (
        umi_tools dedup \
        -I "$FILE" \
        --paired \
        --umi-separator=":" \
        -S "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" \
        )> "logs/spikedeDup/$(basename ${FILE%.BAM}_deDup.log)" &&
        samtools index "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" 
    fi
done

### Generating infoTable
mkdir -p info

if [ ! -s info/infoTable.tsv ]
then
touch info/infoTable.tsv
echo -e Name'\t'RawReads'\t'NonDimerReads'\t'%dimer'\t'insertSize'\t'rRNAreads'\t'%rRNA'\t'passedFilters'\t'\
    bowtieConcordant'\t'bowtieMulti'\t'bowtieUnal'\t'bowtieOverallMap%'\t'bowtieConcordant%'\t'\
    bowtieMulti%'\t'bowtieUnal%'\t'uniqueMapped'\t'uniqueMappedNondup'\t'%PCRdups'\t'uniqueMappedSpikein'\t'\
    uniqueMappedSpikeinNondup'\t'spikeInPCRdups% >> info/infoTable.tsv

for SAMPLE in $(ls BAM/*.BAM | sed 's/.BAM//' | sed 's/BAM\///' )
do
    NAME=${SAMPLE}
    RAW_READS=$(cat logs/fastp/${SAMPLE}_fastp.log |
                grep "total reads:" | head -n 1 | 
                awk '{print $3}')
    TRIMMED_READS=$(cat logs/fastp/${SAMPLE}_fastp.log |
                grep "total reads:" | tail -n 1 |
                awk '{print $3}')
    PER_DIMER=$(echo "(1-"${TRIMMED_READS}"/"${RAW_READS}")*100" | bc -l)%
    INSERT_SIZE=$(cat logs/fastp/${SAMPLE}_fastp.log |
                grep "Insert size peak" |
                awk '{print $8}')
    PASSED_FILTERS=$(cat logs/align/${SAMPLE}_align.log |
                grep "reads; of these:$" |
                awk '{print $1}')
    RRNA=$(echo ${TRIMMED_READS}"-"${PASSED_FILTERS} | bc )
    PER_RRNA=$(echo ${RRNA}"/"${RAW_READS}"*100" | bc -l)%
    B_CONC=$(cat logs/align/${SAMPLE}_align.log |
            grep "aligned concordantly exactly 1 time$" |
            awk '{print $1}')
    B_MULTI=$(cat logs/align/${SAMPLE}_align.log |
            grep "aligned concordantly >1 times$" |
            awk '{print $1}')
    B_UNAL=$(cat logs/align/${SAMPLE}_align.log |
            grep "aligned concordantly 0 times$" |
            awk '{print $1}')
    B_OAP=$(cat logs/align/${SAMPLE}_align.log |
            grep "overall alignment rate$" |
            awk '{print $1}')
    B_CONC_PER=$(echo ${B_CONC}"/"${PASSED_FILTERS}"*100" | bc -l)%
    B_MULTI_PER=$(echo ${B_MULTI}"/"${PASSED_FILTERS}"*100" | bc -l)%
    B_UNAL_PER=$(echo ${B_UNAL}"/"${PASSED_FILTERS}"*100" | bc -l)%
    UNIQ_MAPPED=$(cat logs/deDup/${SAMPLE}_deDup.log |
            grep "Input Reads:" | awk '{print $10}')
    UNIQ_MAPPED_DEDUP=$(cat logs/deDup/${SAMPLE}_deDup.log |
            grep "Number of reads out:" | awk '{print $8}')
    PER_DUPS=$(echo "(1-"${UNIQ_MAPPED_DEDUP}"/"${UNIQ_MAPPED}")*100" | bc -l)%
    UNIQ_MAPPED_SPIKE=$(cat logs/spikedeDup/${SAMPLE}_deDup.log |
            grep "Input Reads:" | awk '{print $10}')
    UNIQ_MAPPED_DEDUP_SPIKE=$(cat logs/spikedeDup/${SAMPLE}_deDup.log |
            grep "Number of reads out:" | awk '{print $8}')
    PER_DUPS_SPIKE=$(echo "(1-"${UNIQ_MAPPED_DEDUP_SPIKE}"/"${UNIQ_MAPPED_SPIKE}")*100" | bc -l)%
 
    echo -e $NAME'\t'\
    $RAW_READS'\t'\
    $TRIMMED_READS'\t'\
    $PER_DIMER'\t'\
    $INSERT_SIZE'\t'\
    $RRNA'\t'\
    $PER_RRNA'\t'\
    $PASSED_FILTERS'\t'\
    $B_CONC'\t'\
    $B_MULTI'\t'\
    $B_UNAL'\t'\
    $B_OAP'\t'\
    $B_CONC_PER'\t'\
    $B_MULTI_PER'\t'\
    $B_UNAL_PER'\t'\
    $UNIQ_MAPPED'\t'\
    $UNIQ_MAPPED_DEDUP'\t'\
    $PER_DUPS'\t'\
    $UNIQ_MAPPED_SPIKE'\t'\
    $UNIQ_MAPPED_DEDUP_SPIKE'\t'\
    $PER_DUPS_SPIKE  >> info/infoTable.tsv
            
done
fi

# Making non-normalized bigWig files
mkdir -p bw
for FILE in BAMdeDuped/*.BAM
do
    if [ ! -s "bw/$(basename ${FILE/.BAM/_fwd.bw})" ]
    then
        bamCoverage \
        --bam $FILE \
        --skipNonCoveredRegions \
        --outFileName bw/$(basename ${FILE/.BAM/_fwd.bw}) \
        --binSize 1 \
        --numberOfProcessors ${THREADS} \
        --normalizeUsing None \
        --Offset 1 \
        --samFlagInclude 82
    fi
    if [ ! -s "bw/$(basename ${FILE/.BAM/_rev.bw})" ]
    then
        bamCoverage \
        --bam $FILE \
        --skipNonCoveredRegions \
        --outFileName bw/$(basename ${FILE/.BAM/_rev.bw}) \
        --binSize 1 \
        --numberOfProcessors ${THREADS} \
        --normalizeUsing None \
        --Offset 1 \
        --samFlagInclude 98
    fi
done

