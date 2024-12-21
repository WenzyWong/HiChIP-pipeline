######################################################
# HiChIP - H3K9me3 in Panc-1 cell line with two groups
# Yunzhe Wang
# Updated: 2024-10-24
#
# Step 01. Traditional ChIP-seq analysis
#
######################################################
# Originally running on server-17
DATA="/data/yzwang/co_work/zxliu/hichip2/Rawdata/"
OUT="/data/yzwang/co_work/zxliu/hichip2/chip_out"
REF="/data/zxwang/ref/hg19/hg19.fa"
THREADS=64

mkdir -p ${OUT}
mkdir -p ${OUT}"/01_cleandata"
ls ${DATA} | while read SAMPLE
do
    TMPOUT=${OUT}"/01_cleandata/"${SAMPLE}
    mkdir -p ${TMPOUT}
    echo -e "\033[31m Trimming for ${SAMPLE} \033[0m"
    trim_galore -q 20 -o ${TMPOUT} --length 50 --paired --gzip ${DATA}/${SAMPLE}/${SAMPLE}"_R1.fq.gz" ${DATA}/${SAMPLE}/${SAMPLE}"_R2.fq.gz"
done

mkdir -p ${OUT}"/02_alignment"
ls ${DATA} | while read SAMPLE
do
    TMPOUT=${OUT}"/02_alignment/"${SAMPLE}
    mkdir -p ${TMPOUT}
    echo -e "\033[31m Aligning for ${SAMPLE} \033[0m" 
    bowtie2 -p ${THREADS} -x ${REF} -1 ${OUT}/01_cleandata/${SAMPLE}/${SAMPLE}"_R1_val_1.fq.gz" -2 ${OUT}/01_cleandata/${SAMPLE}/${SAMPLE}"_R2_val_2.fq.gz" | samtools sort -@ ${THREADS} -o ${TMPOUT}/${SAMPLE}.sorted.bam

    echo -e "\033[31m Filtering for ${SAMPLE} \033[0m"
    samtools index -@ ${THREADS} ${TMPOUT}/${SAMPLE}.sorted.bam 
    samtools view -@ ${THREADS} -q 30 -F 3844 -bo ${TMPOUT}/${SAMPLE}.af.bam ${TMPOUT}/${SAMPLE}.sorted.bam
    samtools index -@ ${THREADS} ${TMPOUT}/${SAMPLE}.af.bam 
    samtools flagstat ${TMPOUT}/${SAMPLE}.af.bam > ${TMPOUT}/${SAMPLE}.af.metrics.txt
    mkdir -p ${TMPOUT}/tmp
    java -jar /data/zxwang/public_tools/ChIP-seq/tools/picard.jar MarkDuplicates I=${TMPOUT}/${SAMPLE}.af.bam O=${TMPOUT}/${SAMPLE}.df.bam METRICS_FILE=${TMPOUT}/${SAMPLE}.df.metrics.txt REMOVE_DUPLICATES=true TMP_DIR=${TMPOUT}/tmp
    samtools index -@ ${THREADS} ${TMPOUT}/${SAMPLE}.df.bam
    samtools flagstat ${TMPOUT}/${SAMPLE}.df.bam > ${TMPOUT}/${SAMPLE}.df.metrics.txt
done

mkdir -p ${OUT}"/03_bigwig"
ls ${DATA} | while read SAMPLE
do
    TMPOUT=${OUT}"/03_bigwig/"${SAMPLE}
    mkdir -p ${TMPOUT}
    echo -e "\033[31m Calculating coverage for ${SAMPLE} \033[0m"
    bamCoverage --bam ${OUT}/02_alignment/${SAMPLE}/${SAMPLE}.df.bam -o ${TMPOUT}/${SAMPLE}.bw -bs 10 --normalizeUsing CPM --effectiveGenomeSize 2862010578 --extendReads 200
done

mkdir -p ${OUT}"/04_macs3_peak"
ls ${DATA} | while read SAMPLE
do
    TMPOUT=${OUT}"/04_macs3_peak/"${SAMPLE}
    mkdir -p ${TMPOUT}
    echo -e "\033[31m Calling peaks for ${SAMPLE} \033[0m"
    macs3 callpeak -t ${OUT}/02_alignment/${SAMPLE}/${SAMPLE}.df.bam -g hg19 -n 1 -B -p 0.001 -f BAMPE --outdir ${TMPOUT}
done