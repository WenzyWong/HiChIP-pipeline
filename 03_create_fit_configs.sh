######################################################
# HiChIP - H3K9me3 in Panc-1 cell line with two groups
# Yunzhe Wang
# Updated: 2024-11-03
#
# Step 03. Create a config file for FitHiChIP
#
######################################################
# Originally running on server-220
TOOL="/data/yzwang/tool/FitHiChIP/FitHiChIP_HiCPro.sh"
INPUT="/home/yzwang/data/cowork/zxliu/hichip/out/hic_results/data/"
CHR="/data/yzwang/cowork/zxliu/hichip/chrom_hg19.sizes"
OUT="/data/yzwang/cowork/zxliu/hichip/out/out_fithichip/output/"
ls ${INPUT} | while read sample
do
    VALIDPAIRS=${INPUT}/${sample}/${sample}".allValidPairs"
    CONFIG="/data/yzwang/cowork/zxliu/hichip/out/out_fithichip/configs/configfile_P2P_CoverageBias_"${sample}
    mkdir ${OUT}/${sample}
    PEAK="/data/yzwang/cowork/zxliu/hichip/chip/"${sample}"/1_summits.bed"
    echo "
#==================================== 
# Sample configuration file for running FitHiChIP
#====================================  

##********
##***** parameters to provide the input HiChIP alignment files
##********

##============
## option 1: provide the valid pairs from HiC-Pro pipeline - can be gzipped as well
##============
ValidPairs=${VALIDPAIRS}

##********
## File containing chromomosome size information corresponding to the reference genome.
##********
ChrSizeFile=${CHR}


##********
## Mandatory parameter - Reference ChIP-seq / HiChIP peaks (in .bed format) - can be gzipped as well
## We recommend using reference ChIP-seq peaks (if available)
## Otherwise, peaks can be computed from HiChIP data. 
## See the documentation: https://ay-lab.github.io/FitHiChIP/usage/Utilities.html#inferring-peaks-from-hichip-data-for-use-in-the-hichip-pipeline
##********
PeakFile=${PEAK}


##********
## Mandatory parameter - Output directory to contain all the results
##********
OutDir=${OUT}/${sample}


##********
## Mandatory parameter - Boolean variable indicating if the reference genome is circular
## 0, by default. If 1 (circular genome), calculation of genomic distance is slightly different
##********
CircularGenome=0


##********
##***** Various FitHiChIP loop calling related parameters
##********

##Interaction type
## 1: peak to peak 
## 2: peak to non peak 
## 3: peak to all (default - both peak-to-peak and peak-to-nonpeak) 
## 4: all to all (similar to Hi-C)
## 5: All of the modes 1 to 4 are computed.
IntType=3

## Bin size, in bases, for the interactions. Default = 5000 (5 Kb).
BINSIZE=5000

## Lower distance threshold of loops - default = 20000 (20 Kb)
LowDistThr=20000

## Upper distance threshold of loops - default = 2000000 (2 Mb)
UppDistThr=2000000

## Values 0/1 - Applicable if IntType = 3 (peak to all output interactions)
## 1 indicates FitHiChIP(S) model - uses only peak to peak loops for background modeling
## 0 corresponds to FitHiChIP(L) - uses both peak to peak and peak to nonpeak loops for background modeling
UseP2PBackgrnd=0

## type of bias - values: 1 / 2
## 1: coverage bias regression
## 2: ICE bias regression
BiasType=1

## if 1 (default), merge filtering (corresponding to either FitHiChIP(L+M) or FitHiChIP(S+M) 
## depending on the parameter UseP2PBackgrnd) is enabled
MergeInt=1

## FDR (q-value) threshold for loop significance
QVALUE=0.05

## prefix string of all the output files (Default = 'FitHiChIP').
PREFIX=FitHiChIP

## Binary variable 1/0: 
## if 1, overwrites any existing output file. 
## otherwise (0), does not overwrite any output file.
OverWrite=0
" > ${CONFIG}
done

# After excecuting this script, 
# there will be a config file for all samples apending analysis.