######################################################
# HiChIP - H3K9me3 in Panc-1 cell line with two groups
# Yunzhe Wang
# Updated: 2024-11-04
#
# Step 03: Run FitHiChIP
#
######################################################
# Originally running on server-220

# Copy "/data/yzwang/tool/Hic-Pro/annotation/chrom_hg19.sizes" to the current directory

TOOL="/data/yzwang/tool/FitHiChIP/FitHiChIP_HiCPro.sh"
INPUT="/data/yzwang/cowork/zxliu/hichip/out/out_fithichip/configs"
ls ${INPUT} | while read CONFIG
do
    bash ${TOOL} -C ${INPUT}/${CONFIG}
done

# This will generate the output files in "./out/out_fithichip/output/"