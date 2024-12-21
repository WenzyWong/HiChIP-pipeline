######################################################
# HiChIP - H3K9me3 in Panc-1 cell line with two groups
# Yunzhe Wang
# Updated: 2024-10-26
#
# Step 02: Run HiC-Pro
#
######################################################
# Originally running on server-220

# Copy "/data/yzwang/tool/Hic-Pro/config-hicpro.txt" and "/data/yzwang/tool/Hic-Pro/annotation/HindIII_resfrag_hg19.be6" to current directory

HiC-Pro -i "Rawdata/" -o "out/" -c "config-hicpro.txt"

# After successfully run this script,
# it will generate "./out/hic_results/" to store all the results.
# "./out/hic_results/data/" containing ".ValidPairs"(the most important files);
# "./out/hic_results/matrix/" containing "raw" and "iced" matrix, and raw bed annotation files;
# "./out/hic_results/stat/" containing statistics to analysing QC.
