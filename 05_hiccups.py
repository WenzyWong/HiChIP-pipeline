######################################################
# HiChIP - H3K9me3 in Panc-1 cell line with two groups
# Yunzhe Wang
# Updated: 2024-11-25
#
# Step 05: Use HICCUPS to calculate loops
#
######################################################
# Originally running on server-220

REP_INDEX = {"1", "2"}
GROUP = {"KO-11-", "wt"}

rule all:
    input:
        expand("241126-hiccups_out/HiChIP-Panc-1-{gr}{rep}/", gr=GROUP, rep=REP_INDEX)

rule hiccups:
    input:
        "241030-juicer_file/HiChIP-Panc-1-{gr}{rep}.allValidPairs.hic"
    output:
        directory("241126-hiccups_out/HiChIP-Panc-1-{gr}{rep}/")
    log:
        "241126-hiccups_out/HiChIP-Panc-1-{gr}{rep}.log"
    params:
        "/data/yzwang/tool/juicer/juicer_tools_1.22.01.jar",
        "hiccups"
    shell:
        "java -jar {params[0]} {params[1]} --cpu --threads 16 -r 100000 --ignore-sparsity {input[0]} {output[0]} > {log} 2>&1 &"

# Run: snakemake -s 05_hiccups.py -p -j 4