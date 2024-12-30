######################################################
# HiChIP - H3K9me3 in Panc-1 cell line with two groups
# Yunzhe Wang
# Updated: 2024-12-29
#
# Step 05: Differential analysis on ChIP results
#
######################################################
# Originally running on server-220

setwd("/data/yzwang/cowork/zxliu/hichip/chip/")

library(DiffBind)
library(paletteer)
library(ggplot2)

sample <- read.csv("./sample_info.csv")
chip <- dba(sampleSheet = "./sample_info.csv")

chip <- dba.count(chip)
saveRDS(chip, file = "./raw_counts.rds")

mycol <- paletteer_d("ggsci::alternating_igv")
dba.plotPCA(chip, label = DBA_ID, vColors = mycol) # sample_pca_rmOutlier, 6.5 *6.5

# Reads count
info <- dba.show(chip)
libsizes <- cbind(LibReads = info$Reads, 
                  FRiP = info$FRiP, # Fraction of Reads in Peaks
                  PeakReads = round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes

# Normalization
chip <- dba.normalize(chip)
norm <- dba.normalize(chip, bRetrieve = T)

normlibs <- cbind(FullLibSize = norm$lib.sizes, 
                  NormFacs = norm$norm.factors,
                  NormLibSize = round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID

save(norm, file = "./norm.data")

# Plotting default heatmap
chip$config$RunParallel <- T
chip$config$cores <- 120
profile.default <- dba.plotProfile(chip, merge = DBA_ID, labels = DBA_ID)
dba.plotProfile(profile.default, labels = DBA_CONDITION) 

profile.condition <- dba.plotProfile(chip, merge = DBA_REPLICATE)
dba.plotProfile(profile.condition) # total_peaks_comparison, 8 * 3

# Differential analysis
ko <- dba.contrast(chip, minMembers = 2, 
                     design = "~Treatment",
                     contrast = c("Treatment", "ko", "WT"),
                     reorderMeta = list(Treatment = "WT"))

ko <- dba.analyze(ko)
save(ko, file = "./ko_contrast.data")

# Export result
ko.db <- dba.report(ko)
ko_res <- as.data.frame(ko.db)
write.csv(ko_res, "./ko_res.csv")

ko_bed <- ko_res[which(ko_res$FDR < 0.05), 
                     c("seqnames", "start", "end", "strand", "Fold")]
write.table(ko_bed, file = "./ko_res.bed", sep = "\t", 
            quote = F, row.names = F, col.names = F)

#####################################################################
# Change the source code in the built-in function of dba.plotVolcano:
# "DiffBind:::pv.DBAplotVolcano"
# Get the original function
modified_fun <- DiffBind:::pv.DBAplotVolcano

# Modify the function definition
body(modified_fun) <- substitute({
  if (missing(contrast)) {
    contrast <- 1:length(pv$contrasts)
  }
  else {
    if (contrast > length(pv$contrasts)) {
      stop("Specified contrast number is greater than number of contrasts", 
           call. = FALSE)
      return(NULL)
    }
  }
  for (con in 1:length(contrast)) {
    conrec <- pv$contrasts[[contrast[con]]]
    name1 <- conrec$name1
    name2 <- conrec$name2
    if (bFlip) {
      name1 <- conrec$name2
      name2 <- conrec$name1
    }
    for (meth in method) {
      res <- pv.DBAreport(pv, contrast = contrast[con], 
                          method = meth, bUsePval = TRUE, th = 100, bNormalized = TRUE, 
                          bFlip = bFlip, precision = 0, lfc = fold)
      if (!is.null(res)) {
        if (bUsePval) {
          vals <- res$"p-value"
          tstr <- "p" 
          res = mutate(res, Legend = case_when(
            vals > th ~ "NS",
            vals <= th & Fold >= fold ~ "Gain",
            vals <= th & Fold <= -fold ~ "Loss",
            TRUE ~ "NS"
          ))
        }
        else {
          vals <- res$FDR
          tstr <- "FDR" 
          res = mutate(res, Legend = case_when(
            vals > th ~ "NS",
            vals <= th & Fold >= fold ~ "Gain",
            vals <= th & Fold <= -fold ~ "Loss",
            TRUE ~ "NS"
          ))
        }
        
        idx <- vals <= th & abs(res$Fold) >= fold
        sigSites <- res[idx, ]
        if (sum(idx) > 0) {
          rownames(sigSites) <- 1:sum(idx)
        }
        res <- cbind(0, res)
        colnames(res)[1] <- "SiteNum"
        if (sum(idx) > 0) {
          res[idx, 1] <- 1:sum(idx)
          sidx <- sum(idx)
        }
        else {
          sidx <- 0
        }
        constr <- pv.getContrastString(conrec, bFlip)
        plotTitle <- sprintf("%s Contrast: %s  [%s %s<=%1.3f", 
                             facname, constr, sidx, tstr, th)
        if (fold > 0) {
          plotTitle <- sprintf("%s & abs(Fold)>=%1.2f]", 
                               plotTitle, 2^fold)
        }
        else {
          plotTitle <- sprintf("%s]", plotTitle)
        }
        xLabel <- "log2 Fold Change"
        yLabel <- sprintf("-log10(%s)", tstr)
        
        p <- ggplot(res, aes(Fold, -log10(vals))) + 
          geom_point(aes(col = Legend), size = dotSize) + 
          scale_color_manual(values = c("Gain" = "#ED2023", 
                                        "NS" = "grey", 
                                        "Loss" = "#324DA0")) + 
          labs(title = plotTitle, x = xLabel, y = yLabel) + 
          theme_bw() +
          theme(element_text(colour = "black"))
        
        if (bLabels) {
          maxLabels <- min(sidx, maxLabels)
          if (maxLabels > 0 && sidx > 0) {
            xx <- which(idx)[1:maxLabels]
            p <- p + geom_text_repel(data = sigSites[1:maxLabels, 
            ], aes(x = Fold, y = -log10(vals[xx]), 
                   label = rownames(sigSites)[1:maxLabels])) + 
              theme_bw() +
              theme(element_text(colour = "black"))
          }
        }
        plot(p)
      }
    }
  }
  if (bReturnPlot) {
    return(p)
  }
  else {
    if (sidx > 0) {
      return(sigSites[, -10])
    }
    else {
      return(NULL)
    }
  }
})

# Assign the modified function back to the package's namespace
assignInNamespace("pv.DBAplotVolcano", modified_fun, "DiffBind")
#####################################################################

# Plot modified volcano
dba.plotVolcano(ko) # diff_volc_ko, 5 * 5.5

# Plot gain/loss peak comparison
ko.gain <- dba.plotProfile(ko)
dba.plotProfile(ko.gain) # diff_peaks_ko, 4 * 10
