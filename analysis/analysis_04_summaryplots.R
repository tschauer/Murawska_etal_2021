



##################################################################################################################################
##################################################################################################################################




rm(list = ls())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




##################################################################################################################################
##################################################################################################################################

########################################################     libraries    ######################################################## 


library(GenomicFeatures)
library(AnnotationDbi)
library(RColorBrewer)
library(DESeq2)
library(affy)
library(pheatmap)
library(scales)
library(genefilter)
library(ggplot2)
library(sva)
library(GenomicFeatures)
library(rtracklayer)
library(grid)
library(gridExtra)
library(dplyr)
library(LSD)
library(tiff)
library(png)
library(gdata)
library(vioplot)
library(matrixStats)


source("functions.R")



##################################################################################################################################
##################################################################################################################################

######################################################   annotation    ########################################################## 



my_chromInfo <- read.table("/Volumes/PromisePegasus/Projects/Magdalena/genome/Schizosaccharomyces_pombe.ASM294v2.chromInfo.txt")



##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################

#################################################           Peaks            ##################################################### 



peaks.H3K9me2_wt_1 <- import("../Homer_reps_Multi/H3K9me2_wt_1.histone.F2.bed")
peaks.H3K9me2_wt_2 <- import("../Homer_reps_Multi/H3K9me2_wt_2.histone.F2.bed")

peaks.H3K9me2_wt <- intersect(peaks.H3K9me2_wt_1, peaks.H3K9me2_wt_2)






##################################################################################################################################
##################################################################################################################################


peaks.H3K9me2_wt <- peaks.H3K9me2_wt[seqnames(peaks.H3K9me2_wt) %in%  c("I","II")]


my_TEL <- makeGRangesFromDataFrame(data.frame(chr = c("I","I","II","II"),
                                              start = c(0,5479000, 0, 4440000),
                                              end= c(100000, 5579000, 100000, 4539000),
                                              strand = "+"))

my_CEN <- makeGRangesFromDataFrame(data.frame(chr = c("I","II"),
                                              start = c(3720000, 1570000),
                                              end= c(3820000, 1670000),
                                              strand = "+"))


my_genome <- makeGRangesFromDataFrame(data.frame(chr = c("I","II","III"),
                                                 start = c(0, 0, 0),
                                                 end= floor(my_chromInfo$V2[1:3]/1000)*1000,
                                                 strand = "+"
))





my_regions <- GRangesList("TEL" = intersect(peaks.H3K9me2_wt, my_TEL),
                          "CEN" = intersect(peaks.H3K9me2_wt, my_CEN),
                          "EU" =  setdiff(my_genome[seqnames(my_genome) %in% c("I","II")], peaks.H3K9me2_wt))



##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################

#################################################      Normalized data       ##################################################### 



my_data_files <- list.files(path = "norm_data//", pattern = ".rda", full.names = TRUE)


for(i in seq_along(my_data_files)){
        
        load(my_data_files[i]) 
        
}


my_data_ses <- ls(pattern = "^data.")



##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################

#################################################      Setup for plots       ##################################################### 



my_colors <- brewer.pal(9, "Set1")[-6]

# new colors
my_colors[1] <- "#E69F00"
my_colors[2] <- "#0072B2"
my_colors[4] <- "#CC79A7"
my_colors[8] <- "#999999"
my_colors[9] <- "#E41A1C"
my_colors[10] <- "#555555"



my_ylims <- list(c(-4,3),
                 c(-2,2),
                 c(-1,6),
                 c(-2,2),
                 c(-4,6),
                 c(-2,4))

i=3

#################################################      


pdf(paste0("plots/boxplots_reps.pdf"), width = 10, height = 7, useDingbats = FALSE)

par(oma=c(4,4,4,2), mar=c(6,2,2,2), mgp=c(1,0.5,0),
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.25, pch=20)



par(mfcol=c(2,6))

for(i in seq_along(my_data_ses)){
        
        
        data.chip <- get(my_data_ses[i])
        data.chip <- data.chip[,colData(data.chip)$chip != "IN"]
        
        plot_data <- list(TEL = assays(data.chip[rowRanges(data.chip) %over% my_regions$TEL])$norm,
                          CEN = assays(data.chip[rowRanges(data.chip) %over% my_regions$CEN])$norm,
                          EU  = assays(data.chip[rowRanges(data.chip) %over% my_regions$EU ])$norm)
        
        for(repi in 1:2){
                
                my_wt_idx <-  which(colData(data.chip)$genotype == "wt" &   colData(data.chip)$replicate == repi)
                my_mut_idx <- which(colData(data.chip)$genotype == "pob3" & colData(data.chip)$replicate == repi)
                
                boxplot(list(plot_data$TEL[,my_wt_idx],
                             plot_data$TEL[,my_mut_idx],
                             plot_data$CEN[,my_wt_idx],
                             plot_data$CEN[,my_mut_idx],
                             plot_data$EU[,my_wt_idx],
                             plot_data$EU[,my_mut_idx]),
                        las=1, main = paste0(gsub("data.", "", my_data_ses[i]),"_", repi),
                        ylim = my_ylims[[i]],
                        col = my_colors[c(8,2)], xaxt="n", cex = 0.2)  
                
                axis(side = 1, at = seq(1.5, 2*length(plot_data),2), 
                     labels = names(plot_data), las=2)
                
                abline(v=c(8.5,16.5), lty=1)
                
                # if(i == 2){
                #         legend("top", legend = c(expression(bold("WT")),expression(paste(bolditalic("pob3"), Delta))), 
                #                fill  = my_colors[c(8,2)], bg = "white")
                # }
                # 
        }
        
        rm(list = "data.chip")
        
}

mtext("log2(ChIP/Input)", side = 2, line = 1, at = 0.5, adj = 0.5,  cex = 1.2, font = 1, outer = TRUE)

par(fig=c(0,1,0,1), oma = c(0,0,0,0), mar = c(1,0,0,0), new = TRUE)

plot.new()

legend("bottom", legend = c(expression(bold("WT")),expression(paste(bolditalic("pob3"), Delta))), 
       fill  = my_colors[c(8,2)], bg = "white", horiz = TRUE, pt.cex=1.2)


dev.off()



##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################






pdf(paste0("plots/boxplots_average.pdf"), width = 12, height = 4, useDingbats = FALSE)

par(oma=c(4,4,4,2), mar=c(4,2,2,2), mgp=c(1,0.5,0),
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.25, pch=20)

par(mfrow=c(1,6))

for(i in seq_along(my_data_ses)){
        
        
        data.chip <- get(my_data_ses[i])
        data.chip <- data.chip[,colData(data.chip)$chip != "IN"]
        
        plot_data <- list(TEL = assays(data.chip[rowRanges(data.chip) %over% my_regions$TEL])$norm,
                          CEN = assays(data.chip[rowRanges(data.chip) %over% my_regions$CEN])$norm,
                          EU  = assays(data.chip[rowRanges(data.chip) %over% my_regions$EU ])$norm)
        
        my_wt_idx <-  which(colData(data.chip)$genotype == "wt")
        my_mut_idx <- which(colData(data.chip)$genotype == "pob3")
        
        boxplot(list(rowMeans(plot_data$TEL[,my_wt_idx]),
                     rowMeans(plot_data$TEL[,my_mut_idx]),
                     rowMeans(plot_data$CEN[,my_wt_idx]),
                     rowMeans(plot_data$CEN[,my_mut_idx]),
                     rowMeans(plot_data$EU[,my_wt_idx]),
                     rowMeans(plot_data$EU[,my_mut_idx])),
                las=1, main = paste0(gsub("data.", "", my_data_ses[i])),
                ylim = my_ylims[[i]],
                col = my_colors[c(8,2)], xaxt="n", cex = 0.2)  
        
        axis(side = 1, at = seq(1.5, 2*length(plot_data),2), 
             labels = names(plot_data), las=2)
        
        abline(v=c(8.5,16.5), lty=1)
        
        # if(i == 2){
        #         legend("top", legend = c(expression(bold("WT")),expression(paste(bolditalic("pob3"), Delta))), 
        #                fill  = my_colors[c(8,2)], bg = "white")
        # }
        
        rm(list = "data.chip")
        
        if(i == 1){
                mtext("log2(ChIP/Input)", side = 2, line = 2.5,  adj = 0.5,  cex = 1.2, font = 1)
                
        }
        
        
}

par(fig=c(0,1,0,1), oma = c(0,0,0,0), mar = c(1,0,0,0), new = TRUE)

plot.new()

legend("bottom", legend = c(expression(bold("WT")),expression(paste(bolditalic("pob3"), Delta))), 
       fill  = my_colors[c(8,2)], bg = "white", horiz = TRUE, pt.cex=1.2)


dev.off()



##################################################################################################################################
##################################################################################################################################






##################################################################################################################################
##################################################################################################################################





pdf(paste0("plots/boxplots_average_wtonly.pdf"), width = 10, height = 4, useDingbats = FALSE)

par(oma=c(4,4,4,2), mar=c(4,2,2,2), mgp=c(1,0.5,0),
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.25, pch=20)

par(mfrow=c(1,6))

for(i in seq_along(my_data_ses)){
        
        
        data.chip <- get(my_data_ses[i])
        data.chip <- data.chip[,colData(data.chip)$chip != "IN"]
        
        plot_data <- list(TEL = assays(data.chip[rowRanges(data.chip) %over% my_regions$TEL])$norm,
                          CEN = assays(data.chip[rowRanges(data.chip) %over% my_regions$CEN])$norm,
                          EU  = assays(data.chip[rowRanges(data.chip) %over% my_regions$EU ])$norm)
        
        my_wt_idx <-  which(colData(data.chip)$genotype == "wt")
        
        boxplot(list(rowMeans(plot_data$TEL[,my_wt_idx]),
                     rowMeans(plot_data$CEN[,my_wt_idx]),
                     rowMeans(plot_data$EU[,my_wt_idx])),
                las=1, main = paste0(gsub("data.", "", my_data_ses[i])),
                ylim = my_ylims[[i]],
                col = my_colors[c(8)], xaxt="n", cex = 0.2)  
        
        axis(side = 1, at = 1:3, 
             labels = names(plot_data), las=2)
        
        abline(v=c(8.5,16.5), lty=1)
        
        # if(i == 2){
        #         legend("top", legend = c(expression(bold("WT"))), 
        #                fill  = my_colors[c(8)], bg = "white")
        # }
        
        rm(list = "data.chip")
        
        if(i == 1){
                mtext("log2(ChIP/Input)", side = 2, line = 2.5,  adj = 0.5,  cex = 1.2, font = 1)
                
        }
        
        
}

par(fig=c(0,1,0,1), oma = c(0,0,0,0), mar = c(1,0,0,0), new = TRUE)

plot.new()

legend("bottom", legend = c(expression(bold("WT"))), 
       fill  = my_colors[c(8)], bg = "white", horiz = TRUE, pt.cex=1.2)



dev.off()



##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################



my_ylims <- list(c(-4,1),
                 c(-1,1),
                 c(-1,7),
                 c(-2,2),
                 c(-3,1),
                 c(-1,0.5))




library(lme4)
library(lmerTest)



pdf(paste0("plots/dotplots_reps.pdf"), width = 12, height = 4, useDingbats = FALSE)

par(oma=c(4,4,4,2), mar=c(4,2,2,2), mgp=c(1,0.5,0),
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.25, pch=20)



par(mfrow=c(1,6))


for(i in seq_along(my_data_ses)){
        
        
        data.chip <- get(my_data_ses[i])
        data.chip <- data.chip[,colData(data.chip)$chip != "IN"]
        
        plot_data <- list(TEL = assays(data.chip[rowRanges(data.chip) %over% my_regions$TEL])$norm,
                          CEN = assays(data.chip[rowRanges(data.chip) %over% my_regions$CEN])$norm,
                          EU  = assays(data.chip[rowRanges(data.chip) %over% my_regions$EU ])$norm)
        
        plot_medians <- unlist(lapply(plot_data, colMedians))
        
        my_df <- data.frame(region = factor(gsub("[1-9]","", names(plot_medians))),
                            sample = factor(gsub("[A-Z]","", names(plot_medians))),
                            medians = plot_medians)
        
        my_df$genotype <- factor(colData(data.chip)$genotype[my_df$sample], levels = c("wt", "pob3"))
        
        my_df$group <- factor(paste(my_df$region, my_df$genotype, sep="_"), 
                              levels = c("TEL_wt","TEL_pob3","CEN_wt","CEN_pob3","EU_wt","EU_pob3"))
        
        
        
        plot(as.numeric(my_df$group)+c(-0.1,+0.1),
             my_df$medians, 
             main = gsub("data.", "", my_data_ses[i]),
             xaxt="n", xlim = c(0.5,6.5), ylim=my_ylims[[i]], 
             pch=20, cex = 1.1, xlab="", ylab="",
             col = my_colors[c(8,2)][my_df$genotype])
        
        points(as.numeric(my_df$group)+c(-0.1,+0.1),
               my_df$medians,
               pch=1, cex = 1.1, col = "#333333", lwd=0.5)
        
        my_labels <- gsub("_.*","",levels(my_df$group))[seq(1,5,2)]
        
        axis(side = 1, at = seq(1.5, 6, 2), labels = my_labels, las=2)
        abline(v=c(2.5,4.5), lty=2)
        
        
        pvals <- sapply(my_labels, function(x){
                my_df$region <- relevel(my_df$region, x)
                coef(summary(lmer(medians ~ genotype * region + (1|sample), data = my_df)))["genotypepob3","Pr(>|t|)"]
        })
        
        
        text(seq(1, 6, 2), my_ylims[[i]][2], adj = c(0.1,1.1),
             labels = paste0("<",ceiling(pvals*100)/100))
        
        if(i == 1){
                mtext("Median log2(ChIP/Input)", side = 2, line = 2.5, cex = 1.2, font = 1)
                
        }
        
        rm(list = "data.chip")
        
        
}  


par(fig=c(0,1,0,1), oma = c(0,0,0,0), mar = c(1,0,0,0), new = TRUE)

plot.new()

legend("bottom", legend = c(expression(bold("WT")),expression(paste(bolditalic("pob3"), Delta))), 
       col  = my_colors[c(8,2)], bg = "white", horiz = TRUE, pch=19, pt.cex=1.2)

legend("bottom", legend = c(expression(bold("WT")),expression(paste(bolditalic("pob3"), Delta))), 
       col  = "#333333", horiz = TRUE, pch=1, pt.cex = 1.2, pt.lwd = 0.5)

dev.off()

##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################




