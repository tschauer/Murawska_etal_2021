


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


gtf <- rtracklayer::import("/Volumes/PromisePegasus/Projects/Magdalena/genome/Schizosaccharomyces_pombe.ASM294v2.37.gtf")


my_genes <- gtf[gtf$type == "gene"]
mcols(my_genes) <- mcols(my_genes)[c(5,6,8)]


my_genes_protein_coding <- my_genes[my_genes$gene_biotype == "protein_coding"]
my_genes_protein_coding <- my_genes_protein_coding[seqnames(my_genes_protein_coding) %in% c("I","II","III")]


my_exons <- gtf[gtf$type == "exon" & gtf$gene_biotype == "protein_coding"]
mcols(my_exons) <- mcols(my_exons)[c(5,6,8)]




my_genes_ncRNA <- gtf[gtf$type == "gene" & gtf$gene_biotype == "ncRNA"]
my_genes_ncRNA <- my_genes_ncRNA[seqnames(my_genes_ncRNA) %in% c("I","II","III")]

my_exons_ncRNA <- gtf[gtf$type == "exon" & gtf$gene_biotype == "ncRNA"]
mcols(my_exons_ncRNA) <- mcols(my_exons_ncRNA)[c(5,6,8)]




my_genes_rRNA <- gtf[gtf$type == "gene" & gtf$gene_biotype == "rRNA"]
my_genes_rRNA <- my_genes_rRNA[seqnames(my_genes_rRNA) %in% c("I","II","III")]

my_exons_rRNA <- gtf[gtf$type == "exon" & gtf$gene_biotype == "rRNA"]
mcols(my_exons_rRNA) <- mcols(my_exons_rRNA)[c(5,6,8)]




my_genes_tRNA <- gtf[gtf$type == "gene" & gtf$gene_biotype == "tRNA"]
my_genes_tRNA <- my_genes_tRNA[seqnames(my_genes_tRNA) %in% c("I","II","III")]

my_exons_tRNA <- gtf[gtf$type == "exon" & gtf$gene_biotype == "tRNA"]
mcols(my_exons_tRNA) <- mcols(my_exons_tRNA)[c(5,6,8)]



my_genes_pseudogene <- gtf[gtf$type == "gene" & gtf$gene_biotype == "pseudogene"]
my_genes_pseudogene <- my_genes_pseudogene[seqnames(my_genes_pseudogene) %in% c("I","II","III")]

my_exons_pseudogene <- gtf[gtf$type == "exon" & gtf$gene_biotype == "pseudogene"]
mcols(my_exons_pseudogene) <- mcols(my_exons_pseudogene)[c(5,6,8)]


##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

#################################################      ChIP Coverage       ##################################################### 



my_covs <- file.path("coverage_chip/", list.files("coverage_chip/", pattern = "coverage.*.rda$"))


for(i in seq_along(my_covs)){
        
        load(my_covs[i])
}


my_covs <- ls(pattern = "coverage\\.")


##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################



#################################################     Average reps     ##################################################### 




my_covs_rep <- my_covs[grep("_1$", my_covs)]

sami=3


for(sami in seq_along(my_covs_rep)){
        
        my_name <- gsub("_1","",my_covs_rep[sami])
        
        my_cov1 <- get(paste0(my_name, "_1"))  
        my_cov2 <- get(paste0(my_name, "_2"))  
        
        my_cov3 <- (my_cov1 + my_cov2)/2
        
        assign(paste0(my_name, "_3"), my_cov3)
}





##################################################################################################################################
##################################################################################################################################














##################################################################################################################################
##################################################################################################################################

#################################################     Peaks     ##################################################### 




peaks.H3K9me2_wt_1 <- import("../Homer_reps_Multi/H3K9me2_wt_1.histone.F2.bed")
peaks.H3K9me2_wt_2 <- import("../Homer_reps_Multi/H3K9me2_wt_2.histone.F2.bed")

peaks.H3K9me2_wt <- intersect(peaks.H3K9me2_wt_1, peaks.H3K9me2_wt_2)




##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################



system("mkdir coverage_plots")




my_colors <- brewer.pal(9, "Set1")[-6]

# new colors
my_colors[1] <- "#E69F00"
my_colors[2] <- "#0072B2"
my_colors[4] <- "#CC79A7"
my_colors[8] <- "#999999"
my_colors[9] <- "#E41A1C"
my_colors[10] <- "#555555"

#################################################     coverage plots old style     ##################################################### 


source("functions.R")


my_regions <- list(
        chrI_TELL = c("I", 0, 102000),
        chrI_CEN = c("I", 3720000, 3820000),
        chrI_TELR = c("I", 5479000, 5579000),
        chrII_TELL = c("II", 0, 102000), 
        chrII_CEN = c("II", 1570000,1670000),
        chrII_TELR = c("II", 4440000, 4540000),
        chrIII_TELL = c("III", 0, 102000),
        chrIII_CEN =  c("III", 1055000, 1155000),
        chrIII_TELR = c("III", 2352000, 2452000)
)


my_region <- my_regions[[1]]

my_bin_size = 250
my_xticks = 20*10^3


regi=1
sami=5
repi = 1
ordi = 1

my_ylims <- list(c(-4,3),
                 c(-3,2),
                 c(-1,7),
                 c(-1,3),
                 c(-4,4),
                 c(-1,3))

my_covs <- ls(pattern = "coverage\\.")


for(repi in 1:3){
        
        
        
        my_covs_rep <- my_covs[grep(paste0("_", repi, "$"), my_covs)]
        my_covs_rep_wt <-  my_covs_rep[grep("_wt", my_covs_rep)]
        my_covs_rep_mut <-  my_covs_rep[grep("_pob3", my_covs_rep)]
        
        if(!identical(gsub("_.*","", my_covs_rep_wt), gsub("_.*","", my_covs_rep_mut))){next()}
        
        for(ordi in 1:2){
                
                if(repi == 3){
                        pdf(paste0("coverage_plots/coverage_aveplot",ordi,".oldstyle.pdf"), width = 20, height = 11.69, useDingbats = FALSE)
                        
                } else {
                        pdf(paste0("coverage_plots/coverage_rep",repi, "plot",ordi,".oldstyle.pdf"), width = 20, height = 11.69, useDingbats = FALSE)
                }
                
                par(oma=c(2,2,2,2), mar=c(0,0,0,0), mgp=c(1,0.5,0),
                    cex.axis = 1, cex.main = 1.25, cex.lab=1.25, pch=19)
                
                
                for(sami in seq_along(my_covs_rep_wt)){
                        
                        for(regi in seq_along(my_regions)){
                                
                                
                                my_region <- my_regions[[regi]]
                                
                                
                                if(ordi == 1){
                                        
                                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), (0.125+(sami-1)/7), (0.125+(sami-1)/7+0.12)), new = TRUE)
                                        
                                        plotCoverage(my_coverage = get(my_covs_rep_wt[sami]), 
                                                     my_region = my_region, 
                                                     my_color = my_colors[8], 
                                                     my_title = "",  
                                                     my_bin_size = my_bin_size,
                                                     my_ylims = my_ylims[[sami]], 
                                                     log_scale = FALSE, cex_yaxis = 1,
                                                     smoother = 11, line_lwd = 2,
                                                     my_peaks = peaks.H3K9me2_wt)
                                        
                                        title(main = gsub(".*\\.|_.*", "",(my_covs_rep_mut[sami])), line = 0)
                                        
                                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), (0.125+(sami-1)/7), (0.125+(sami-1)/7+0.12)), new = TRUE)
                                        
                                        plotCoverage(my_coverage = get(my_covs_rep_mut[sami]), 
                                                     my_region = my_region, 
                                                     my_color = my_colors[2], 
                                                     my_title = "",  
                                                     my_bin_size = my_bin_size,
                                                     my_ylims = my_ylims[[sami]], 
                                                     log_scale = FALSE,
                                                     smoother = 11, 
                                                     y_axis = F, line_lwd = 2)
                                } else {
                                        
                                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), (0.125+(sami-1)/7), (0.125+(sami-1)/7+0.12)), new = TRUE)
                                        
                                        plotCoverage(my_coverage = get(my_covs_rep_mut[sami]), 
                                                     my_region = my_region, 
                                                     my_color = my_colors[2], 
                                                     my_title = "",  
                                                     my_bin_size = my_bin_size,
                                                     my_ylims = my_ylims[[sami]], 
                                                     log_scale = FALSE, cex_yaxis = 1,
                                                     smoother = 11, line_lwd = 2,
                                                     my_peaks = peaks.H3K9me2_wt)
                                        
                                        title(main = gsub(".*\\.|_.*", "",(my_covs_rep_mut[sami])), line = 0)
                                        
                                        
                                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), (0.125+(sami-1)/7), (0.125+(sami-1)/7+0.12)), new = TRUE)
                                        
                                        plotCoverage(my_coverage = get(my_covs_rep_wt[sami]), 
                                                     my_region = my_region, 
                                                     my_color = my_colors[8], 
                                                     my_title = "",  
                                                     my_bin_size = my_bin_size,
                                                     my_ylims = my_ylims[[sami]], 
                                                     log_scale = FALSE,
                                                     smoother = 11,
                                                     y_axis = F, line_lwd = 2)
                                }
                                
                        }
                        
                }
                
                for(regi in seq_along(my_regions)){
                        
                        my_region <- my_regions[[regi]]
                        
                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), 0.030, 0.120), new = TRUE)
                        plotAnnotation(my_genes = my_genes_protein_coding, my_exons = my_exons, my_region = my_region, my_xticks = my_xticks, gene_name_text = FALSE, rect_color = my_colors[10], fill_color = my_colors[10], cex_xaxis = 1)
                        
                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), 0.030, 0.120), new = TRUE)
                        plotAnnotation(my_genes = my_genes_ncRNA, my_exons = my_exons_ncRNA, my_region = my_region, my_xticks = NULL, gene_name_text = FALSE, rect_color = my_colors[1], fill_color = my_colors[1])
                        
                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), 0.030, 0.120), new = TRUE)
                        plotAnnotation(my_genes = my_genes_rRNA, my_exons = my_exons_rRNA, my_region = my_region, my_xticks = NULL, gene_name_text = FALSE, rect_color = my_colors[4], fill_color = my_colors[4])
                        
                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), 0.030, 0.120), new = TRUE)
                        plotAnnotation(my_genes = my_genes_tRNA, my_exons = my_exons_tRNA, my_region = my_region, my_xticks = NULL, gene_name_text = FALSE, rect_color = my_colors[3], fill_color = my_colors[3])
                        
                        par(fig = c((0+(regi-1)/9), (0+(regi-1)/9+0.1111), 0.030, 0.120), new = TRUE)
                        plotAnnotation(my_genes = my_genes_pseudogene, my_exons = my_exons_pseudogene, my_region = my_region,  my_xticks = NULL, gene_name_text = FALSE, rect_color = my_colors[6], fill_color = my_colors[6])
                        
                }
                
                mtext(c("chrI","chrII", "chrIII"), side = 1, las=1, line = -2.5, at = c(0.200,0.500,0.800),  cex = 1.2, font = 1, outer = TRUE)
                mtext("log2(ChIP/Input)", side = 2, line = -0.5, at = 0.5, adj = 0.475,  cex = 1.2, font = 1, outer = TRUE)
                
                if(repi == 3){
                        mtext("average", side = 3, las=1, line = 0, at = 0.5,  cex = 1.2, font = 1, outer = TRUE)
                } else {
                        mtext(paste("replicate #", repi), side = 3, las=1, line = 0, at = 0.5,  cex = 1.2, font = 1, outer = TRUE)
                }
                
                par(oma = c(0,0,0,0), mar = c(1,1,1,1), fig = c(0,0.20,0.75,1.0), new = TRUE)
                
                plot.new()
                
                legend("topleft", legend = c(expression(bold("WT")),expression(paste(bolditalic("pob3"), Delta))), 
                       col = my_colors[c(8,2)], lwd = 5)
                
                par(oma = c(0,0,0,0), mar = c(0.5,1,1,1), fig = c(0,1.00,0.00,0.25), new = TRUE)
                
                plot.new()
                
                legend("bottom", legend = c("mRNA", "ncRNA","rRNA","tRNA","pseudo"), 
                       col = my_colors[c(10,1,4,3,6)], lwd = 5, horiz = T)
                
                dev.off()
                
        }
        
}


##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################



my_colors <- brewer.pal(9, "Set1")[-6]

# new colors
my_colors[1] <- "#E69F00"
my_colors[2] <- "#0072B2"
my_colors[4] <- "#CC79A7"
my_colors[8] <- "#999999"
my_colors[9] <- "#E41A1C"
my_colors[10] <- "#555555"

#################################################     coverage plots new style     ##################################################### 


source("functions.R")


my_regions <- list(
        chrI_TELL = c("I", 0, 105000),
        chrI_CEN = c("I", 3720000, 3825000),
        chrI_TELR = c("I", 5479000, 5580000),
        chrII_TELL = c("II", 0, 105000), 
        chrII_CEN = c("II", 1570000,1675000),
        chrII_TELR = c("II", 4440000, 4542000),
        chrIII_TELL = c("III", 0, 105000),
        chrIII_CEN =  c("III", 1055000, 1160000),
        chrIII_TELR = c("III", 2352000, 2454000),
        chrIII_rRNA_L  = c("III", 0,27000),
        chrIII_rRNA_R  = c("III",2432000,2453000),
        mat = c("II", 2100000, 2160000),
        mcp7 = c("I", 576000, 583000),
        mei4 = c("II", 1471000, 1478000)
)


my_region <- my_regions[[2]]

my_bin_size = 250
my_xticks = 20*10^3

n_of_regions <- length(my_regions)
size_of_snapshot <- 1/n_of_regions

regi=1
sami=1
repi = 1
ordi = 1




for(repi in 1:3){
        
        
        
        my_covs_rep <- my_covs[grep(paste0("_", repi, "$"), my_covs)]
        my_covs_rep_wt <-  my_covs_rep[grep("_wt", my_covs_rep)]
        my_covs_rep_mut <-  my_covs_rep[grep("_pob3", my_covs_rep)]
        
        if(!identical(gsub("_.*","", my_covs_rep_wt), gsub("_.*","", my_covs_rep_mut))){next()}
        
        for(ordi in 1:2){
                
                
                if(repi == 3){
                        pdf(paste0("coverage_plots/coverage_aveplot",ordi,".newstyle.pdf"), width = 9, height = 7, useDingbats = FALSE)
                        
                } else {
                        pdf(paste0("coverage_plots/coverage_rep",repi, "plot",ordi,".newstyle.pdf"), width = 9, height = 7, useDingbats = FALSE)
                }    
                
                par(oma=c(3,3,3,3), mar=c(0,0,0,0), mgp=c(2,1.5,0.5),
                    cex.axis = 1, cex.main = 1.25, cex.lab=1.25, pch=19, lwd=2)
                
                
                for(sami in seq_along(my_covs_rep_wt)){
                        
                        for(regi in seq_along(my_regions)){
                                
                                my_region <- my_regions[[regi]]
                                
                                # if(grepl("TEL|CEN", names(my_regions[regi]))){
                                #         my_bin_size <- ifelse(grepl("H3_", my_covs_rep_wt[sami]), 500, 500) 
                                # } else if(grepl("rRNA", names(my_regions[regi]))){
                                #         my_bin_size <- 50
                                # } else {
                                #         my_bin_size <- 30
                                # }
                                
                                if(ordi == 1){
                                        
                                        par(fig = c(0,1,0.25,0.75))
                                        
                                        plotCoverage(my_coverage = get(my_covs_rep_wt[sami]), 
                                                     my_region = my_region, 
                                                     my_color = my_colors[8], 
                                                     my_title = "",  
                                                     my_bin_size = my_bin_size, 
                                                     smoother = 5, line_lwd = 3,
                                                     my_ylims = my_ylims[[sami]], 
                                                     log_scale = FALSE,
                                                     my_peaks = peaks.H3K9me2_wt)
                                        
                                        par(fig = c(0,1,0.25,0.75), new = TRUE)
                                        
                                        plotCoverage(my_coverage = get(my_covs_rep_mut[sami]), 
                                                     my_region = my_region, 
                                                     my_color = my_colors[2], 
                                                     my_title = "",  
                                                     my_bin_size = my_bin_size,
                                                     smoother = 5,  line_lwd = 3,
                                                     my_ylims = my_ylims[[sami]], 
                                                     log_scale = FALSE, 
                                                     y_axis = F)
                                } else {
                                        par(fig = c(0,1,0.25,0.75))
                                        
                                        plotCoverage(my_coverage = get(my_covs_rep_mut[sami]), 
                                                     my_region = my_region, 
                                                     my_color = my_colors[2], 
                                                     my_title = "",  
                                                     my_bin_size = my_bin_size,
                                                     smoother = 5,  line_lwd = 3,
                                                     my_ylims = my_ylims[[sami]], 
                                                     log_scale = FALSE)
                                        
                                        
                                        par(fig = c(0,1,0.25,0.75), new = TRUE)
                                        
                                        plotCoverage(my_coverage = get(my_covs_rep_wt[sami]), 
                                                     my_region = my_region, 
                                                     my_color = my_colors[8], 
                                                     my_title = "",  
                                                     my_bin_size = my_bin_size,
                                                     smoother = 5, line_lwd = 3,
                                                     my_ylims = my_ylims[[sami]], 
                                                     log_scale = FALSE, 
                                                     y_axis = F,
                                                     my_peaks = peaks.H3K9me2_wt)
                                        
                                }
                                
                                par(fig = c(0,1,0,0.9), new = TRUE)
                                
                                plot.new()
                                
                                legend("topright", legend = c(expression(bold("WT")),expression(paste(bolditalic("pob3"), Delta))),
                                       title = "Genotype", bty="n", cex=1.25,
                                       col = my_colors[c(8,2)], lwd = 5)
                                
                                
                                legend("bottom", legend = c("mRNA", "ncRNA","rRNA","tRNA","pseudo"), 
                                       text.width = c(0.06,0.06,0.06,0.06,0.06), 
                                       fill =  my_colors[c(10,1,4,3,6)], border = NA, bty="n", horiz = T)
                                
                                
                                if(grepl("TEL|CEN", names(my_regions[regi]))){
                                        my_xticks <- 2e4
                                } else if(grepl("rRNA|mat", names(my_regions[regi]))){
                                        my_xticks <- 5e3
                                } else {
                                        my_xticks <- 2e3
                                }
                                
                                
                                par(fig = c(0, 1, 0.0725, 0.250), new = TRUE)
                                plotAnnotation(my_genes = my_genes_protein_coding, my_exons = my_exons, my_region = my_region, my_xticks = my_xticks, gene_name_text = FALSE, rect_color = my_colors[10], my_bin_size = my_bin_size, lwd_exon = 3)
                                
                                par(fig = c(0, 1, 0.0725, 0.250), new = TRUE)
                                plotAnnotation(my_genes = my_genes_ncRNA, my_exons = my_exons_ncRNA, my_region = my_region, my_xticks = NULL, gene_name_text = FALSE, rect_color = my_colors[1], my_bin_size = my_bin_size, lwd_exon = 3)
                                
                                par(fig = c(0, 1, 0.0725, 0.250), new = TRUE)
                                plotAnnotation(my_genes = my_genes_rRNA, my_exons = my_exons_rRNA, my_region = my_region, my_xticks = NULL, gene_name_text = FALSE, rect_color = my_colors[4], my_bin_size = my_bin_size, lwd_exon = 3)
                                
                                par(fig = c(0, 1, 0.0725, 0.250), new = TRUE)
                                plotAnnotation(my_genes = my_genes_tRNA, my_exons = my_exons_tRNA, my_region = my_region, my_xticks = NULL, gene_name_text = FALSE, rect_color = my_colors[3], my_bin_size = my_bin_size, lwd_exon = 3)
                                
                                par(fig = c(0, 1, 0.0725, 0.250), new = TRUE)
                                plotAnnotation(my_genes = my_genes_pseudogene, my_exons = my_exons_pseudogene, my_region = my_region, my_xticks = NULL, gene_name_text = FALSE, rect_color = my_colors[6], my_bin_size = my_bin_size, lwd_exon = 3)
                                
                                
                                
                                mtext(text = paste(gsub(".*\\.|_.*", "",(my_covs_rep_wt[sami])),"ChIP"), side = 3, line = -6, cex = 1.5, font = 2, outer = TRUE)
                                
                                mtext(text = names(my_regions[regi]), side = 1, line = -3.25, cex = 1.25, font = 2, outer = TRUE)
                                
                                if(repi == 3){
                                        mtext("average", side = 3, las=1, line = 0, at = 0.5,  cex = 1.2, font = 1, outer = TRUE)
                                } else {
                                        mtext(paste("replicate #", repi), side = 3, las=1, line = 0, at = 0.5,  cex = 1.2, font = 1, outer = TRUE)
                                }
                                
                                mtext("log2(ChIP/Input)", side = 2, line = 0, at = 0.5, adj = 0.475,  cex = 1.5, font = 1, outer = TRUE)
                                
                        }
                }
                
                dev.off()
        }
}



##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################



