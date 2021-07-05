



############################################################################################################################
############################################################################################################################
############################################################################################################################



rm(list = ls())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




############################################################################################################################
############################################################################################################################
############################################################################################################################


########################################################  libraries ######################################################## 


library(rtracklayer)
library(GenomicFeatures)
library(RColorBrewer)
library(matrixStats)
library(topGO)
library(grid)
library(gridBase)
library(gridExtra)
library(tsTools)
library(BiocParallel)
library(GenomicAlignments)
library(LSD)
library(csaw)
library(pheatmap)
library(Vennerable)
library(csaw)
library(edgeR)






############################################################################################################################
############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################
############################################################################################################################


#####################################################    Load BAM      ###################################################### 


system("mkdir data")


my_dir <- "../BAM_Multi/"

bam_path <- list.files(path = my_dir, pattern = "_pob3.*.bam$|_wt.*.bam$", full.names = TRUE)    

my_chips <- unique(gsub("_.*|.*\\/","",bam_path))
my_chips <- my_chips[my_chips != "IN"]


###################################################################

i=1

for(i in seq_along(my_chips)){
        
        
        data.chip <- windowCounts(bam.files = c(grep(paste0(my_chips[i],"_"), bam_path, value = TRUE),
                                                grep("IN_", bam_path, value = TRUE)),
                                  bin=TRUE, width=250,
                                  param = readParam(minq = 0),
                                  BPPARAM = MulticoreParam(16))
        
        colData(data.chip)$bam.files <- gsub("_[G,A,T,C][G,A,T,C][G,A,T,C].*|.*\\/", "", colData(data.chip)$bam.files)
        colData(data.chip)$bam.files <- gsub("\\-", "_", colData(data.chip)$bam.files)
        colData(data.chip)$bam.files <- gsub("wt1", "wt_1", colData(data.chip)$bam.files)
        colData(data.chip)$bam.files <- gsub("wt2", "wt_2", colData(data.chip)$bam.files)
        
        my_col_data <- data.frame(Reduce(rbind, strsplit(colData(data.chip)$bam.files, "_")), stringsAsFactors = FALSE)
        colnames(my_col_data) <-  c("chip","genotype","replicate")
        rownames(my_col_data) <- 1:nrow(my_col_data)
        
        colData(data.chip) <- cbind(colData(data.chip), my_col_data)
                
        data.chip <- data.chip[seqnames(rowRanges(data.chip)) %in% c("I","II","III")]
        
        my_name <- paste0("data.",my_chips[i])
        assign(my_name, data.chip)
        
        save(list = my_name,file = paste0("data/",my_name,".rda"))
        
}




###################################################################



for(i in seq_along(my_chips)){
        
        
        bin.chip <- windowCounts(bam.files = c(grep(paste0(my_chips[i],"_"), bam_path, value = TRUE),
                                               grep("IN_", bam_path, value = TRUE)),
                                 bin=TRUE, width=10000,
                                 param = readParam(minq = 0),
                                 BPPARAM = MulticoreParam(16))
        
        colData(bin.chip)$bam.files <- gsub("_[G,A,T,C][G,A,T,C][G,A,T,C].*|.*\\/", "", colData(bin.chip)$bam.files)
        colData(bin.chip)$bam.files <- gsub("\\-", "_", colData(bin.chip)$bam.files)
        colData(bin.chip)$bam.files <- gsub("wt1", "wt_1", colData(bin.chip)$bam.files)
        colData(bin.chip)$bam.files <- gsub("wt2", "wt_2", colData(bin.chip)$bam.files)
        
        bin.chip <- bin.chip[seqnames(rowRanges(bin.chip)) %in% c("I","II","III")]
        
        my_name <- paste0("bin.",my_chips[i])
        assign(my_name, bin.chip)
        
        save(list = my_name, 
             file = paste0("data/",my_name,".rda"))
        
}


############################################################################################################################
############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################
############################################################################################################################

