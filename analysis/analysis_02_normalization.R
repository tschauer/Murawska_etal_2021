



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
library(matrixStats)
library(tsTools)
library(BiocParallel)
library(GenomicAlignments)
library(csaw)
library(edgeR)



############################################################################################################################
############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################
############################################################################################################################

##############################################       Load Counts      ###################################################### 



system("mkdir norm_data")



my_data_files <- list.files(path = "data/", pattern = ".rda", full.names = TRUE)


for(i in seq_along(my_data_files)){
        
        load(my_data_files[i]) 
        
}


my_data_ses <- ls(pattern = "^data.")
my_bin_ses <- ls(pattern = "^bin.")


i=1

if(identical(gsub("data.","", my_data_ses),gsub("bin.","", my_bin_ses))){
        
        
        for(i in seq_along(my_data_ses)){
                
                
                data.chip <- get(my_data_ses[i])
                bin.chip <-  get(my_bin_ses[i])
                
                
                if(identical(colnames(data.chip), colnames(bin.chip))){
                        
                        
                        data.chip$norm.factors <- normFactors(bin.chip, se.out = FALSE)
                        
                        
                        assays(data.chip)$adjc <- cpm(asDGEList(data.chip), normalized.lib.sizes = TRUE, log=TRUE)
                        assays(data.chip)$scac <- assays(data.chip)$adjc
                        
                        
                        my_IN_idx <- colData(data.chip)$chip == "IN"
                        my_IP_idx <- colData(data.chip)$chip != "IN"
                        
                        my_IN_assay <- assays(data.chip[,my_IN_idx])$adjc
                        my_IN_means <- colMeans(my_IN_assay)
                        
                        assays(data.chip)$scac[,my_IN_idx] <- t(t(my_IN_assay) - my_IN_means)
                        
                        
                        if(all(identical(colData(data.chip)$genotype[my_IP_idx], 
                                         colData(data.chip)$genotype[my_IN_idx]),
                               identical(colData(data.chip)$replicate[my_IP_idx],
                                         colData(data.chip)$replicate[my_IN_idx]))){
                                
                                
                                assays(data.chip)$norm <- assays(data.chip)$scac
                                
                                assays(data.chip)$norm[,my_IP_idx] <- assays(data.chip)$scac[,my_IP_idx] - assays(data.chip)$scac[,my_IN_idx]
                                assays(data.chip)$norm[,my_IN_idx] <- assays(data.chip)$scac[,my_IN_idx] - assays(data.chip)$scac[,my_IN_idx]
                                
                                
                        }
                        
                        
                        assign(my_data_ses[i], data.chip)
                        
                        save(list = my_data_ses[i], file = paste0("norm_data/", my_data_ses[i],".rda"))
                        
                        
                        rm(list = "data.chip")
                        rm(list = "bin.chip")
                        
                }
        }
}



############################################################################################################################
############################################################################################################################
############################################################################################################################




system("mkdir plots")

pdf(file = paste0("plots/density_plots.pdf"), width = 10, height = 6)

par(mfrow=c(2,4))


for(i in seq_along(my_data_ses)){

        
        data.chip <- get(my_data_ses[i])
        
        for(i in 1:ncol(assays(data.chip)$adjc)){
                
                plot(density(assays(data.chip)$adjc[,i], from = -2, to = 10),
                     main = colData(data.chip)$bam.files[i], xlab = "log2 norm counts")
        }
        
        for(i in 1:ncol(assays(data.chip)$scac)){
                
                plot(density(assays(data.chip)$scac[,i], from = -2, to = 10),
                     main = colData(data.chip)$bam.files[i], xlab = "log2 norm counts")
        }
        
        for(i in 1:ncol(assays(data.chip)$norm)){
                
                plot(density(assays(data.chip)$norm[,i], from = -2, to = 10),
                     main = colData(data.chip)$bam.files[i], xlab = "log2 norm counts")
        }
        
        rm(list = "data.chip")
        
}

dev.off()



############################################################################################################################
############################################################################################################################
############################################################################################################################



system("mkdir coverage_chip")


i=1


for(i in seq_along(my_data_ses)){
        
        
        data.chip <- get(my_data_ses[i])
        
        gr <- rowRanges(data.chip)
        
        for(i in which(colData(data.chip)$chip != "IN")){
                
                gr$score <- assays(data.chip)$norm[,i]     
                
                my_cov <- coverage(gr, weight = "score")
                
                my_name <- paste0("coverage.",colData(data.chip)$bam.files[i])
                
                assign(my_name, my_cov)
                
                save(list = my_name,   file = paste0("coverage_chip/", my_name,".rda"))
                export.bw(get(my_name), con = paste0("coverage_chip/", my_name,".bw"))
                
        }
        
        rm(list = "data.chip")
        
}




############################################################################################################################
############################################################################################################################
############################################################################################################################






