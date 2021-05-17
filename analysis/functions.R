



library(tsTools)
library(ShortRead)
library(IRanges)




#####################################################################################################################
#####################################################################################################################



normalize.cov <- function(cov){
        total.cov <- sum(as.numeric(unlist(cov)))
        norm_cov <- round(((cov * 1000000) / total.cov), 4)
        return(norm_cov)
}





#####################################################################################################################
#####################################################################################################################



bedsub2dyad <- function (file.id, type = c("SINGLE", "PAIRED"), width = 1, length.quantile = 0.25) 
{
        if (!requireNamespace("data.table", quietly = TRUE)) {
                stop("library 'data.table' is needed for this function to work. Please install it.", 
                     call. = FALSE)
        }
        options(scipen = 999)
        shift = 73
        bed <- data.frame(data.table::fread(file.id, showProgress = T))
        if (type == "SINGLE") {
                bedR <- GenomicRanges::GRanges(bed[, 1], IRanges::IRanges(bed[,2], bed[, 3]), strand = bed[, 4])
                # bug in tstools esize -> resize
                dyads <- GenomicRanges::resize(bedR, width)
                amount <- strand(dyads)
                runValue(amount) <- ifelse(runValue(amount) == "-", -shift, shift)
                dyads <- GenomicRanges::shift(dyads, as.vector(amount))
        } else {
                bedR <- GenomicRanges::GRanges(bed[, 1], IRanges::IRanges(bed[,2], bed[, 3]))
                bedR <- bedR[width(bedR) <= quantile(width(bedR), prob = length.quantile)]
                dyads <- GenomicRanges::resize(bedR, width, fix = "center")
        }
        cov <- GenomicRanges::coverage(dyads)
        cov
}

#####################################################################################################################
#####################################################################################################################


bam2cov <- function(bam_file, type = c("SINGLE","PAIRED"), stranded = FALSE, width = 150, MAPQ = 255){
        
        if(type == "SINGLE"){
                
                my_bam <- readGAlignments(bam_file, param = ScanBamParam(mapqFilter = MAPQ))
                grs <- as(my_bam, "GRanges")
                grsr <- resize(grs, width)
                
                if(stranded == TRUE){
                        # swap sign
                        covn <- coverage(grsr[strand(grsr) == "+"])
                        covp <- coverage(grsr[strand(grsr) == "-"])
                        return(list(plus = covp, min = covn))
                        
                } else {
                        cov <- coverage(grsr)
                        return(cov)
                }
                
        } else if(type == "PAIRED"){
                
                my_bam <- readGAlignmentPairs(bam_file,param = ScanBamParam(mapqFilter = MAPQ))
                grs <- as(my_bam, "GRanges")
                cov <- coverage(grs)
                
                if(stranded == TRUE){
                        # swap sign
                        covn <- coverage(grsr[strand(grsr) == "+"])
                        covp <- coverage(grsr[strand(grsr) == "-"])
                        return(list(plus = covp, min = covn))
                        
                } else {
                        cov <- coverage(grsr)
                        return(cov)
                }
                
        }
        
}



#####################################################################################################################





#####################################################################################################################
#####################################################################################################################





#####################################################################################################################


vector.resizing <- function(x,final.len){
        y <- vector()
        len <- length(x)
        y <-spline(1:len,x,n=final.len)$y
        return(y)
}


#####################################################################################################################








coverageGeneBodyScaled <- 
        
        function( 
                coverage, 
                coords, 
                margin_outer = 750,
                margin_inner = 500, 
                scale=1000,
                chr.include
        )
        {
                
                coverage <- coverage[ names(coverage) %in% chr.include]
                
                coords <- coords[coords$chr %in% names(coverage), ]
                
                
                result <- lapply(names(coverage), function(x) {
                        
                        my.cov <- coverage[[x]]
                        my.coords <- coords[coords$chr == x, ]
                        
                        mw.views.up <- IRanges::Views(my.cov, start = my.coords$start - (margin_outer-1), end = my.coords$start + (margin_inner))
                        
                        mw.views.body <- IRanges::Views(my.cov, start = my.coords$start + (margin_inner), my.coords$end - (margin_inner))
                        
                        mw.views.down <- IRanges::Views(my.cov, start = my.coords$end - (margin_inner), end = my.coords$end + (margin_outer-1))
                        
                        flt1 <- start(mw.views.up) > 0 &  end(mw.views.down) < length(my.cov)
                        
                        mw.views.up <-   mw.views.up[flt1, ]
                        mw.views.body <- mw.views.body[flt1, ]
                        mw.views.down <- mw.views.down[flt1, ]
                        
                        my.coords <- my.coords[flt1, ]
                        
                        flt2 <- 
                                !(is.na(mean(mw.views.up))) &
                                !(is.na(mean(mw.views.body))) &
                                !(is.na(mean(mw.views.down)))
                        
                        mw.views.up <-   mw.views.up[flt2, ]
                        mw.views.body <- mw.views.body[flt2, ]
                        mw.views.down <- mw.views.down[flt2, ]
                        
                        my.coords <- my.coords[flt2, ]
                        
                        if (length(mw.views.body) > 0) {
                                
                                
                                mat.up <- matrix(unlist(lapply(mw.views.up, as.numeric)), nrow = length(mw.views.up), byrow = TRUE)
                                mat.down <- matrix(unlist(lapply(mw.views.down, as.numeric)), nrow = length(mw.views.down), byrow = TRUE)
                                
                                mat.body <- matrix(nrow = length(mw.views.body), ncol = scale)
                                
                                for(i in seq_along(mw.views.body)){
                                        mat.body[i,] <- vector.resizing(as.vector(mw.views.body[[i]]) ,scale)
                                }
                                
                                mat <- cbind(mat.up,mat.body,mat.down)
                                rownames(mat) <- rownames(my.coords)
                                
                                return(mat)
                        }
                        else {
                                return(NULL)
                        }
                })
                
                mat <- Reduce(rbind, result)
                coords <- coords[rownames(coords) %in% rownames(mat),]
                o <- match(rownames(coords), rownames(mat))
                mat <- mat[o, ]
                mat[coords$strand == "-", ] <- t(apply(mat[coords$strand == "-", ], 1, rev))
                #colnames(mat) <- seq(-ceiling(window.size/2), ceiling(window.size/2))
                mat
        }









#####################################################################################################################




#####################################################################################################################
#####################################################################################################################




#####################################################################################################################

# sense_mat <-  log2Diff.sense.204_201
# antis_mat <-  log2Diff.antis.204_201
# anno_df <- as.data.frame(my_genes)
# gene_list <- rownames(my_res_top200)
# color_limits = c(-3,3)
# title = "204 vs. 201"
# split = rep(1, length(gene_list))
# names(hm_split) <- gene_list


asHeatmap <- function(title,
                      group_name,
                      sense_mat, 
                      antis_mat, 
                      anno_df,
                      gene_list,
                      color_limits = c(-3,3),
                      my_bin_size = 10,
                      split,
                      show_names = FALSE
)
{
        
        sense_mat <- sense_mat[rownames(sense_mat) %in% gene_list,]
        sense_mat <- sense_mat[order(rownames(sense_mat)),]
        
        antis_mat <- antis_mat[rownames(antis_mat) %in% gene_list,]
        antis_mat <- antis_mat[order(rownames(antis_mat)),]
        
        anno_df <- anno_df[anno_df$gene_id %in% rownames(antis_mat),]
        anno_df <- anno_df[order(anno_df$gene_id),]
        
        if(!(is.null(split))){
                split <- split[names(split) %in% rownames(antis_mat)]
                split <- split[order(names(split))]
        }
        
        
        my_order <- order(anno_df$width)
        my_width <- anno_df$width[my_order]
        
        sense_mat <- sense_mat[my_order,]
        antis_mat <- antis_mat[my_order,]
        
        if(!(is.null(split))){
                split <- split[my_order]
        }
        
        top_anno_sense <- HeatmapAnnotation(
                points = anno_points(colMeans(sense_mat), border=FALSE, size = unit(1, "mm"), axis=TRUE, gp = gpar(col="darkgrey")))
        
        top_anno_antis <- HeatmapAnnotation(
                points = anno_points(colMeans(antis_mat), border=FALSE, size = unit(1, "mm"), axis=TRUE, gp = gpar(col="darkgrey")))
        
        
        sense_mat[sense_mat > color_limits[2]] <- color_limits[2]
        sense_mat[sense_mat < color_limits[1]] <- color_limits[1]
        
        antis_mat[antis_mat > color_limits[2]] <- color_limits[2]
        antis_mat[antis_mat < color_limits[1]] <- color_limits[1]
        
        my_col_names <- rep("", ncol(sense_mat))
        my_col_names[seq(1000/my_bin_size,length(my_col_names),1000/my_bin_size)] <- c("TSS",1000,2000,3000,4000)[1:(length(my_col_names)/(1000/my_bin_size))]
        
        
        heatmap_anno = HeatmapAnnotation(
                text = anno_text(my_col_names, rot = 0, just = "center", offset = unit(5, "mm"), gp = gpar(cex=1.1)),
                show_legend = FALSE
        )
        
        
        
        my_heatmap_names <- c("log2FC (S)","log2FC (A)")
        
        
        
        ht_diff_sense <- Heatmap(sense_mat,
                                 cluster_columns = FALSE, 
                                 cluster_rows = FALSE,
                                 name = my_heatmap_names[1], 
                                 column_title = "sense",
                                 split = split, gap = unit(2, "mm"),
                                 bottom_annotation = heatmap_anno, bottom_annotation_height = unit(12, "mm"),
                                 top_annotation = top_anno_sense,
                                 top_annotation_height = unit(15, "mm"),
                                 show_row_names = FALSE, show_column_names = FALSE, 
                                 show_heatmap_legend = FALSE,
                                 col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
        
        
        ht_diff_antis <- Heatmap(antis_mat,
                                 cluster_columns = FALSE, 
                                 cluster_rows = FALSE,
                                 name = my_heatmap_names[2], 
                                 column_title = "anti-sense",
                                 split = split, gap = unit(2, "mm"),
                                 heatmap_legend_param = 
                                         list(legend_direction = "horizontal", 
                                              color_bar = "continuous", 
                                              title_position = "lefttop",
                                              title = "log2FC",
                                              title_gp = gpar(fontsize = 14, fontface = "bold"),
                                              legend_width = unit(4,"cm"),
                                              labels_gp = gpar(fontsize = 12)),
                                 bottom_annotation = heatmap_anno, bottom_annotation_height = unit(12, "mm"),
                                 top_annotation = top_anno_antis,
                                 top_annotation_height = unit(15, "mm"),
                                 show_heatmap_legend = TRUE, row_names_gp = gpar(fontsize = (1/nrow(antis_mat)*500)),
                                 show_row_names = show_names, show_column_names = FALSE, 
                                 col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
        
        
        png(paste("log2Diff", gsub(" ","",title),group_name, "png", sep="."),width = 8.27, height = 8.27, units = "in", res = 300)
        
        
        
        draw(ht_diff_sense + ht_diff_antis, 
             heatmap_legend_side = "bottom",
             column_title = title,
             column_title_gp = gpar(fontsize = 18, fontface = "bold"),
             padding = unit(c(1, 2.25, 1, 2.25), "cm"),
             gap = unit(1, "cm"))
        
        
        if(!(is.null(split))){
                my_slices <- unique(split)
        } else {
                my_slices <- 1
        }
        
        for(my_heatmap_name in my_heatmap_names){
                
                for(my_slice in my_slices){
                        
                        
                        decorate_heatmap_body(my_heatmap_name, {
                                
                                i = (1000/my_bin_size)
                                x = i/ncol(sense_mat)
                                grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 3))
                                
                                
                                if(my_slice == max(my_slices)){
                                        for(i in 1:floor(ncol(sense_mat) / 100)){
                                                grid.lines(c(i*x, i*x), unit(c(-5, -2), "mm"), gp = gpar(lwd = 2))
                                                
                                        }
                                }
                                
                                i = (1000/my_bin_size)+(my_width/my_bin_size)
                                x = rev(i/ncol(sense_mat))
                                
                                x[x > 1] <- 1
                                
                                y = (1:nrow(sense_mat))/nrow(sense_mat)
                                
                                grid.circle( x,  y, r = unit(0.3, "mm"), gp = gpar(fill="black") )
                                
                        }, slice = my_slice)
                }
                
                
        }
        
        
        dev.off()
        
        
}




#####################################################################################################################





#####################################################################################################################



genebodyHeatmap <- function(title,
                            group_name,
                            sense_mat, 
                            antis_mat, 
                            gene_list,
                            color_limits = c(-3,3),
                            my_bin_size = 10,
                            order = FALSE
)
{
        
        sense_mat <- sense_mat[rownames(sense_mat) %in% gene_list,]
        antis_mat <- antis_mat[rownames(antis_mat) %in% gene_list,]
        
        # random order
        set.seed(1)
        my_order <- sample(1:nrow(sense_mat), nrow(sense_mat), replace = FALSE)
        
        sense_mat <- sense_mat[my_order,]
        antis_mat <- antis_mat[my_order,]
        
        # ordered my mean
        if(order){
                my_order <- order(rowMeans(antis_mat), decreasing = TRUE)
                
                sense_mat <- sense_mat[my_order,]
                antis_mat <- antis_mat[my_order,]
        }
        
        top_anno_sense <- HeatmapAnnotation(
                points = anno_points(colMeans(sense_mat), border=FALSE, size = unit(1, "mm"), axis=TRUE, gp = gpar(col="darkgrey")))
        
        top_anno_antis <- HeatmapAnnotation(
                points = anno_points(colMeans(antis_mat), border=FALSE, size = unit(1, "mm"), axis=TRUE, gp = gpar(col="darkgrey")))
        
        sense_mat[sense_mat > color_limits[2]] <- color_limits[2]
        sense_mat[sense_mat < color_limits[1]] <- color_limits[1]
        
        antis_mat[antis_mat > color_limits[2]] <- color_limits[2]
        antis_mat[antis_mat < color_limits[1]] <- color_limits[1]
        
        my_axis_marks <- c(1,50,100,200,250,300)
        
        my_col_names <- rep("", ncol(sense_mat))
        my_col_names[my_axis_marks] <- c("-500","TSS","+500","-500","TTS","+500")
        
        heatmap_anno = HeatmapAnnotation(
                text = anno_text(my_col_names, rot = 0, just = "center", offset = unit(5, "mm"), gp = gpar(cex=1.1)),
                show_legend = FALSE
        )
        
        
        my_heatmap_names <- c("log2FC (S)","log2FC (A)")
        
        
        
        ht_diff_sense <- Heatmap(sense_mat,
                                 cluster_columns = FALSE, 
                                 cluster_rows = FALSE,
                                 name = my_heatmap_names[1], 
                                 column_title = "sense",
                                 gap = unit(2, "mm"),
                                 bottom_annotation = heatmap_anno, bottom_annotation_height = unit(12, "mm"),
                                 top_annotation = top_anno_sense,
                                 top_annotation_height = unit(15, "mm"),
                                 show_row_names = FALSE, show_column_names = FALSE, 
                                 show_heatmap_legend = FALSE,
                                 col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
        
        
        ht_diff_antis <- Heatmap(antis_mat,
                                 cluster_columns = FALSE, 
                                 cluster_rows = FALSE,
                                 name = my_heatmap_names[2], 
                                 column_title = "anti-sense",
                                 gap = unit(2, "mm"),
                                 heatmap_legend_param = 
                                         list(legend_direction = "horizontal", 
                                              color_bar = "continuous", 
                                              title_position = "lefttop",
                                              title = "log2FC",
                                              title_gp = gpar(fontsize = 14, fontface = "bold"),
                                              legend_width = unit(4,"cm"),
                                              labels_gp = gpar(fontsize = 12)),
                                 bottom_annotation = heatmap_anno, bottom_annotation_height = unit(12, "mm"),
                                 top_annotation = top_anno_antis,
                                 top_annotation_height = unit(15, "mm"),
                                 show_heatmap_legend = TRUE,
                                 show_row_names = FALSE, show_column_names = FALSE, 
                                 col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
        
        
        png(paste("log2Diff", gsub(" ","",title),group_name, "png", sep="."),width = 8.27, height = 8.27, units = "in", res = 300)
        
        
        draw(ht_diff_sense + ht_diff_antis, 
             heatmap_legend_side = "bottom",
             column_title = title,
             column_title_gp = gpar(fontsize = 18, fontface = "bold"),
             padding = unit(c(1, 2, 1, 2), "cm"),
             gap = unit(2, "cm"))
        
        
        
        for(my_heatmap_name in my_heatmap_names){
                
                
                decorate_heatmap_body(my_heatmap_name, {
                        
                        for(x in  ( (my_axis_marks+1) / ncol(sense_mat))[c(2,5)]){
                                grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 3))
                        }
                        
                        for(x in  ( (my_axis_marks+1) / ncol(sense_mat))){
                                grid.lines(c(x, x), unit(c(-5, -2), "mm"), gp = gpar(lwd = 2))
                        }
                        
                })
        }
        
        
        dev.off()
        
        
}






#####################################################################################################################



#####################################################################################################################
#####################################################################################################################





#####################################################################################################################




zoomHeatmap <- function(title,
                        group_name,
                        sense_mat, 
                        antis_mat, 
                        anno_df,
                        gene_list,
                        color_limits = c(-3,3),
                        my_bin_size = 10,
                        total_log2_cutoff = 0
)
{
        
        sense_mat <- sense_mat[rownames(sense_mat) %in% gene_list,]
        sense_mat <- sense_mat[order(rownames(sense_mat)),]
        
        antis_mat <- antis_mat[rownames(antis_mat) %in% gene_list,]
        antis_mat <- antis_mat[order(rownames(antis_mat)),]
        
        my_sub <- rowMeans(antis_mat) > total_log2_cutoff
        
        sense_mat <- sense_mat[my_sub,]
        antis_mat <- antis_mat[my_sub,]
        
        my_order <- order((rowMeans(antis_mat[,1:(ncol(antis_mat)/2)]) - (rowMeans(antis_mat[,(ncol(antis_mat)/2):ncol(antis_mat)]))))
        
        sense_mat <- sense_mat[my_order,]
        antis_mat <- antis_mat[my_order,]
        
        
        
        
        sense_mat[sense_mat > color_limits[2]] <- color_limits[2]
        sense_mat[sense_mat < color_limits[1]] <- color_limits[1]
        
        antis_mat[antis_mat > color_limits[2]] <- color_limits[2]
        antis_mat[antis_mat < color_limits[1]] <- color_limits[1]
        
        my_axis_marks <- c(  ((ncol(sense_mat)/2)-(500/my_bin_size)+1), (ncol(sense_mat)/2), (ncol(sense_mat)/2)+(500/my_bin_size)  )
        
        my_col_names <- rep("", ncol(sense_mat))
        my_col_names[my_axis_marks] <- c(-500,"TSS",500)
        
        heatmap_anno = HeatmapAnnotation(
                text = anno_text(my_col_names, rot = 0, just = "center", offset = unit(5, "mm"), gp = gpar(cex=1.1)),
                show_legend = FALSE
        )
        
        
        my_heatmap_names <- c("log2FC (S)","log2FC (A)")
        
        
        
        ht_diff_sense <- Heatmap(sense_mat,
                                 cluster_columns = FALSE, 
                                 cluster_rows = FALSE,
                                 name = my_heatmap_names[1], 
                                 column_title = "sense",
                                 gap = unit(2, "mm"),
                                 bottom_annotation = heatmap_anno, bottom_annotation_height = unit(12, "mm"),
                                 show_row_names = FALSE, show_column_names = FALSE, 
                                 show_heatmap_legend = FALSE,
                                 col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
        
        
        ht_diff_antis <- Heatmap(antis_mat,
                                 cluster_columns = FALSE, 
                                 cluster_rows = FALSE,
                                 name = my_heatmap_names[2], 
                                 column_title = "anti-sense",
                                 gap = unit(2, "mm"),
                                 heatmap_legend_param = 
                                         list(legend_direction = "horizontal", 
                                              color_bar = "continuous", 
                                              title_position = "lefttop",
                                              title = "log2FC",
                                              title_gp = gpar(fontsize = 14, fontface = "bold"),
                                              legend_width = unit(4,"cm"),
                                              labels_gp = gpar(fontsize = 12)),
                                 bottom_annotation = heatmap_anno, bottom_annotation_height = unit(12, "mm"),
                                 show_heatmap_legend = TRUE,
                                 show_row_names = FALSE, show_column_names = FALSE, 
                                 col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
        
        
        png(paste("log2Diff", gsub(" ","",title),group_name, "png", sep="."),width = 8.27, height = 8.27, units = "in", res = 300)
        
        
        draw(ht_diff_sense + ht_diff_antis, 
             heatmap_legend_side = "bottom",
             column_title = title,
             column_title_gp = gpar(fontsize = 18, fontface = "bold"),
             padding = unit(c(1, 2.25, 1, 2.25), "cm"),
             gap = unit(1, "cm"))
        
        
        
        for(my_heatmap_name in my_heatmap_names){
                
                
                decorate_heatmap_body(my_heatmap_name, {
                        
                        x =  (ncol(sense_mat)/2) / 100
                        grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 3))
                        
                        for(i in my_axis_marks / 100){
                                grid.lines(c(i, i), unit(c(-5, -2), "mm"), gp = gpar(lwd = 2))
                        }
                        
                })
        }
        
        
        dev.off()
        
        
}



#####################################################################################################################





#####################################################################################################################
#####################################################################################################################









#####################################################################################################################

myRLE <- function (m,  ...) 
{
        # m = log transformed matrix
        
        rle <- m - rowMedians(m)
        return(rle)
}



#####################################################################################################################

myDensLogFC <- function(my_list, title = NULL, xlim = x(-2, 2), ylim = c(0,4), my_colors = my_colors){
        
        
        plot(1, 1, xlim = xlim, ylim = ylim, cex.lab = 1.5, cex.main = 1.5,
             type = 'n', ylab = "Density", xlab="log2 fold change")
        
        for(i in 1:length(my_list)){
                lines(density(my_list[[i]]), col = my_colors[i], lwd=3)
        }
        abline(v = 0, lty=2)
        title(title, line = 0.5, cex.main=1.5)
}

#####################################################################################################################


# similar to Package ‘DAMisc’

# lmobj <- my_lm
# loessobj <- my_loess
# alpha = 0.05

myTestLoess <- function (lmobj, loessobj, alpha = 0.05) 
{
        n <- nrow(model.matrix(lmobj))
        if (n != loessobj$n) {
                stop("Models estimated on different numbers of observations")
        }
        rss0 <- sum(lmobj$residuals^2)
        rss1 <- sum(loessobj$residuals^2)
        d1a <- loessobj$one.delta
        d2a <- loessobj$two.delta
        dfdenom <- d1a^2/d2a
        dfnum <- (n - dfdenom) - lmobj$rank
        F0 <- ((rss0 - rss1)/dfnum)/(rss1/dfdenom)
        cat("F = ", round(F0, 2), "\n", sep = "")
        pval <- pf(F0, dfnum, dfdenom, lower.tail = F)
        cat("Pr( > F) = ", round(pval, 2), "\n", sep = "")
        if (pval < alpha) {
                cat("LOESS preferred to alternative\n")
        }
        else {
                cat("LOESS not statistically better than alternative\n")
        }
        return(pval)
}



#####################################################################################################################






############# common MNase genes #############


commonMNase <- function(my_mats){
        
        for(i in seq_along(my_mats)){
                
                #print(i)
                
                if(i == 1){
                        my_MNase_genes <- rownames(get(my_mats[i]))
                } else {
                        my_MNase_genes <- intersect(my_MNase_genes, rownames(get(my_mats[i])))
                }
        }
        
        for(i in seq_along(my_mats)){
                
                #print(i)
                
                my_mat <- get(my_mats[i])
                
                my_mat <- my_mat[rownames(my_mat) %in% my_MNase_genes,]
                
                my_mat <- my_mat[order(rownames(my_mat)),]
                
                assign(my_mats[i], my_mat, envir = .GlobalEnv)
        }
}


#####################################################################################################################





#####################################################################################################################
#####################################################################################################################








plotComposite <- function(my_sample_mats, 
                          my_sub_range = 1:400,
                          ylims = c(1,5),
                          my_binning = 10,
                          my_colors_composite,
                          my_title = "",
                          yaxt ="s",
                          #site_label = 0,
                          add_line = FALSE,
                          line_lwd = 2,
                          smoother = 11,
                          log_scale = FALSE, 
                          x_axis =TRUE,
                          zscore = FALSE){
        
        
        x_range <- ncol(get(my_sample_mats[1], envir = .GlobalEnv)[,my_sub_range])
        y_data <- colMeans(get(my_sample_mats[1], envir = .GlobalEnv), na.rm = TRUE)[my_sub_range] 
        
        
        plot(1:x_range, y_data, xaxt = "n", yaxt = yaxt,
             main = my_title, xlab = "", ylab = "",
             type="n", ylim = ylims)
        
        #abline(h=1)
        
        
        # axis(side = 1, at = seq(1, x_range , length.out = 3),
        #      labels =  c( paste("-",round((x_range/2)*my_binning/1000), "kb", sep="") ,
        #                   site_label,
        #                   paste("+",round((x_range/2)*my_binning/1000), "kb", sep=""))
        # )
        
        if(x_axis){
                axis(side = 1, at = c(250, 1000, 1750), labels = c("-750", "0", "+750"))
        }
        
        
        
        for(i in seq_along(my_sample_mats)){
                
                my_sample_mat <- get(my_sample_mats[i])[,my_sub_range] 
                
                if(log_scale){
                        my_sample_mat <- log2(my_sample_mat+0.001)
                }
                
                x_range <- ncol(my_sample_mat)
                y_data <- colMeans(my_sample_mat)
                
                if(zscore){y_data <- scale(y_data)[,1]}
                
                my_yerror <- apply(my_sample_mat, 2, function(x){ qt(0.975, df = length(x)-1)*sd(x, na.rm = TRUE)/sqrt(length(x)) })
                
                y_data <- zoo::rollmean(y_data, smoother)
                my_yerror <- zoo::rollmean(my_yerror, smoother)
                
                if(!(zscore)){
                        
                        xx <- c(((smoother-1)/2+1):(x_range-((smoother-1)/2)), (x_range-((smoother-1)/2)):((smoother-1)/2+1))
                        yy <- c(y_data-my_yerror, rev(y_data+my_yerror))
                        
                        polygon(xx, yy, border = NA, col = paste(my_colors_composite[i], "55",sep=""))
                        
                        
                }
                
                
                
                
                lines(((smoother-1)/2+1):(x_range-((smoother-1)/2)), y_data, col = my_colors_composite[i], lwd=line_lwd)
                
                
        }
        
        if(add_line){
                abline(h = c(min(y_data) , max(y_data)), lty=2)
        }
        
}



############################################################################################################################
############################################################################################################################



plotCoverage <- function(my_coverage, 
                         my_region = c("II", 1772954, 1773992),
                         my_bin_size = 10,
                         smoother = NULL,
                         my_color = "#D55E00",
                         my_ylims = c(0,0.6),
                         my_title = "Track",
                         my_peaks = NULL,
                         peak_title = NULL,
                         log_scale = FALSE,
                         y_axis = TRUE,
                         line_lwd = 2,
                         cex_yaxis = 1.25){
        
        par(mar = c(0,3,1,0),mgp=c(2,1,0))
         
        # my_region[2] <- as.numeric(my_region[2]) - 10000
        # my_region[3] <- as.numeric(my_region[3]) + 10000
        # 
        # if(as.numeric(my_region[2]) < 1){my_region[2] <- 1}
        
        my_signal <- as.numeric(my_coverage[[my_region[1]]])[my_region[2]:my_region[3]]
        
        my_cov_subset <- binMeans(y = my_signal, x = seq_along(my_signal), bx = seq(1, length(my_signal), my_bin_size))
        
        if(log_scale){
                my_cov_subset <- log2(my_cov_subset+0.001)
        }
        
        if(!(is.null(smoother))){
                my_cov_subset <- zoo::rollmean(c(rep(my_cov_subset[1],(smoother-1)/2), my_cov_subset, rep(my_cov_subset[length(my_cov_subset)],(smoother-1)/2)), k = smoother)         
        }
        
        
        xrange <- seq(as.numeric(my_region[2])+(my_bin_size/2), as.numeric(my_region[3])-(my_bin_size/2), length.out = length(my_cov_subset))

        plot(xrange,
             rep(0, length(xrange)), 
             type="n", 
             ylim = c(my_ylims[1],my_ylims[2]),
             xaxt = "n", yaxt = "n",ylab = "", main = "", bty = "n", xlab = "")
        
        if(y_axis){
                axis(side = 2, at = c(my_ylims[1],mean(c(my_ylims[1],mean(my_ylims))),mean(my_ylims),mean(c(my_ylims[2],mean(my_ylims))),my_ylims[2]),
                     c(round(my_ylims[1]),"","","",my_ylims[2]), cex.axis = cex_yaxis)
                axis(side = 2, at = mean(my_ylims), round(mean(my_ylims)), col = NA, cex.axis = cex_yaxis)
                #axis(side = 2, at = -3:6, labels = -3:6)
                
        }
        
        
        title(main = my_title, line = 0, col.main = my_color)
        
        
        polygon(c(rev(xrange), xrange), 
                c(rep(my_ylims[1]-0.5, length(xrange)), ((my_cov_subset))), 
                col = paste0(my_color, "11"), border = NA)
        
        lines(xrange, my_cov_subset, col = my_color, lwd = line_lwd)
        
        if(!(is.null(my_peaks))){
                
                my_region_gr <- makeGRangesFromDataFrame(data.frame(chr = my_region[1], start = my_region[2], end = my_region[3]))
                
                my_peaks_subset <- subsetByOverlaps(my_peaks, my_region_gr)
                
                if(length(my_peaks_subset) != 0){
                        
                        my_starts <- start(my_peaks_subset)
                        my_starts[my_starts < as.numeric(my_region[2])] <-  as.numeric(my_region[2])
                        
                        my_ends   <- end(my_peaks_subset)
                        my_ends[my_ends > as.numeric(my_region[3])] <-  as.numeric(my_region[3])
                        
                        rect(xleft = my_starts, xright = my_ends, 
                             ytop = my_ylims[1] + ((my_ylims[2] - my_ylims[1])*0.95), 
                             ybottom = my_ylims[1] + ((my_ylims[2] - my_ylims[1])*0.90), 
                             border = NA, col = "#555555", lwd=1)    
                        
                        if(!(is.null(peak_title))){
                                text(x = my_starts[1], 
                                     y = my_ylims[1] + ((my_ylims[2] - my_ylims[1])*0.975), 
                                     labels = peak_title,adj = 0)
                        }
                        
                        
                }
                
        }
        
        
        
}






######################################################


plotAnnotation <- function(my_genes = my_genes, 
                           my_exons = my_exons,
                           my_region = c("II", 1, 1e5),
                           gene_name_color = "#EEEEEE",
                           rect_color = "#555555",
                           text_shift_vertical = 0,
                           text_shift_horizontal = 0,
                           my_xticks = 2*10^3,
                           gene_name_text = TRUE,
                           my_bin_size = 10,
                           fill_color = "white",
                           cex_xaxis = 1.25,
                           lwd_exon = 0.75 ){
        
        par(mar = c(2,3,0,0),mgp=c(1,0.5,0))
        
        my_region_gr <- makeGRangesFromDataFrame(data.frame(chr = my_region[1], 
                                                            start = my_region[2], 
                                                            end = my_region[3]))
        
        xrange <- seq(as.numeric(my_region[2])+(my_bin_size/2), as.numeric(my_region[3])-(my_bin_size/2), my_bin_size)
        
        
        plot(xrange, rep(0, length(xrange)), ylim = c(-0.5,1), 
             type="n", xaxt = "n", yaxt = "n", ylab = "", xlab = "" , bty = "n")
        
        if(!(is.null(my_xticks))){
                axis(side = 1, at =  seq(xrange[1], xrange[length(xrange)], my_xticks), line = -1, col = NA, col.ticks = "black", cex.axis = cex_xaxis,
                     labels =  round(seq(xrange[1], xrange[length(xrange)], my_xticks)/10^3))  
        }
        
        
        my_genes_subset <- subsetByOverlaps(my_genes, my_region_gr)

        my_starts <- start(my_genes_subset)
        my_starts[my_starts < as.numeric(my_region[2])] <-  as.numeric(my_region[2])
        
        my_ends <- end(my_genes_subset)
        my_ends[my_ends > as.numeric(my_region[3])] <-  as.numeric(my_region[3])
        
        if(length(my_starts) > 0 & sum(strand(my_genes_subset) == "+") > 0){
                rect(xleft = my_starts[as.logical(strand(my_genes_subset) == "+")], 
                     xright = my_ends[as.logical(strand(my_genes_subset) == "+")], 
                     ytop =  0.7875, ybottom = 0.7375, border = rect_color, lwd=3)
        }
                
        if(length(my_starts) > 0 & sum(strand(my_genes_subset) == "-") > 0){
                rect(xleft = my_starts[as.logical(strand(my_genes_subset) == "-")], 
                     xright = my_ends[as.logical(strand(my_genes_subset) == "-")], 
                     ytop =  0.2625, ybottom = 0.2125, border = rect_color, lwd=3)
        }
        
        
        my_exons_subset <- subsetByOverlaps(my_exons, my_region_gr)

        my_exon_starts <- start(my_exons_subset)
        my_exon_starts[my_exon_starts < as.numeric(my_region[2])] <-  as.numeric(my_region[2])
        
        my_exon_ends <- end(my_exons_subset)             
        my_exon_ends[my_exon_ends > as.numeric(my_region[3])] <-  as.numeric(my_region[3])
        
        
        if(length(my_exon_starts) > 0 & sum(strand(my_exons_subset) == "+") > 0){
                rect(xleft = my_exon_starts[as.logical(strand(my_exons_subset) == "+")], 
                     xright = my_exon_ends[as.logical(strand(my_exons_subset) == "+")], 
                     ytop =  1.00, ybottom = 0.525, col = fill_color, border = rect_color, lwd=lwd_exon)
        }
        
        if(length(my_exon_starts) > 0 & sum(strand(my_exons_subset) == "-") > 0){
                rect(xleft = my_exon_starts[as.logical(strand(my_exons_subset) == "-")], 
                     xright = my_exon_ends[as.logical(strand(my_exons_subset) == "-")], 
                     ytop =  0.475, ybottom = 0.00, col = fill_color, border = rect_color, lwd=lwd_exon)
        }
        
        if(gene_name_text){
                text(x = rowMeans(cbind(my_starts[!(is.na(my_genes_subset$gene_name))], my_ends[!(is.na(my_genes_subset$gene_name))])) + text_shift_horizontal, 
                     y = 0.55 + text_shift_vertical, 
                     labels = my_genes_subset$gene_name[!(is.na(my_genes_subset$gene_name))], 
                     col =  gene_name_color, cex = 0.8)
                
        }
        
}




######################################################









############################################################################################################################
############################################################################################################################



binMeans2 <- function(my_signal, my_bin_size){
        my_cov_subset <- data.frame(signal = my_signal[1 : (length(my_signal) -  (length(my_signal) %% my_bin_size)  )],
                                    binner = rep(1:floor((length(my_signal)/my_bin_size)), each = my_bin_size))
        
        my_cov_subset <- aggregate(my_cov_subset$signal, by = list(my_cov_subset$binner), FUN = mean, na.rm = TRUE )
        my_cov_subset <- my_cov_subset$x
        
        my_cov_subset
}




############################################################################################################################
############################################################################################################################



