library(heatmaply)
get_single_plotly <- function(x){
    # browser()
    if (nrow(x) == 0){
        return(NULL)
    }
    x1 <- x[,c(3,13:(ncol(x)-3))]
    rownames(x1) <- x1$teID
    x2 <- x1[,-1]
    x3 <- sapply(x2, as.numeric)
    if(nrow(x2) == 1){
        plot.double.marker <- 1
        x31 <- cbind(c(1:length(x3)),x3)
        colnames(x31) <- c('cellLine', 'readCount')
        x32 <- data.frame(x31)
        # browser()
        p1 <- plot_ly(data = x32,
                      x = ~cellLine,
                      y = ~readCount,
                      hoverinfo = "text",
                      hovertext = paste("cell line :", rownames(x32),
                                        "<br> readCount :", x32$readCount))
        p2 <- p1 %>% layout(xaxis = list(title = paste("cell lines for <br>", x1$teID, sep = '')))
        return(p2)
    }
    rownames(x3) <- rownames(x1)
    x4 <- log10(x3 +1 )
    dentree.count <- nrow(x4) -1
    if (dentree.count >3){
        dentree.count <- 3
        }
    hp1 <- heatmaply(x3, k_row = dentree.count, k_col = dentree.count,
                     showticklabels = c(FALSE,FALSE),
                     margins = c(10,10,NA,0))
    # ?heatmaply
    return(hp1)
}

heatmap_geneID <- function(geneID){
    source('./1_transcript_coordinate_scale/2_data_for_geneID.r')
    
    ## get geneID ====
    geneID.a <- format_geneID(geneID)
    df1 <- datasets$cft
    plot.double.marker <- 0
    list1 <- list(p1 = NULL, p2 = NULL)
    
    ## plot for geneID ====
    df.circ <- df1[df1$geneID == geneID.a,]
    list1$p1 <- get_single_plotly(df.circ)
    if(plot.double.marker){
        return(list1)
    }
    
    ## plot for teID if possble ====
    if(str_detect(geneID, 'chr')){
        teID <- geneID
        # teID <- 'chr19_18855925_18856089_+'
        df.circ.v <- df1[df1$teID == teID,]
        list1$p2  <- get_single_plotly(df.circ.v)
        return(list1)
    }
    return(list1)
}

## 
