library(heatmaply)
library(dendextend)

global.order <- 0
plot.stop.marker <- 'normal'

get_single_plotly <- function(x){
    if (nrow(x) == 0){
        return(NULL)
    }
    x1 <- x[,c(3,13:(ncol(x)-3))]
    rownames(x1) <- x1$teID
    x2 <- x1[,-1]
    x3 <- sapply(x2, as.numeric)
    if(nrow(x2) == 1){
        df.hml$type == 'scp'
        x31 <- cbind(c(1:length(x3)),x3)
        colnames(x31) <- c('cellLine', 'readCount')
        x32 <- data.frame(x31)
        x32$readCount <- round(x32$readCount, 3)
        
        if(plot.stop.marker == 'order'){
            x32 <- x32[match(global.order,rownames(x32)),]
        }
        rownames(x32) <- str_match(rownames(x32), '\\dD_(.*)_\\d')[,2]
        x32$cellLine <- c(1:nrow(x32))
        p1 <- plot_ly(data = x32,
                      x = ~cellLine,
                      y = ~readCount,
                      hoverinfo = "text",
                      hovertext = paste("cell line :", rownames(x32),
                                        "<br> readCount :", x32$readCount))
        # p1 <- p1 %>% layout(xaxis = list(title = paste("cell lines for <br>", x1$teID, sep = '')))
        p1 <- p1 %>%
            layout(
                xaxis = list(tickangle = -45,
                             title = paste("cell lines for <br>", x1$teID, sep = ''),
                             ticktext = as.list(rownames(x32)), 
                             tickvals = as.list(1:nrow(x32)),
                             tickmode = "array"
                ))
        
        plot.stop.marker <<- 'stop'
        return(p1)
    }
    rownames(x3) <- rownames(x1)
    x3 <- round(x3, 3)
    
    # x4 <- log10(x3 +1 )
    # dentree.count <- nrow(x4) -1
    # if (dentree.count >3){
    #     dentree.count <- 3
    # }
    
    
    hc.row <- hclust(dist(x3, method="euclidean"), method="complete")
    hc.col <- hclust(dist(t(x3), method="euclidean"), method="complete")
    oc <- hc.col$order
    or <- hc.row$order
    x3 <- x3[or,oc]
    
    global.order <<- colnames(x3)
    colnames(x3) <- str_match(colnames(x3), '\\dD_(.*)_\\d')[,2]
    plot.stop.marker <<- 'order'
    
    width <- 1200
    height <- nrow(x3)/ncol(x3)*width+300
    
    df.hml$width <<- width
    df.hml$height.hp <<- height
    
    hp1 <- heatmaply(x3, dendrogram = "none",
                     Rowv = NULL,
                     Colv = NULL,
              # k_row = dentree.count, k_col = dentree.count,
              show_dendrogram = c(FALSE, FALSE),
              row_text_angle = 45,
              column_text_angle = 45,
              # showticklabels = c(FALSE,FALSE),
              margins = c(10,10,NA,0)
              )

    return(hp1)
}

heatmap_geneID <- function(geneID){
    source('./1_transcript_coordinate_scale/2_data_for_geneID.r')
    ## get geneID ====
    geneID.a <- format_geneID(geneID)
    df1 <- datasets$cft
    list1 <- list(p1 = NULL, p2 = NULL)
    
    ## plot for geneID ====
    df.circ <- df1[df1$geneID == geneID.a,]
    list1$p1 <- get_single_plotly(df.circ)
    if(plot.stop.marker == 'stop'){
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
