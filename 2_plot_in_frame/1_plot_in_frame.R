plot_circSpan_1 <- function(datasets, geneName = 'ENSG00000135776.4'){
    ## 2) get data for a geneID ====
    data.geneID <- get_geneID_data(datasets, geneName)
    print(paste('nrow(data.geneID$df.geneID.transcripts) : ',nrow(data.geneID$df.geneID.transcripts)))
    print(paste('nrow(data.geneID$df.geneID.circ) : ',nrow(data.geneID$df.geneID.circ)))
    ## 3) plot circSpan ====
    data.geneID <- plot_circSpan(data.geneID)
    return(data.geneID)
}
# x <- plot_circSpan_1(datasets, 'ENSG00000198691.11')
# x <- plot_circSpan_1(datasets)
# remove(x)


plot_circSpan <- function(data.geneID){
    ## 0.1) functions ====
    seg2 <- function(x,y,hight,width){
        segments(x,y,x,y-hight, col= 'green')
        segments(x,y-hight,x+width,y-hight, col= 'green')
        segments(x+width,y-hight,x+width,y, col= 'green')
        segments(x+width,y,x,y, col= 'green')
    }
    
    ## 0.2) basic values ====
    df.geneID.circ <- data.geneID$df.geneID.circ
    # browser()
    df.geneID.transcripts <- data.geneID$df.geneID.transcripts
    max.value <- data.geneID$max.value
    base <- data.geneID$base
    
    dir.current <- paste(dir.base, '2_plot_in_frame', sep = '/')
    r1 <- nrow(df.geneID.transcripts)
    r2 <- nrow(df.geneID.circ)
    r.all <- r1+r2+1
    print(paste('r.all: ', r.all))
    res <- 300
    p.width <- 5
    p.width.px <- res*p.width
    padding <- 30
    
    df.th <- data.frame(
        ID=character(), 
        LT=character(), 
        stringsAsFactors=FALSE) 
    
    ## 0.3) advanced values ====
    plot.width <- p.width.px - padding*2
    cell <- plot.width/100 
    row.height.in.cell <- 3
    row.height <- row.height.in.cell*cell
    plot.height <-row.height*r.all
    wh.ratio <- plot.width/plot.height
    p.height.px <- plot.height+2*padding
    p.height <- p.height.px/res
    
    
    ## 0.3.1) percentage modifier ====
    perc2px.x <- function(i){
        i <- i*cell + padding
        return(i)
    }
    perc2px.y <- function(i){
        i <- p.height.px - (i*cell + padding)
        return(i)
    }
    co.start <- perc2px.x(10)
    
    ## 2) plot in frame ====
    plot_in_frame <- function(){
        ## comman values ====
        row.accu <- 0
        ## x coordination ===
        co.transID <- perc2px.x(0)
        co.end <- perc2px.x(17)
        co.exon <- perc2px.x(24)
        co.exon.title <- perc2px.x(50)
        co.exp <- perc2px.x(85)
        co.exp.title <- perc2px.x(92)
        
        font.size.title <- 1.5
        font.size.txt <- 1.2
        font.size.txt.circ <- 1
        font.size.exp.value <- 1
        
        
        
        ### the distance from height-line-per-line to exon_rect_position
        hight.exon.box <- 2*cell
        hight.exp.box <- 1.5*cell
        
        
        distance.circ.span <- 0.5*cell
        distance.exp.value.v <- hight.exp.box/2
        distance.exp.value.h <- 0.5*cell
        width.exon <- co.exp - co.exon -1*cell
        width.exp  <- plot.width - co.exp - 1*cell
        
        
        plot(c(0, p.width.px), c(0, p.height.px), ann = F, type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i", yaxs="i")
        # abline(v = perc2px.x(seq(0, 100, by=10)), col = "black", lty =2)
        # abline(h = perc2px.y(seq(0, 100, by=3)), col = "black", lty =2)
        
        
        
        ## plot title ====
        text(x = co.transID, y =perc2px.y(0), "Transcript Name", 
             cex = font.size.title, font = 2, col = "black", adj = c(0,1)
        )
        
        text(x = co.start, y =perc2px.y(0), "Start", 
             cex = font.size.title, font = 2, col = "black", adj = c(0,1)
        )
        text(x = co.end, y =perc2px.y(0), "End", 
             cex = font.size.title, font = 2, col = "black", adj = c(0,1)
        )
        text(x = co.exon.title, y =perc2px.y(0), "Exon structure", 
             cex = font.size.title, font = 2, col = "black", adj = c(0,1)
        )
        # abline(v = co.exp, col = "red")
        # abline(v = co.exon, col = "red")
        
        text(x = co.exp.title, y =perc2px.y(0), "Expression", 
             cex = font.size.title, font = 2, col = "black", adj = c(0,1)
        )
        
        ## 2.1) plot for each transID ====
        for (i in 1:nrow(df.geneID.transcripts)){
            row.accu <- row.accu+1
            transID <- df.geneID.transcripts[i,]$transID
            print(paste('now print for transID:', i))
            row.horizon.position <- perc2px.y(0)-row.accu*row.height
            text(x = co.transID, y = row.horizon.position, transID, 
                 cex = font.size.txt, font = 1, col = "black", adj = c(0,1)
            )
            df.th[nrow(df.th) + 1,] = list(transID, df.geneID.transcripts[i,]$LT)
            
            exon.index <- df.geneID.transcripts[i,]$exon_index[[1]]
            exon.index <- as.numeric(unlist(sapply(exon.index, strsplit, ',', USE.NAMES = FALSE)))
            exon.start <- df.geneID.transcripts[i,]$exons_start[[1]]
            exon.end <- df.geneID.transcripts[i,]$exons_end[[1]]
            
            text(x = co.start, y = row.horizon.position, 
                 paste(exon.start[1], '\nExon: ', exon.index[1], sep = ''),
                 cex = font.size.txt, font = 1, col = "black", adj = c(0,1)
            )
            
            text(x = co.end, y = row.horizon.position, 
                 paste(exon.end[length(exon.end)], '\nExon: ', exon.index[length(exon.index)], sep = ''),
                 cex = font.size.txt, font = 1, col = "black", adj = c(0,1)
            )
            
            exon.start.base <- df.geneID.transcripts[i,]$exons_start_base[[1]]
            exon.end.base <- df.geneID.transcripts[i,]$exons_end_base[[1]]
            ## 2.1.1) plot intron line ====
            transID.border.s <- exon.end.base[1]/max.value*width.exon + co.exon
            transID.border.e <- exon.end.base[length(exon.end.base)]/max.value*width.exon + co.exon
            segments(x0 = transID.border.s,
                     y0 = row.horizon.position - hight.exon.box/2,
                     x1 = transID.border.e,
                     y1 = row.horizon.position - hight.exon.box/2,
                     col= 'black')
            ## 2.1.2) plot exon boxes ====
            for (j in 1:length(exon.start.base)){
                
                co.point.s <- exon.start.base[j]/max.value*width.exon + co.exon
                co.point.e <- exon.end.base[j]/max.value*width.exon + co.exon
                
                rect(xleft   = co.point.s, 
                     xright  = co.point.e, 
                     ytop    = row.horizon.position, 
                     ybottom = row.horizon.position - hight.exon.box, 
                     density = NULL, 
                     col = 'black', border = NULL
                )
            }
            ## 2.1.3) plot circ span/exp ====
            df.geneID.circ.sub <- df.geneID.circ[df.geneID.circ$transID == transID,]
            print(transID)
            print(paste('row.accu: ',row.accu))
            print(paste('circRNA, transID, nrow for sub;', nrow(df.geneID.circ.sub)))
            if (nrow(df.geneID.circ.sub) > 0){
                for(m in 1:nrow(df.geneID.circ.sub)){
                    df <- df.geneID.circ.sub
                    row.accu <- row.accu + 1
                    row.horizon.position <- perc2px.y(0)-row.accu*row.height
                    
                    ## 2.1.31) plot circ span ====
                    text(x = co.transID, y = row.horizon.position, df[m,]$teID, 
                         cex = font.size.txt.circ, font = 1, col = "blue", adj = c(0,1)
                    )
                    df.th[nrow(df.th) + 1,] <- list(as.character(df[m,]$teID), df.geneID.transcripts[i,]$LT)
                    
                    text(x = co.start, y = row.horizon.position, 
                         paste(df[m,]$start, '\nExon: ', df[m,]$exon_index_start, sep = ''),
                         cex = font.size.txt, font = 1, col = "black", adj = c(0,1)
                    )
                    
                    text(x = co.end, y = row.horizon.position, 
                         paste(df[m,]$end, '\nExon: ', df[m,]$exon_index_end, sep = ''),
                         cex = font.size.txt, font = 1, col = "black", adj = c(0,1)
                    )
                    co.circ.s <- (df[m,]$start-base+1)/max.value*width.exon + co.exon
                    co.circ.e <- (df[m,]$end-base+1)/max.value*width.exon + co.exon
                    print(paste('co.circ.s, co.circ.e : ', co.circ.s, ':', co.circ.e, sep = ' '))
                    segments(x0 = co.circ.s,
                             y0 = row.horizon.position - distance.circ.span,
                             x1 = co.circ.e,
                             y1 = row.horizon.position - distance.circ.span,
                             col= 'green',
                             lwd=2
                    )
                    ## 2.1.32) plot exp ====
                    max.circ.exp <- max(df$readCountSum)
                    co.circ.exp.s <- co.exp
                    co.circ.exp.e <- df[m,]$readCountSum/max.circ.exp*width.exp + co.exp
                    rect(xleft   = co.circ.exp.s, 
                         xright  = co.circ.exp.e, 
                         ytop    = row.horizon.position, 
                         ybottom = row.horizon.position - hight.exp.box, 
                         density = NULL, 
                         col = 'grey', border = NULL
                    )
                    ## 2.1.33) plot exp value ====
                    text(x = co.circ.exp.e + distance.exp.value.h, 
                         y = row.horizon.position - distance.exp.value.v, 
                         round(df[m,]$readCountSum,3), 
                         cex = font.size.exp.value, font = 1, col = "black", adj = c(0,1)
                    )
                }
                ## 2.1.4) exp baseline ====
                row.horizon.position.bottom <- perc2px.y(0)-row.accu*row.height
                row.horizon.position.top <- perc2px.y(0)-(row.accu-nrow(df.geneID.circ.sub)+1)*row.height
                segments(x0 = co.exp,
                         y0 = row.horizon.position.top + 0.5*cell,
                         x1 = co.exp,
                         y1 = row.horizon.position.bottom - row.height,
                         col= 'black')
            }
        }
        ## 2.2) exon markline ====
        row.horizon.position.top <- perc2px.y(0)-row.height*2
        row.horizon.position.bottom <- row.height + padding
        exon.start.base <- df.geneID.transcripts[1,]$exons_start_base[[1]]
        exon.end.base <- df.geneID.transcripts[1,]$exons_end_base[[1]]
        for (j in 1:length(exon.start.base)){
            co.point.s <- exon.start.base[j]/max.value*width.exon + co.exon
            co.point.e <- exon.end.base[j]/max.value*width.exon + co.exon
            
            segments(x0 = co.point.s,
                     y0 = row.horizon.position.top,
                     x1 = co.point.s,
                     y1 = row.horizon.position.bottom,
                     col= 'red',
                     lty = 3,
                     lwd = 0.3
            )
            
            
            segments(x0 = co.point.e,
                     y0 = row.horizon.position.top,
                     x1 = co.point.e,
                     y1 = row.horizon.position.bottom,
                     col= 'blue',
                     lty = 3,
                     lwd = 0.3
            )
        }
        
        saveRDS(df.th, file = paste(dir.base, '2_plot_in_frame/df.th.rdata', sep = '/'))
        saveRDS(co.start, file = paste(dir.base, '2_plot_in_frame/co.start.rdata', sep = '/'))
        saveRDS(p.height.px, file = paste(dir.base, '2_plot_in_frame/p.height.px.rdata', sep = '/'))
        print(df.th)
    }
    
    ## 3.1) plot and  save as svg ====
    file.output <- paste(tempfile('a'),'.png', sep = '')
    if(nrow(df.geneID.transcripts) > 0){
        print(paste('width:', p.width))
        print(paste('wh.ratio:', wh.ratio))
        png(file.output,
            width     = p.width,
            height    = p.height,
            units     = "in",
            res       = res,
            pointsize = p.width/2
        )
        par(mar = c(0, 0, 0, 0))
        plot_in_frame()
        graphics.off()
    }
    # # install.packages('rsvg')
    # library(magick)
    # library(rsvg)
    # img <- image_read_svg('./test.svg')
    # plot(img) # or print(img)
    ## return ====
    data.geneID$df.th <- df.th
    data.geneID$file.output <- file.output
    # data.geneID$padding <- padding
    data.geneID$max.value <- NULL
    return(data.geneID)
    
}


# plot_circSpan