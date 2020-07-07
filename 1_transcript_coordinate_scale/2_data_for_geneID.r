# ## 0) functions ====
# l.x <- function(x,p){
#     print(length(unique(as.character(factor(x)))))
#     if(p == 1){
#         print(unique(as.character(factor(x))))
#     }else if(p == 2){
#         print(head(unique(as.character(factor(x)))))
#     }
# }
# a <- df.te.circGene$geneID
# b <- unlist(sapply(a, strsplit, '\\.', USE.NAMES = FALSE))
# c <- b[seq(1, length(b), by=2 )]
# length(c)
# l.x(c,2)
# l.x(df.te.circGene$transID,2)
# dim(df.te.circGene)

##0) functions ====
split_vector <- function(x, sep = ',', type = 'char'){
    # x <- df.geneID.transcripts$exons_start[1]
    
    if (type == 'char'){
        x <- unlist(sapply(x, strsplit, sep, USE.NAMES = FALSE))
    }else if (type == 'num'){
        x <- as.numeric(unlist(sapply(x, strsplit, sep, USE.NAMES = FALSE)))
    }
    return(x)
}


format_geneID <- function(geneID){
    m.git <- datasets$df.te.circGene.git
    if(str_detect(geneID, 'ENSG')){
        geneID <- unlist(sapply(geneID, strsplit, '\\.', USE.NAMES = FALSE))[1]
    }else if(str_detect(geneID, 'ENST')){
        transID <- unlist(sapply(geneID, strsplit, '\\.', USE.NAMES = FALSE))[1]
        geneID <- m.git[m.git$transID == transID,]$geneID
    }else if(str_detect(geneID, 'chr')){
        m <- datasets$cft.lite
        geneID <- m[m$teID == geneID,]$geneID
    }else{
        if(geneID == 'TBCE'){
            geneID <- 'ENSG00000285053'
            type <- 'geneID'
        }else if (geneID == 'PINX1'){
            geneID <- 'ENSG00000258724'
            type <- 'geneID'
        }else{
            m <- m.git[,c(1,2)]
            m1 <- m[!duplicated(m[c(1,2)]),]
            geneID <- m1[m1$geneName == geneID, ]$geneID
        }
    }
    cft <- datasets$cft
    df.geneID.circ.cellLines <- cft[cft$geneID == geneID,]
    cft.lite <- datasets$cft.lite
    df.geneID.circ <- cft.lite[cft.lite$geneID == geneID,]
    saveRDS(df.geneID.circ.cellLines, file = './0_data/df.geneID.circ.cellLines.rdata')
    saveRDS(df.geneID.circ, file = './0_data/df.geneID.circ.rdata')
    return(geneID)
}

get_geneID_data <- function(datasets, geneID){
    ##2) get circ and transID for each geneID ====
    # geneID <- 'ENSG00000135776.4'
    geneID.a <- format_geneID(geneID)
    cft.lite <- datasets$cft.lite
    df.te.circGene <- datasets$df.te.circGene
    ### 2.1) circ: ====
    
    df.geneID.circ <- cft.lite[cft.lite$geneID == geneID.a,]
    df.geneID.circ <- df.geneID.circ[order(df.geneID.circ$start),]


    ## 1.2) transIDs ====
    
    # df.geneID.transcripts <- df.te.circGene[df.te.circGene$geneID == geneID,]
    df.geneID.transcripts <- df.te.circGene[df.te.circGene$geneID == geneID.a,]
    df.geneID.transcripts$exons_start <- sapply(df.geneID.transcripts$exons_start, split_vector, ',','num',USE.NAMES = FALSE)
    df.geneID.transcripts$exons_end <- sapply(df.geneID.transcripts$exons_end, split_vector, sep = ',', type = 'num', USE.NAMES = FALSE)
    df.geneID.transcripts <- df.geneID.transcripts[order(df.geneID.transcripts$LT),]
    ## 1.2.1) 
    x <- df.geneID.transcripts$exons_start
    base <- min(unlist(x))
    minxf <- function(x){
        base <<- base
        x <- x-base+1
    }
    
    df.geneID.transcripts$exons_start_base <- sapply(df.geneID.transcripts$exons_start, minxf,USE.NAMES = FALSE)
    df.geneID.transcripts$exons_end_base <- sapply(df.geneID.transcripts$exons_end, minxf,USE.NAMES = FALSE)
    max.value <- max(unlist(df.geneID.transcripts$exons_end)) - base + 1
    
    ## return data:
    list1 <- list(df.geneID.circ = df.geneID.circ, 
                  df.geneID.transcripts = df.geneID.transcripts, 
                  max.value = max.value,
                  base = base)
    return(list1)
}
