library(stringr)
format_ID_lite <- function(i, IDtype){
    i1 <- i[i[[IDtype]] != '',]
    a <- as.character(i1[[IDtype]])
    b <- a[str_detect(a, "EN.*")]
    c <- unlist(sapply(b, strsplit, '\\.', USE.NAMES = FALSE))
    d <- c[seq(1, length(c), by=2)]
    i1[[IDtype]] <- d
    length(d)
    return(i1)
}



Prepare_dataset <- function(){
    ## 1) input database ====
    ## 1.1) circRNA full table, normalized, unconvoluted, average [cft.lite, circGenes]====
    circRNA.final.table <- paste(dir.base, '0_data', 'cirRNA.readCount.normalized.average2.rdata', sep = '/')
    cft <- readRDS(circRNA.final.table)
    cft <- format_ID_lite(cft, 'transID')
    cft <- format_ID_lite(cft, 'geneID')
    cft$readCountMean <- rowMeans(cft[,13:72])
    cft$readCountMedian <- apply(cft[,13:72], 1, median, na.rm = T)
    cft$readCountSum <- rowSums(cft[,13:72])
    cft.lite <- cft[,c(1:12,75)]
    
    circGenes <-unique(as.character(factor(cft$geneID)))
    
    ## 1.2) transcripts [df.te.circGene]  ====
    transcript.exons <- '~/yan150/reference_data/genome_data_UCSC/hg38/gtf/hg38_allTranscriptID_exons/all_transcriptIDs_exons_LT'
    df.transcript.exons <- read.table(transcript.exons, sep = '\t', header = TRUE)
    colnames(df.transcript.exons) <- c('chrID', 'geneID','transID', 'strand', 'exons_start', 'exons_end', 'exon_index', 'LT')
    geneNames <- '~/yan150/reference_data/genome_data_UCSC/hg38/gtf/geneName.txt'
    geneNames <- read.table(geneNames, sep = '\t', header = FALSE)[,c(1,2)]
    colnames(geneNames) <- c('geneID', 'geneName')
    df.transcript.exons.2 <- merge(df.transcript.exons, geneNames, by.x = 'geneID', by.y = 'geneID', all.x = TRUE)
    df.te.circGene <- format_ID_lite(df.transcript.exons.2, 'geneID')
    df.te.circGene <- format_ID_lite(df.te.circGene, 'transID')
    df.te.circGene.2 <- df.te.circGene[df.te.circGene$geneID %in% circGenes,]
    # head(df.te.circGene.2$transID)
    
    

    ## 1.3) circRNA convoluted ====
    
    df.circ.conv <- readRDS('./0_data/df.circRNA.normalized.average.convoluted.rdata')
    df.circ.conv <- format_ID_lite(df.circ.conv, 'geneID')
    x <- rowSums(df.circ.conv[,c(5:ncol(df.circ.conv))])
    df.circ.conv$circRNA.sumReadCount <- round(x,2)
    df.circ.conv.lite <- df.circ.conv[,c(1:4, ncol(df.circ.conv))]
    
    
    ## update geneName ====
    source('./4_go_search/1_go_search.r')
    df.go <- load_df.go()
    df.te.circGene.2 <- update_geneName(df.te.circGene.2, df.go)
    cft <- update_geneName(cft, df.go)
    cft.lite <- update_geneName(cft.lite, df.go)
    df.circ.conv <- update_geneName(df.circ.conv, df.go)
    df.circ.conv.lite <- update_geneName(df.circ.conv.lite, df.go)
    
    ## 121 m.git ====
    m.git <- df.te.circGene.2[,c(1,2,4)]
    m.git2 <- m.git[!duplicated(m.git$transID),]
    
    
    ##  return result ====
    
    list1 <- list(df.te.circGene = df.te.circGene.2,
                  df.te.circGene.git = m.git2,
                  cft = cft, 
                  cft.lite = cft.lite, 
                  df.circ.conv = df.circ.conv,
                  df.circ.conv.lite = df.circ.conv.lite,
                  df.go = df.go)
    saveRDS(list1, file = paste(dir.base, '0_data', 'datasets.rdata', sep = '/'))
}
# Prepare_dataset()

load_dataset<- function(){
    list1 <- readRDS(file = paste(dir.base, '0_data', 'datasets.rdata', sep = '/'))
    x <- list1$df.go
    list1$df.go <- x[x$GO.term.name != '',]
    return(list1)
}

# dir.base <- '~/yan150/report_work_weekly/week2020_16_online_panel/2_shiny/3_rstudio_shiny_project'
dir.base <- '/home/lyuya/shiny.server/circ2go'
# setwd(dir.base)
datasets <- load_dataset()
