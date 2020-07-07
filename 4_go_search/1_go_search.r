get_df.go <- function(){
    df.go <- './4_go_search/mart_export.csv'
    saveRDS(df.go, file = './4_go_search/df.go.rdata')
}

load_df.go <- function(){
    df.go <- './4_go_search/mart_export.csv'
    df.go <- read.csv(file = df.go, sep = '\t', header = TRUE)
    return(df.go)
}

update_geneName <- function(cft, df.go){
    df <- df.go[,c(1,7)]
    df.1 <- df[!duplicated(df[c(1,2)]),]
    cft2 <- merge(cft, df.1, by.x = 'geneID', by.y = 'Gene.stable.ID', all.x = TRUE)
    cft3 <- cft2
    cft3[is.na(cft3$Gene.name),]$Gene.name <- cft3[is.na(cft3$Gene.name),]$geneName
    cft4 <- cft3[, !(colnames(cft3) %in% c("geneName"))]
    colnames(cft4)[ncol(cft4)] <- 'geneName'
    cft4 <- cft4[,c(ncol(cft4),1:(ncol(cft4)-1))]
    return(cft4)
}

