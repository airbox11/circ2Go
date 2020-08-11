source('./3_shiny_frame/3.1_functions.r')
source('./3_shiny_frame/3.1_ui_modules.r')
    
# 1.3) ui ====
ui <- dashboardPage(
    title = "circ2GO",
    dashboardHeader(title = img(src = "image/c2g.title.png", height = 60)),
    sidebar,
    body
)
## 1.4) server ====
server <- function(input, output, session) {
    ## load data: ====
    ## common vals ====
    source('./5_heatmap/1_heatmap_circ.r')
    vals <- reactiveValues(co.default = c(0,0), 
                           geneName = 'ENSG00000135776.4', 
                           tb.term = data.frame(),
                           tb.info = data.frame(),
                           GoInfoTable.IDselected = '',
                           plotlys = list(p1 = NULL, p2 = NULL),
                           mpV1 = 'ENSG00000135776',
                           boxInfo = list(geneName='',
                                          alias = '',
                                          descrp = '',
                                          geneID = '',
                                          chrID= '',
                                          strand= '',
                                          start='',
                                          end=''),
                           go.term.defi = '',
                           miId = ''
                           )
    ## common function ====
    renderUiInfobox <- function(x){
      uit <- renderUI({
        if(isolate(length(x$geneID)) != 0){
          coordinate <- paste(x$chrID, 
                              x$start,
                              x$end,
                              x$strand,
                              sep = '')
        }else{
          coordinate <- ''
        }
        HTML(paste(sep = '',
                   '<span><strong>Gene name: &nbsp</strong> ', x$geneName, '</span><br>
                      <span><strong>Gene ID: &nbsp</strong> ',x$geneID, '</span><br>
                      <span><strong>Alias: &nbsp</strong> ', x$alias, '</span><br>
                      <span><strong>Description: &nbsp</strong> ', x$descrp, '</span><br>
                      <span><strong>Genomic Location: &nbsp</strong> ',coordinate, '</span><br>'
        ))
      })
      return(uit)
    }
    
    renderUItermtb <- function(x){
      if (isolate(nrow(x))==0){
        uif <- renderUI(fluidRow())
        return(uif)
      }
      output$GoTermTable <- renderDT(x[,c(-3)], selection = 'single', rownames = FALSE)
      
      uif <- renderUI({fluidRow(box(width = 12,
                                    fluidRow(DTOutput("GoTermTable"))
                                    ))})
      return(uif)
    }
    
    renderUIverb <- function(x){
      if (isolate(nchar(x))==0){
        uif <- renderUI(fluidRow())
        return(uif)
      }
      output$goTermDefi <- renderText(x)
      
      uif <- renderUI({
        fluidRow(box(width = 12,
                     tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
                     HTML(paste('<h5><strong>Definition for selected GO term:&nbsp&nbsp&nbsp</strong><i style = "color: #0000a0">',vals$go.term.name,'</i></h5>'),sep= ''),
                     verbatimTextOutput(outputId = 'goTermDefi')
                     ))
        })
      return(uif)
    }
    
    renderUIinfotb <- function(tb.info){
      if (isolate(nrow(tb.info))==0){
        uif <- renderUI(fluidRow())
        return(uif)
      }
      output$GoInfoTable <- renderDT(tb.info, selection = 'single', rownames = FALSE)
      output$DL_goInfo <- downloadHandler(
        filename = function() {
          paste("data_", vals$go.term.name, ".csv", sep="")
        },
        content = function(file) {
          df <- vals$tb.info
          df <- df[order(df$`circRNA sumReadCount`, decreasing = TRUE),]
          write.csv(df, file, row.names = FALSE)
        })
      
      uif <- renderUI(fluidRow(box(width = 12,
                                   HTML(paste('<h4 ><strong>Genes for GO term/accession ID:&nbsp&nbsp&nbsp</strong><i style = "color: #0000a0">',vals$go.term.name,'</i></h4>'),sep= ''),
                                   box(width=12,
                                       fluidRow(DTOutput("GoInfoTable")),
                                       fluidRow(downloadButton('DL_goInfo', 'Download data'))
                                       )
                                   )))
      return(uif)
    }
    
    
    renderUIplotly <- function(x){
      if (isolate(is.null(x$p1))){
        uif <- renderUI(fluidRow())
        return(uif)
      }
      output$plotlyHp1 <- renderPlotly(x$p1)
      
      width <- df.hml$width
      height.hp <- df.hml$height.hp
      height.scp <- 300
      
      if (df.hml$type == 'hp'){
        height2 <- height.hp
      }else{
        height2 <- height.scp
      }
      
      if (isolate(!is.null(x$p2))){
        output$plotlyHp2 <- renderPlotly(x$p2)
        uif <- renderUI(fluidRow(box(width = 12,
          fluidRow(div(plotlyOutput('plotlyHp1', width = width, height = height.hp), align = 'center')),
          fluidRow(div(plotlyOutput('plotlyHp2',width = width, height = height.scp),align = 'center'))
        )))
      }else{
        uif <- renderUI(fluidRow(box(width = 12,
                                     tags$style(HTML("#plotlyHp1 {text-align: center;}")),
                                     div(plotlyOutput('plotlyHp1', width = width, height = height2),align = 'center')))
                        )
      }
      return(uif)
    }
    
    
    ## renderUImi 
    source('./6_miRNA/miTable.r')
    renderUImi <- function(id, df.tm){
      df <- miTable(id, df.tm)
      vals$miDf <- df
      if(isolate(nrow(df) == 0)){
        output$nofound <- renderUI(fluidRow(
          HTML("<div id='div1'>
          <i style='font-size:14px;color: #0000a0;'>
          No circRNA with a predicted microRNA binding site was found.
          </i></div>
               ")
          # edit
          ))
        uif <- renderUI(fluidRow())
      }else{
        output$nofound <- renderUI(fluidRow())
        output$miTb <- renderDT(df, selection = 'single', rownames = FALSE)
        uif <- renderUI(fluidRow(box(width = 12,
                                     DTOutput('miTb'),
                                     downloadButton('miDl', 'Download data')
                                     ))
        )
      }
      return(uif)
    }
    
    renderUiCircGotermtb <- function(id){
      
      if(length(id)==0){
        uif <- renderUI(fluidRow())
        return(uif)
      }else{
        d <- datasets$df.go
        e <- d[d$Gene.stable.ID == id,c(9:12)]
        e <- e[order(e$GO.term.name),]
        e <- e[e$GO.term.name !='',]
        e <- e[!duplicated(e$GO.term.name),]
        colnames(e) <- gsub('\\.', ' ', colnames(e))
        vals$goTable2 <- e
        output$GoTermTableInfobox2 <- renderDT(e[,c(-2)], selection = 'single', rownames = FALSE)
        
        uif <- renderUI(fluidRow(box(width = 12,
                                     # tags$style(HTML("hr {border-top: 1px solid #000000;}")),
                                     # hr(),
                                     fluidRow(DTOutput("GoTermTableInfobox2"))
                                     # br()
        )))
        return(uif)
      }
    }
    
    renderUiCircGotermDef <- function(list1){
      if (length(list1) == 0 ){
        uif <- renderUI(fluidRow())
        return(uif)
      }else{
        go.term.index <- list1[[1]]
        output$goTermDefi2 <- renderText(as.character(vals$goTable2[go.term.index, 2]))
        
        uif <- renderUI(fluidRow(box(width = 12,
                                     tags$h4("Definition for selected GO term", style = "font-weight: bold;"),
                                     # tags$hr(),
                                     verbatimTextOutput(outputId = 'goTermDefi2')
                                     # br()
        )))
        return(uif)
      }}
    
    renderUiCircMap <- function(id){
      # if(id == ''){
      #   uif <- renderUI(fluidRow())
      #   return(uif)
      # }
      data.geneID <- plot_circSpan_1(datasets, id)
      if (nrow(data.geneID$df.geneID.transcripts) == 0){
        uif <- renderUI(fluidRow())
      }else{
        outfile <- data.geneID$file.output
        png1 <- readPNG(outfile)
        width1 <- dim(png1)[2]
        height1 <- dim(png1)[1]
        output$myImage <- renderImage({
          list(src = outfile,
               # contentType = 'image/svg+xml',
               contentType = 'image/png',
               width = width1,
               height = height1,
               alt = "This is alternate text")
        }, deleteFile = TRUE)
        uif <- renderUI(fluidRow(style = paste("height:",(height1+300),"px", sep = ''),
                                 box(width = 12,
                                     style = paste("height:",(height1+200),"px", sep = ''),
                                 tags$style('#div1 {margin-left: 10px;}'),
                                 HTML("<div id='div1'>
                                 <h3><strong>circRNA Transcript Map</strong><h4>
                                 <p>(<span style='color: #000000'><strong>Black</strong></span>: linear transcripts)</p>
                                 <p>(<span style='color: #0000a0'><strong>Blue</strong></span>: circular RNAs)</p>
                                 <i style='font-size:14px;'>Click <span style='color: #0000a0; font-weight: bold;'>circRNA ID</span> on map below to get gene's heatmap and circRNA's scatter plot.</i>
                                 </div>
                                      "),
                                 tags$style(HTML("#myImage {text-align: center;}")),
                                 imageOutput("myImage",
                                             height = '1000px',
                                             click = "image_click")
        ))
        )
        return(uif)
      }
    }
    
    
    ## 140) mainPage ====
    output$logoImage <- renderImage({
        outfile <- './www/image/c2g1.png'
        list(src = outfile,
             contentType = 'image/png'
             # width = 400
             # height = 100
             )
        }, deleteFile = FALSE)

    
    ## 141) circSpan ====
    ## 1411) button: submit ====
    observeEvent(input$ButtonGeneName,{
      withProgress(message = 'In progress',
                   detail = 'This may take a while...',
                   value = 0.7, {
                     
                     ## info box
                     vals$mpV1 <- isolate(input$textInputGeneName)
                     boxInfo <- get_box_info(vals$mpV1)
                     output$uiCircInfobox <- renderUiInfobox(boxInfo)
                     output$uiCircGotermtb <- renderUiCircGotermtb(boxInfo$geneID)
                     ## circMap
                     output$uiCircMap <- renderUiCircMap(vals$mpV1)
                   })
    })
    
    ## 1413) box info table click ====
    observeEvent(input$GoTermTableInfobox2_cell_clicked, {
      list.tb1.select <- isolate(input$GoTermTableInfobox2_cell_clicked)
      output$uiCircGotermDef <- renderUiCircGotermDef(list.tb1.select)
      })
      
    
    ## 1414) image click ====
    observeEvent(input$image_click, {
        x <- isolate(input$image_click$x)
        y <- isolate(input$image_click$y)
        vals$co.default <- c(x,y)
        teID <- get_transID(vals$co.default)
        if (length(teID) == 0){return()}
        showModal(modalDialog(
            easyClose = TRUE,
            tagList(teID), 
            title="Plot heatmap for selected ID?",
            footer = tagList(actionButton("confirmPlot", "Plot"),
                             modalButton("Cancel")
            )
        ))
    })
    observeEvent(input$confirmPlot, {
      withProgress(message = 'In progress',
                   detail = 'This may take a while...',
                   value = 0.7, {
        removeModal()
        ## info box in heatmap page
        boxInfo <- get_box_info(vals$mpV1)
        output$uiInfoboxHp <- renderUiInfobox(boxInfo)
        
        plotlys <- heatmap_geneID(get_transID(vals$co.default))
        output$uiPlotly <- renderUIplotly(plotlys)
        
        updateTabItems(session, "sideBarID", 'TabcircHp')
                   })
    })
    
    
    

    
    ## 142) Go Table ====
    ## 1421) basic go search ====
    format_goterm <- function(x){
      x <- trimws(x, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      
      gon <- unique(df.go$GO.term.accession)
      gon <- gon[!is.na(gon)]
      gon1 <- as.numeric(str_match(gon, 'GO:(\\d+)')[,2])
      df.gon <- data.frame(gon1,gon)
      df.gon <- df.gon[!is.na(df.gon$gon1),]
      
      
      
      gom <- unique(df.go$GO.term.name)
      gom <- gom[!is.na(gom)]
      gomu <- toupper(gom)
      df.gom <- data.frame(gom,gomu)
      
      xup <- toupper(x)
      x1 <- as.numeric(str_match(xup, 'GO:(\\d+)')[2])
      x2 <- as.numeric(str_match(xup, '(\\d+)')[1])
      
      
      df.circ <- datasets$df.circ.conv.lite
      
      if(!is.na(x1)){
        x <- df.gon[df.gon$gon1 == x1,]$gon
        vals$go.term.name <- x
        id <- df.go[df.go$GO.term.accession == x,]$Gene.stable.ID
        id <- unique(id)
        id <- id[!is.na(id)]
        df <- df.circ[df.circ$geneID %in% id,]
        return(df)
      }else if(!is.na(x2)){
        x <- df.gon[df.gon$gon1 == x2,]$gon
        vals$go.term.name <- x
        id <- df.go[df.go$GO.term.accession == x,]$Gene.stable.ID
        id <- unique(id)
        id <- id[!is.na(id)]
        df <- df.circ[df.circ$geneID %in% id,]
        return(df)
      }else{
        x <- df.gom[df.gom$gomu == xup,]$gom
        vals$go.term.name <- x
        id <- df.go[df.go$GO.term.name == x,]$Gene.stable.ID
        id <- unique(id)
        id <- id[!is.na(id)]
        df <- df.circ[df.circ$geneID %in% id,]
        return(df)
      }
    }
    
    observeEvent(input$actBGoTermBasicSubmit,{
      vals$go.term.name <- isolate(input$goSearchBasic)
      vals$tb.info <- format_goterm(vals$go.term.name)
      colnames(vals$tb.info) <- gsub('\\.', ' ', colnames(vals$tb.info))
      output$uiInfotb <- renderUIinfotb(vals$tb.info)
    })
    
    observeEvent(input$actBGoTermBasicClear,{
      updateTextInput(session,inputId = "goSearchBasic", value="")
      vals$tb.info <- data.frame()
      output$uiInfotb <- renderUIinfotb(vals$tb.info)
    })
    
    ## 1422) advanced go search ====
    ## 14221) gotermTb to goInfoTb: buttons/table click ====
    observeEvent(input$actBGoTermSubmit, {
        k1 <- isolate(input$autoC1)
        k2 <- isolate(input$autoC2)
      
        if (k1 == '  ' & k2 == '  '){
          x.s1 <- ''
        }else if(k1 == '  ' | k2 == '  '){
          if(k1 == '  '){
            k3 <- k2
          }else{
            k3 <- k1
          }
          x.s1 <- grep(k3, x.long, value=TRUE)
        }else{
          x.s <- grep(k1, x.long, value=TRUE)
          x.s1 <- grep(k2, x.s, value=TRUE)
        }
        vals$tb.term <- y4[y4$GO.term.name %in% x.s1,]
        colnames(vals$tb.term) <- gsub('\\.', ' ', colnames(vals$tb.term))
        output$uitermtb <- renderUItermtb(vals$tb.term)
        })
    
    observeEvent(input$actBGoTermClear,{
      updateSelectInput(session,inputId = "autoC1", selected = "  ")
      updateSelectInput(session,inputId = "autoC2", selected = "  ")
      vals$tb.term <- data.frame()
      vals$go.term.defi <- ''
      vals$tb.info <- data.frame()
      
      output$uitermtb <- renderUItermtb(vals$tb.term)
      output$uiDefi <- renderUIverb(vals$go.term.defi)
      output$uiInfotb <- renderUIinfotb(vals$tb.info)
    })
  
    observeEvent(input$GoTermTable_cell_clicked, {
        list.tb1.select <- isolate(input$GoTermTable_cell_clicked)
        if (length(list.tb1.select) >0 ){
            go.term.index <- list.tb1.select[[1]]
            vals$go.term.name <- vals$tb.term[go.term.index,2]
            vals$go.term.defi <- vals$tb.term[go.term.index,3]
            output$uiDefi <- renderUIverb(vals$go.term.defi)
            
            x1 <- df.go[df.go$GO.term.name == vals$go.term.name,]
            x11 <- unique(x1$Gene.stable.ID)
            x2 <- datasets$df.circ.conv.lite
            x3 <- x2[x2$geneID %in% x11,]
            vals$tb.info <- x3
            colnames(vals$tb.info) <- gsub('\\.', ' ', colnames(vals$tb.info))
            output$uiInfotb <- renderUIinfotb(vals$tb.info)
            }
        })
    
    
    ##1423) table GO heatmap ====
    
    observeEvent(input$GoInfoTable_cell_clicked,{
        list.tb2.select <- isolate(input$GoInfoTable_cell_clicked)
        if(length(list.tb2.select) > 0){
            go.info.index <- list.tb2.select[[1]]
            vals$GoInfoTable.IDselected <- vals$tb.info[go.info.index,1]
            }
        output$go_vb_geneID_hp <- renderText(vals$GoInfoTable.IDselected)
        })
    
    
    observeEvent(input$GoGeneIDHp,{
      withProgress(message = 'In progress',
                   detail = 'This may take a while...',
                   value = 0.7, {
        boxInfo <- get_box_info(vals$GoInfoTable.IDselected)
        output$uiInfoboxHp <- renderUiInfobox(boxInfo)
        
        plotlys <- heatmap_geneID(vals$GoInfoTable.IDselected)
        output$uiPlotly <- renderUIplotly(plotlys)
        
        updateTabItems(session, "sideBarID", 'TabcircHp')
                   })
        })
    ##1424) table GO circMap ====
    observeEvent(input$GoGeneIDcm,{
      withProgress(message = 'In progress',
                   detail = 'This may take a while...',
                   value = 0.7, {
      boxInfo <- get_box_info(vals$GoInfoTable.IDselected)
      vals$mpV1 <- boxInfo$geneID
      output$uiCircInfobox <- renderUiInfobox(boxInfo)
      output$uiCircGotermtb <- renderUiCircGotermtb(boxInfo$geneID)
      output$uiCircGotermDef <- renderUiCircGotermDef(list())
      output$uiCircMap <- renderUiCircMap(vals$GoInfoTable.IDselected)
      
      updateTabItems(session, "sideBarID", 'circSpan')
                   })
    })
    

    

    
    ## 143) heatmap geneID ====
    ## 1430) infobox/goTable ====
    observeEvent(input$ButtonGeneNameHeatmap, {
      withProgress(message = 'In progress',
                   detail = 'This may take a while...',
                   value = 0.7, {
      namehp <- isolate(input$textInput_geneName_heatmap)
      boxInfo <- get_box_info(namehp)
      output$uiInfoboxHp <- renderUiInfobox(boxInfo)
      
      plotlys <- heatmap_geneID(namehp)
      output$uiPlotly <- renderUIplotly(plotlys)
                   })

        
    })
    
    
    
    ## 1431) downloader ====
    output$DL_circ <- downloadHandler(
      filename = function() { 
        paste("dataset_", isolate(input$textInput_geneName_heatmap), ".csv", sep="")
      },
      content = function(file) {
        df <- readRDS('./0_data/df.geneID.circ.rdata')
        write.csv(df, file)
      })
    
    output$DL_circ_cellLines <- downloadHandler(
      filename = function() { 
        paste("dataset_", isolate(input$textInput_geneName_heatmap), ".csv", sep="")
      },
      content = function(file) {
        df <- readRDS('./0_data/df.geneID.circ.cellLines.rdata')
        write.csv(df, file)
      })
    
    ## 144) circRNA-microRNA ====
    ## 1441) submit/clear button, render table =====
    df.tm <- readRDS('./0_data/df.tr.mr.4.rdata')
    # edit
    
    observeEvent(input$miSearchSubmit,{
      withProgress(message = 'In progress',
                   detail = 'This may take a while...',
                   value = 0.7, {
      vals$miId <- isolate(input$miSearch)
      output$uiMiTb <- renderUImi(vals$miId, df.tm)
                   })
    })
    
    observeEvent(input$miSearchClear,{
      updateTextInput(session,inputId = "miSearch", value="")
      output$uiMiTb <- renderUI(fluidRow())
      output$nofound <- renderUI(fluidRow())
    })
    
    
    ## 1442) table selected id   ====
    observeEvent(input$miTb_cell_clicked,{
      list.tb.select <- isolate(input$miTb_cell_clicked)
      if(length(list.tb.select) > 0){
        index <- list.tb.select[[1]]
        vals$miSelectId <- vals$miDf[index,2]
      }
      output$miSelectId <- renderText(vals$miSelectId)
    })
    
    ## 1443) button Hp   ====
    observeEvent(input$miIdHp,{
      withProgress(message = 'In progress',
                   detail = 'This may take a while...',
                   value = 0.7, {
      boxInfo <- get_box_info(vals$miSelectId)
      output$uiInfoboxHp <- renderUiInfobox(boxInfo)
      
      plotlys <- heatmap_geneID(vals$miSelectId)
      output$uiPlotly <- renderUIplotly(plotlys)
      
      updateTabItems(session, "sideBarID", 'TabcircHp')
                   })
    })
    ## 1444) button Cm   ====
    observeEvent(input$miIdCm,{
      withProgress(message = 'In progress',
                   detail = 'This may take a while...',
                   value = 0.7, {
      boxInfo <- get_box_info(vals$miSelectId)
      vals$mpV1 <- boxInfo$geneID
      output$uiCircInfobox <- renderUiInfobox(boxInfo)
      output$uiCircGotermtb <- renderUiCircGotermtb(boxInfo$geneID)
      output$uiCircGotermDef <- renderUiCircGotermDef(list())
      output$uiCircMap <- renderUiCircMap(vals$miSelectId)
      
      updateTabItems(session, "sideBarID", 'circSpan')
                   })
    })
    
    
    ## 1445) button Download  ====
    output$miDl <- downloadHandler(
      filename = function() { 
        paste('dataset_', vals$miId, '.csv', sep="")
      },
      content = function(file) {
        df <- miTable(vals$miId, df.tm)
        write.csv(df, file, row.names = FALSE)
      })
    
    
    
    ## 145) bulky downlaod ====
    ## 1451) bulky circRNA ====
    output$dlbulk1 <- downloadHandler(
      filename = function() { 
        paste("dataset_circRNA_readCounts.csv", sep="")
      },
      content = function(file) {
        celllines <- isolate(input$cellLineSelect)
        all <- 'all cell lines'
        if(all %in% celllines)celllines<- cellLineNames[-1]
        conv <- isolate(input$radioConv)
        # gotermIn <- isolate(input$GOtermIn)
        if (conv == 'circl'){
          e <- datasets$cft[,c(1:72)]
          colnames(e)[13:72] <- str_match(colnames(e)[13:72], '\\dD_(.*)_\\d')[,2]
          e1 <- e[, celllines]
          e1 <- as.data.frame(e1)
          if (ncol(e1) == 1)colnames(e1) <- celllines
          e2 <- cbind(e[,c(1:12)], e1)
        }else{
          e <- datasets$df.circ.conv[,c(1:64)]
          colnames(e)[5:64] <- str_match(colnames(e)[5:64], '\\dD_(.*)_\\d')[,2]
          e1 <- e[, celllines]
          e1 <- as.data.frame(e1)
          if (ncol(e1) == 1)colnames(e1) <- celllines
          e2 <- cbind(e[,c(1:4)], e1)
        }
        # write.csv(e2[1:100,], file)
        write.csv(e2, file)
      })
    output$dlbulk2 <- downloadHandler(
      filename = function() { 
        paste("dataset_circRNA_readCounts.csv", sep="")
      },
      
      content = function(file) {
        a <- isolate(input$textAi1)
        a <- strsplit(a, split = '\\s|;|,')
        a <- unlist(lapply(a, function(x) x[nchar(x) >= 1]))
        
        conv <- isolate(input$radioConv1)
        if (conv == 'circl'){
          e <- datasets$cft[,c(1:72)]
          colnames(e)[13:72] <- str_match(colnames(e)[13:72], '\\dD_(.*)_\\d')[,2]
          e1 <- e[e$geneName %in% a,]
          e1 <- as.data.frame(e1)
          
        }else{
          e <- datasets$df.circ.conv[,c(1:64)]
          colnames(e)[5:64] <- str_match(colnames(e)[5:64], '\\dD_(.*)_\\d')[,2]
          e1 <- e[e$geneName %in% a,]
          e1 <- as.data.frame(e1)
          
        }

        write.csv(e1, file)
      }
      
    )
    
    ## 1452) bulky miRNA ====
    output$miDlAll <- downloadHandler(
      filename = function() { 
        paste("dataset_circRNA_miRNA_pairs_all.csv", sep="")
      },
      content = function(file) {
        write.csv(df.tm[1:100,], file)
      })
    
    output$miDlIds <- downloadHandler(
      filename = function() { 
        paste("dataset_circRNA_miRNA_pairs_by_Ids.csv", sep="")
      },
      content = function(file) {
        ids <- isolate(input$miMultText)
        df <- miTbs(ids, df.tm)
        write.csv(df, file)
      })
    
} # end of server
