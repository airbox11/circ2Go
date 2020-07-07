
## 0) function: ====
get_box_info <- function(geneID){
    x <- format_geneID(geneID)
    gene.symbol <- readRDS('./0_data/gene_alias.rdata')
    
    
    df <- datasets$df.te.circGene
    df <- df[df$LT == 'LT',]
    geneName <- df[df$geneID == x,]$geneName
    alias <- gene.symbol[gene.symbol$Symbol == geneName,]$Aliases
    descrp <- gene.symbol[gene.symbol$Symbol == geneName,]$description
    chrID <- df[df$geneID == x,]$chrID
    strand <- df[df$geneID == x,]$strand
    start <- split_vector(df[df$geneID == x,]$exons_start)[1]
    end <- split_vector(df[df$geneID == x,]$exons_end)
    end <- end[length(end)]
    list1 <- list(geneName=geneName,
                  alias = alias,
                  descrp = descrp,
                  geneID = x,
                  chrID=chrID,
                  strand=strand,
                  start=start,
                  end=end)
    return(list1)
}

get_transID <- function(co.default){
    ## fetch df.th : ====
    # df.th <- readRDS(file = '~/yan150/report_work_weekly/week2020_16_online_panel/2_shiny/3_rstudio_shiny_project/2_plot_in_frame/df.th.rdata')
    df.th <- readRDS(file = paste(dir.base, '2_plot_in_frame/df.th.rdata', sep = '/'))
    co.start <- readRDS(file = paste(dir.base, '2_plot_in_frame/co.start.rdata', sep = '/'))
    p.height.px <- readRDS(file = paste(dir.base, './2_plot_in_frame/p.height.px.rdata', sep = '/'))
    padding <- 30
    row.px <- (p.height.px - padding*2)/(nrow(df.th)+1)
    x <- co.default[1]
    y <- co.default[2]
    row.index <- (y - padding)%/%row.px
    ID <- df.th[row.index,]$ID
    
    if (co.start > x & x> padding){
        return(ID)
    }else{
        return('empty')
    }
}

## 1) func: get_shiny_modules ====
get_shiny_modules <- function(){
    ## 1.1) common values ====
    source('./4_go_search/1_go_search.r')
    df.go <- datasets$df.go
    x <- df.go$GO.term.name
    x.short <- unique(split_vector(x, sep = ' '))
    x.long <- unique(x)
    
    y <- df.go[,c(8:12)]
    y2 <- y[y$GO.term.name != '',]
    y4 <- y2[!duplicated(y2$GO.term.name),]
    
    ## 1.11) common value : genen name, geneID, transID 
    # df1 <- datasets$df.te.circGene.git
    # geneName <- unique(as.character(factor(df1$geneName)))
    # geneID <- unique(as.character(factor(df1$geneID)))
    # transID <- unique(as.character(factor(df1$transID)))
    # ggt <- c(geneName)
    # length(ggt)
    
    ## 1.2) modules for panel ====
    ## 1.2.1) sidebar ====
    sidebar <- dashboardSidebar(
        sidebarSearchForm(textId = "searchText", buttonId = "searchButton",
                          label = "Search..."),
        sidebarMenu(id = "sideBarID",
            menuItem(text = "circ2GO main page", tabName = "mainPage", selected = TRUE),
            menuItem(text = "circRNA search", icon = icon("circle"),startExpanded = TRUE,
                     menuSubItem(text = 'circRNA Transcript Map',
                                 tabName = 'circSpan'),
                     menuSubItem(text = 'circRNA Heatmap',
                                 tabName = 'TabcircHp')
                     ),
            menuItem(text = "GO Search", tabName = 'goSearch', icon = icon("circle"), 
                     badgeLabel = "new", badgeColor = "green"),
            menuItem(text = "link", icon = icon("file-code-o"), 
                     href = "https://github.com/rstudio/shinydashboard/"),
            menuItem(text = 'About', tabName = 'tabAbout')
        )
    )
    ## 1.2.2) body ====
    ## 1.2.21) flowR for GO panel ====
    flowR1.GoSelectButton <- fluidRow(
        column(3,
               selectizeInput(inputId = 'autoC1',
                              label = 'Search',
                              choices = x.short,
                              selected = 'shouldNotInTermsDict',
                              size = 10,
                              multiple = FALSE, # allow for multiple inputs
                              options = list(create = FALSE) # if TRUE, allows newly created inputs
                              )
        ),
        
        
        column(3,
               selectizeInput(inputId = 'autoC2',
                              label = 'Search',
                              choices = x.short,
                              selected = 'shouldNotInTermsDict',
                              size = 10,
                              multiple = FALSE, # allow for multiple inputs
                              options = list(create = FALSE)
               ),
        ),
        column(3,
               actionButton(inputId = 'actB.goTerm', label = 'submit' )
        )
    )
    flowR2.table.goTerm <- fluidRow(
        dataTableOutput("GoTermTable")
    )
    flowR2.verbtext.goTerm.defi <- fluidRow(
        tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
        tags$h4("Selected GO term definition"),
        verbatimTextOutput(outputId = 'goTermDefi'),
        tags$hr(style="border-color: black;")
    )
    flowR3.0.goHP <- fluidRow(
       actionButton(inputId = 'GoGeneIDHp', label = 'heatmap geneID')
    )
    
    flowR301.goHP <- fluidRow(
        tags$h4("Selected gene for heatmap:"),
        verbatimTextOutput(outputId = 'go.vb.geneID.hp')
    )
    
    flowR3.table.goInfo <- fluidRow(
        dataTableOutput("GoInfoTable")
    )
    ## 1.2.211) flowR for heatmap panel ====
    flowR.hp.vbtext <- fluidRow(
        tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
        # tags$h4("Selected gene:"),
        # verbatimTextOutput(outputId = 'hp.vb.geneID'
        #                    # placeholder = 'Go term defination'
        # ),
        tags$hr(style="border-color: black;")
    )
    
    flowR.hp.box <- fluidRow(
        tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
        # tags$h4("Selected gene:"),
        box(title = 'Gene Information',
            width = 12,
            uiOutput('ui2')
        )


    )
    
    ## 1.2.22) body itself ====
    body <- dashboardBody(
        tabItems(
            ## tab main page ====
            tabItem(tabName = "mainPage",
                    h2("first tab should be shown.."),
                    imageOutput('logoImage')
            ),
            tabItem(tabName = 'tabAbout',
                    h2("first tab should be shown."),
                    box( width = NULL, status = "warning", solidHeader = TRUE, title="About RBP2GO",
                         fluidRow(column(12,h4(HTML("<strong>How to cite</strong>")))),
                         fluidRow(column(12, h4(HTML("Manuscript in preparation"), style = "color: #000000"))),
                         fluidRow(column(12, h4(HTML("RBP2GO: A comprehensive, pan-species database on RNA-binding proteins, their interactions and functions. Maïwen Caudron-Herger, Ralf E. Jansen, Elsa Wassmer and Sven Diederichs"), style = "color: #000000;"))),  
                         br(),
                         fluidRow(column(12,h4(HTML("<strong>Related publications</strong>")))),
                         fluidRow(column(12, h4("Caudron-Herger et al., R-DeeP: proteome-wide and quantitative identification of RNA-dependent proteins by density gradient ultracentrifugation"))),
                         fluidRow(column(12, h4("2019, Molecular Cell 75, 1–16"))),
                         a(href="https://www.cell.com/molecular-cell/fulltext/S1097-2765(19)30310-7", "doi.org/10.1016/j.molcel.2019.04.018", style="font-size: 17px; color: #2950a3", target="_blank"),
                         br(),
                         br(),
                         br(),
                         fluidRow(column(12, h4("Caudron-Herger et al., Identification, quantification and bioinformatic analysis of RNA-dependent proteins by RNase treatment and density gradient ultracentrifugation using R-DeeP"))),
                         fluidRow(column(12, h4("2020, Nature Protocols 15, 1338–1370"))),
                         a(href="https://www.nature.com/articles/s41596-019-0261-4", "doi.org/10.1038/s41596-019-0261-4", style="font-size: 17px; color: #2950a3", target="_blank"),
                         br(),
                         br(),
                         br(),
                         fluidRow(column(12,h4(HTML("<strong>Related database</strong>")))),
                         fluidRow(column(12, h4(HTML("<a href='http://R-DeeP.dkfz.de' target ='_blank'>R-DeeP: Database for RNA-dependent Proteins</a>")))),
                         br(),
                         br(),
                         fluidRow(column(12,h4(HTML("<strong>Contact</strong>")))),
                         fluidRow(column(12, h4("The authors of the RBP2GO database can be contacted at the following email address") )),
                         a(href="mailto:databases.diederichslab@dkfz.de","databases.diederichslab(at)dkfz.de",style="font-size: 17px; color: #2950a3"),
                         br(),
                         br()
                    ) # box
                    
                    ),
            ## tab circ span ====
            tabItem(tabName = "circSpan",
                    h2("Dashboard tab content"),
                    tags$p("p creates a paragraph of text."),
                    # fluidRow(tags$head(includeScript("./www/returnClick_mainPage.js")),
                    #         textInput("mpTextInput1", "", placeholder = "Enter text then hit return@", width = "100%"),
                    #         actionButton("mpActB1", "Submit"),
                    #         verbatimTextOutput("mpOutput1")
                    # ),
                    
                    
                    fluidRow(
                        tags$head(includeScript("./www/returnClick_circSpan.js")),
                        textInput(inputId = "textInput.geneName",
                                  width = '250px',
                                  # width = "100%",
                                  label = "gene name/gene ID/transcript ID:",
                                  # value = "ENSG00000135776",
                                  placeholder = "eg: ENSG00000135776, ABCB10 "
                        ),
                        # selectizeInput(inputId = 'textInput.geneName',
                        #                label = "gene name/gene ID/transcript ID:",
                        #                choices = ggt,
                        #                selected = 'shouldNotInTermsDict',
                        #                size = 10,
                        #                multiple = FALSE, # allow for multiple inputs
                        #                options = list(create = FALSE) # if TRUE, allows newly created inputs
                        # ),
                        
                        actionButton(inputId = 'Button.geneName',
                                     label = 'Submit'
                                     ),
                        tags$h3(),
                        tags$hr(style="border-color: black;")
                    ),
                    fluidRow(
                      box(title = 'Gene Information',
                          width = 12,
                          uiOutput('ui1')
                          )
                      ),
                    tags$h3(),
                    fluidRow(
                        # actionButton(inputId = 'GoGeneIDHp.circSpan', label = 'circRNA Heatmap'),
                        # helpText('gene heatmap based on selected transcript ID/circRNA ID'),
                        tags$h4("circRNA Transcript Map")
                        # verbatimTextOutput(outputId = 'ID_text')
                        # textOutput("co_xy"),
                        # textOutput("ID_text")
                    ),
                    box(width = 12,
                        imageOutput("myImage",
                                    height = '1000px',
                                    click = "image_click",
                                    hover = hoverOpts(
                                    id = "image_hover",
                                    delay = 0,
                                    delayType = "throttle"),)
                    ),
            ), # tabItem circSpan
            ## tab circ heatmap ====
            tabItem(tabName = 'TabcircHp',
                    h2("Widgets tab content"),
                    fluidRow(
                        tags$head(includeScript("./www/returnClick_circHp.js")),
                        # textInput(inputId = "textInput.geneName.heatmap", 
                        #           label = "gene name/gene ID/transcript ID",
                        #           width = '250px',
                        #           value = "ENSG00000135776.4"
                        # ),
                        textInput(inputId = "textInput.geneName.heatmap", 
                                  width = '250px',
                                  # width = "100%",
                                  label = "gene name/gene ID/transcript ID:",
                                  # value = "ENSG00000135776",
                                  placeholder = "eg: ENSG00000135776, ABCB10 "
                        ),
                        actionButton(inputId = 'Button.geneName.heatmap',
                                     label = 'Submit'
                        )
                    ),
                    flowR.hp.vbtext,
                    flowR.hp.box,
                    fluidRow(
                      downloadButton('DL_circ', 'Download data circ sum'),
                      downloadButton('DL_circ_cellLines', 'Download data circ cellLines')
                    ),
                    fluidRow(plotlyOutput('plotlyHp1')),
                    fluidRow(plotlyOutput('plotlyHp2'))
                    
                    
            ),
            ## tab go search ====
            tabItem(tabName = "goSearch",
                    h2("Widgets tab content"),
                    flowR1.GoSelectButton,
                    flowR2.table.goTerm,
                    flowR2.verbtext.goTerm.defi,
                    flowR3.table.goInfo,
                    flowR3.0.goHP,
                    flowR301.goHP,
                    downloadButton('DL_goInfo', 'Download data')
            )
        )
    )
    
    
    # 1.3) ui ====
    ui <- dashboardPage(
        # dashboardHeader(title = "circ2GO"),
        dashboardHeader(title = HTML('<p style="font-size:100%;"><b style="font-size:160%;">circ2GO</b></p>')),
        sidebar,
        body
    )
    ## 1.4) server ====
    server <- function(input, output, session) {
        ## common vals ====
        source('./5_heatmap/1_heatmap_circ.r')
        vals <- reactiveValues(co.default = c(0,0), 
                               geneName = 'ENSG00000135776.4', 
                               tb.term = data.frame(),
                               tb.info = data.frame(),
                               GoInfoTable.IDselected = 'Nothing select yet',
                               plotlys = list(p1 = NULL, p2 = NULL),
                               mpV1 = 'ENSG00000135776',
                               boxInfo = list(geneName='',
                                              alias = '',
                                              descrp = '',
                                              geneID = '',
                                              chrID= '',
                                              strand= '',
                                              start='',
                                              end='')
        )
        
        
        ## 140) mainPage ====
        output$logoImage <- renderImage({
            outfile <- './www/Selection_048.png'
            list(src = outfile,
                 contentType = 'image/png',
                 width = 400,
                 height = 300,
                 alt = "This is alternate text")
            }, deleteFile = FALSE)
        
        
        ## 141) circSpan ====

        ## 1412) get new geneName ====
        observeEvent(input$Button.geneName,{
            vals$mpV1 <- isolate(input$textInput.geneName)
            vals$boxInfo <- get_box_info(vals$mpV1)
            output$textOutput.geneName <- renderText(vals$mpV1)
            output
            print('what happened?')
        })
            
            
        observeEvent(input$Button.geneName, {
            vals$geneName <- isolate(input$textInput.geneName)
            data.geneID <- plot_circSpan_1(datasets, vals$geneName)
            if (nrow(data.geneID$df.geneID.transcripts) == 0){
                return(NULL)
            }
            # outfile <- paste(dir.base, '2_plot_in_frame/test.png', sep = '/')
            outfile <- data.geneID$file.output
            png1 <- readPNG(outfile)
            width1 <- dim(png1)[2]
            height1 <- dim(png1)[1]
            output$myImage <- renderImage({
                print(paste('height: ',height1))
                print(width1)
                # Return a list containing the filename
                list(src = outfile,
                     # contentType = 'image/svg+xml',
                     contentType = 'image/png',
                     width = width1,
                     height = height1,
                     alt = "This is alternate text")
            }, deleteFile = TRUE)
        })
            
        ##14121) box info ====
        renderUIx <- function(){
            ui1 <- renderUI({
                if(length(vals$boxInfo$geneID) != 0){
                    coordinate <- paste(vals$boxInfo$chrID, ':',
                                        vals$boxInfo$start,'-',
                                        vals$boxInfo$end, ',(',
                                        vals$boxInfo$strand, ')Strand(GRCh38)',
                                        sep = '')
                }else{
                    coordinate <- ''
                }
                HTML(paste(sep = '',
                           '<span><strong>Gene name: &nbsp</strong> ', vals$boxInfo$geneName, '</span><br>
                          <span><strong>Gene ID: &nbsp</strong> ',vals$boxInfo$geneID, '</span><br>
                          <span><strong>Alias: &nbsp</strong> ', vals$boxInfo$alias, '</span><br>
                          <span><strong>Description: &nbsp</strong> ', vals$boxInfo$descrp, '</span><br>
                          <span><strong>Coordiante: &nbsp</strong> ',coordinate, '</span><br>'
                ))
            })
            return(ui1)
        }
        output$ui1 <- renderUIx()
        output$ui2 <- renderUIx()
        
        ## 1413) output click coordination ====
        observeEvent(input$image_click, {
            x <- isolate(input$image_click$x)
            y <- isolate(input$image_click$y)
            vals$co.default <- c(x,y)
            teID <- get_transID(vals$co.default)
            # browser()
            if (teID == 'empty'){return()}
            showModal(modalDialog(
                easyClose = TRUE,
                tagList(
                    # textInput("newfilename", label = "Filename", placeholder = "my_file.txt")
                    # textOutput(outputId = 'ID_text')
                    teID
                ), 
                title="Plot selected circRNA?",
                footer = tagList(actionButton("confirmPlot", "Plot"),
                                 modalButton("Cancel")
                )
            ))
            # browser()

        })
        observeEvent(input$confirmPlot, {
            print('popup confirm')
            removeModal()
            vals$plotlys <- heatmap_geneID(get_transID(vals$co.default))
            updateTabItems(session, "sideBarID", 'TabcircHp')
            # output$hp.vb.geneID <- renderText(get_transID(vals$co.default))
        })
        
        output$co_xy <- renderText({
            vals$co.default
        })
        output$ID_text <- renderText(get_transID(vals$co.default))
        
        
        ## 1414) plot circSpan ====
        output$image_clickinfo <- renderPrint({
            cat("Click:\n")
            str(input$image_click)
        })
        
        ## 142) Go Table ====
        observeEvent(input$actB.goTerm, {
            x.s <- grep(isolate(input$autoC1), x.long, value=TRUE)
            x.s1 <- grep(isolate(input$autoC2), x.s, value=TRUE)
            # browser()
            vals$tb.term <- y4[y4$GO.term.name %in% x.s1,]
            })
        
        output$GoTermTable <- renderDataTable({vals$tb.term[,c(-3)]}, selection = 'single', rownames = FALSE)
      
        observeEvent(input$GoTermTable_cell_clicked, {
            list.tb1.select <- isolate(input$GoTermTable_cell_clicked)
            if (length(list.tb1.select) >0 ){
                go.term.index <- list.tb1.select[[1]]
                go.term.name <- vals$tb.term[go.term.index,2]
                go.term.defi <- vals$tb.term[go.term.index,3]
                x1 <- df.go[df.go$GO.term.name == go.term.name,]
                x11 <- unique(x1$Gene.stable.ID)
                x2 <- datasets$df.circ.conv.lite
                x3 <- x2[x2$geneID %in% x11,]
                vals$go.term.defi <- go.term.defi
                vals$tb.info <- x3
                }
            })
        output$goTermDefi <- renderText(vals$go.term.defi)
        output$GoInfoTable <- renderDataTable({vals$tb.info}, selection = 'single', rownames = FALSE)
        
        ##1423) flow3: table GO heatmap ====
        
        observeEvent(input$GoInfoTable_cell_clicked,{
            list.tb2.select <- isolate(input$GoInfoTable_cell_clicked)
            if(length(list.tb2.select) > 0){
                go.info.index <- list.tb2.select[[1]]
                vals$GoInfoTable.IDselected <- vals$tb.info[go.info.index,2]
                }
            })
        
        output$go.vb.geneID.hp <- renderText(vals$GoInfoTable.IDselected)
        
        observeEvent(input$GoGeneIDHp,{
            print('GoGeneIDHpButton:')
            # browser()
            vals$boxInfo <- get_box_info(vals$GoInfoTable.IDselected)
            vals$plotlys <- heatmap_geneID(vals$GoInfoTable.IDselected)
            # output$hp.vb.geneID <- renderText(vals$GoInfoTable.IDselected)
            updateTabItems(session, "sideBarID", 'TabcircHp')
            })
        
        # observeEvent(input$GoGeneIDHp.circSpan,{
        #     print('GoGeneIDHpButton.circSpan:')
        #     # browser()
        #     vals$plotlys <- heatmap_geneID(get_transID(vals$co.default))
        #     updateTabItems(session, "sideBarID", 'TabcircHp')
        #     # output$hp.vb.geneID <- renderText(get_transID(vals$co.default))
        #     })
        
        ## 1424) download table ====
        output$DL_goInfo <- downloadHandler(
          filename = function() { 
            paste("dataset-", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            write.csv(vals$tb.info, file)
          })
        
        ## 143) heatmap geneID ====
        observeEvent(input$Button.geneName.heatmap, {
            # browser()
            vals$mpV1 <- isolate(input$textInput.geneName.heatmap)
            vals$boxInfo <- get_box_info(vals$mpV1)
            
            vals$plotlys <- heatmap_geneID(vals$mpV1)
            # output$hp.vb.geneID <- renderText(vals$mpV1)
        })
        output$plotlyHp1 <- renderPlotly(vals$plotlys$p1)
        output$plotlyHp2 <- renderPlotly(vals$plotlys$p2)
        
        output$DL_circ <- downloadHandler(
          filename = function() { 
            paste("dataset_", isolate(input$textInput.geneName.heatmap), ".csv", sep="")
          },
          content = function(file) {
            df <- readRDS('./0_data/df.geneID.circ.rdata')
            write.csv(df, file)
          })
        
        output$DL_circ_cellLines <- downloadHandler(
          filename = function() { 
            paste("dataset_", isolate(input$textInput.geneName.heatmap), ".csv", sep="")
          },
          content = function(file) {
            df <- readRDS('./0_data/df.geneID.circ.cellLines.rdata')
            write.csv(df, file)
          })
        
    } # end of server
    ## return ui/server ====
    list1 <- list(ui = ui, 
                  server = server)
}