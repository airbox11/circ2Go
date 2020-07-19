
## 0) function: ====
wait <- function(){
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Please wait", value = 0)
  progress$inc(0.95)
}


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
    
    ## database download 
    d <- datasets$cft
    cellLineNames <- str_match(colnames(d)[13:72], '\\dD_(.*)_\\d')[,2]
    cellLineNames <- c('all cell lines', cellLineNames)
    ## 1.2) modules for panel ====
    ## 1.2.1) sidebar ====
    sidebar <- dashboardSidebar(
        sidebarMenu(id = "sideBarID",
            menuItem(text = "circ2GO main page", tabName = "mainPage", selected = TRUE),
            menuItem(text = "circRNA search", icon = icon("bars"),startExpanded = TRUE,
                     menuSubItem(text = 'circRNA Transcript Map',
                                 tabName = 'circSpan'),
                     menuSubItem(text = 'circRNA Heatmap',
                                 tabName = 'TabcircHp')
                     ),
            menuItem(text = "GO Search", tabName = 'goSearch', icon = icon("bars")),
            menuItem(text = 'Data Download', tabName = 'tabDownload', icon = icon('bars')),
            menuItem(text = 'Help', tabName = 'tabHelp'),
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
        )
        ) # flowR1.GoSelectButton

    
    flowR2.table.goTerm <- fluidRow(DTOutput("GoTermTable"))
    flowR2.verbtext.goTerm.defi <- fluidRow(
        tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
        tags$h4("Definition for selected GO term", style = "font-weight: bold;"),
        verbatimTextOutput(outputId = 'goTermDefi'),
        tags$hr(style="border-color: black;")
    )
    flowR3.0.goHP <- fluidRow(
       actionButton(inputId = 'GoGeneIDHp', label = 'Heatmap')
    )
    
    flowR301.goHP <- fluidRow(
        tags$h4("Heatmap for selected gene:", style = "font-weight: bold;padding-left:10px"),
        tags$style("#go_vb_geneID_hp {
    			margin-left:10px;
    			}"),
        verbatimTextOutput(outputId = 'go_vb_geneID_hp')
    )
    flowR3.table.goInfo <- fluidRow(DTOutput("GoInfoTable"))
    ## 1.2.211) flowR for heatmap panel ====
    flowR.hp.vbtext <- fluidRow(
        tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
        tags$hr(style="border-color: black;")
    )
    
    flowR.hp.box <- fluidRow(
        tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
        # tags$h4("Selected gene:"),
        box(title = 'Gene Information',
            width = 12,
            uiOutput('ui2'),
            tags$br(),
            tags$h4(strong('GO terms')),
            tags$style(HTML("hr {border-top: 1px solid #000000;}")),
            tags$hr(),
            fluidRow(DTOutput("GoTermTableInfobox1")),
            tags$br(),
            tags$h4("Definition for selected GO term", style = "font-weight: bold;"),
            tags$hr(),
            verbatimTextOutput(outputId = 'goTermDefi1')
        )
        
    )
    footer.all.pages <- div(class="navbar navbar-default navbar-fixed-bottom", 
                             div(
                               fluidRow(column(12,div(style = "height:10px; background-color: grey100"))),
                               fluidRow(
                                 column(8,
                                        class="container",
                                        tags$strong( 
                                          HTML('&ensp;&ensp;&ensp;&ensp;'),
                                          (HTML("&copy; 2020 - ")), 
                                          "(",
                                          a(href="http://www.dkfz.de/en/molekulare-rna-biologie/projects.html", "Diederichs Lab/Projects - DKFZ Heidelberg", style = "color: #2950a3", target="_blank"),
                                          ")"
                                        )
                                 ),
                                 column(4,
                                        # offset = 2,
                                        tags$a(
                                          img(height = 30, src = "image/DKFZ_Logo_Small.png",class = "pull-right"),
                                          # imageOutput('logoImage2'),
                                          href = "http://www.dkfz.de",
                                          target = "_blank"
                                        )		
                                 )	   	   
                               )
                             ))
    
    ## 1.2.22) body itself ====
    body <- dashboardBody(
        tags$script(HTML('
        $(document).ready(function() {
          $("header").find("nav").append(\'<span style="font-size: 20px;padding-top:0px; background: inherit;" class="myClass"><strong><font color="white"> circ2go Database  </font></strong></span>\');
        })
       ')),
        tags$head(tags$style(HTML('.form-group, .selectize-control {margin: 0px}
                                               .box-body {padding-left: 30px;padding-right: 30px;}'))),
        tabItems(
            ## tab main page ====
            tabItem(tabName = "mainPage", #1
                    fluidPage(
                      fluidRow(column(12,
                                      br(),
                                      fluidRow(column(12, h1(strong("Welcome to circ2GO"), style = "color: #0073B7;font-size:40px")), align = "left"),
                                      br(),
                                      fluidRow(column(12,h4(HTML("
                  Thousands of circular RNAs(circRNAs) have recently been shown to be expressed in eukaryotic cells. Here we provide a dataset of circRNA expression profiles for more than 50 LUAD(Lung adenocarcinoma) cell lines. You can explore circRNA datasets, view visualized circRNAs expressed on each gene, check exression level of circRNAs for each gene, and download the custom data to discover circRNAs interested. "))), align = "left"),
                                      br(),
                                      br(),
                                      fluidRow(column(5, offset = 2,  
                                                      fluidRow(
                                                        valueBox(
                                                          tags$p(HTML("<h5>Cell lines</h5>")),
                                                          subtitle = div(
                                                            fluidRow(
                                                              column(6, h3("60", align = "left"))
                                                            ),	
                                                            style = "text-align: right;"
                                                          ), 
                                                          color = "blue",
                                                          width = 8
                                                        )
                                                      )),
                                               column(5, 
                                                      fluidRow(
                                                        valueBox(
                                                          tags$p(HTML("<h5>circRNA candidates</h5>")),
                                                          subtitle = div(
                                                            fluidRow(
                                                              column(6, h3("148806", align = "left"))
                                                            ),	
                                                            style = "text-align: right;"
                                                          ), 
                                                          color = "green",
                                                          width = 8
                                                        )
                                                      ))
                                      ), # fluidRow with valueBoxes
                                      br(),
                                      br(),
                                      br(),
                                      br(),
                                      
                                      fluidRow(column(12, h4(HTML("<strong>References</strong>"), style = "color: #000000"))),
                                      fluidRow(column(12,h4(HTML("The Circular RNA Landscape of Non-Small Cell Lung Cancer Cells. <a href='https://www.mdpi.com/2072-6694/12/5/1091' target ='_blank'>Nele Van Der Steen <i>et al.</i> Cancers 2020</a>.")))),
                                      br(),
                                      br(),
                                      br(),
                                      br(),
                                      br(),                    
                                      imageOutput('logoImage'),
                                      imageOutput('logoImage2'),
                                      img(src = "image/Selection_046.png"),
                                      img(src = "image/Selection_042.png")
                                      # img(src = "image/circ2go2.png", height = 60)
                      ))
                    ) # end of fluidPage
            ), # end tabItem #1
            
  
            tabItem(tabName = 'tabAbout',
                    fluidRow(column(12, h1(strong("About circ2GO"), style = "color: #0073B7;font-size:35px")), align = "left"),
                    box( width = NULL, solidHeader = TRUE, 
                         # title="About circ2GO",
                         fluidRow(column(12,h4(HTML("<strong>How to cite</strong>")))),
                         fluidRow(column(12,h4(HTML("The Circular RNA Landscape of Non-Small Cell Lung Cancer Cells. <a href='https://www.mdpi.com/2072-6694/12/5/1091' target ='_blank'>Nele Van Der Steen <i>et al.</i> Cancers 2020</a>.")))),
                         br(),
                         fluidRow(column(12,h4(HTML("<strong>Related publications</strong>")))),
                         br(),
                         br(),
                         fluidRow(column(12,h4(HTML("<strong>Related database</strong>")))),
                         br(),
                         br(),
                         fluidRow(column(12,h4(HTML("<strong>Contact</strong>")))),
                         fluidRow(column(12, h4("The authors of the circ2GO database can be contacted at the following email address") )),
                         a(href="mailto:databases.diederichslab@dkfz.de","databases.diederichslab(at)dkfz.de",style="font-size: 17px; color: #2950a3"),
                         br(),
                         br()
                    ) # box
                    
                    ),
            ## tab circ span ====
            tabItem(tabName = "circSpan",
                    h2("circRNA Transcript Map", style = "font-weight: bold;"),
                    fluidRow(
                        box(width = 12,
                            textInput(inputId = "textInput.geneName",
                                      width = '250px',
                                      # width = "100%",
                                      label = "gene name/gene ID/transcript ID:",
                                      # value = "ENSG00000135776",
                                      placeholder = "eg: ENSG00000135776, ABCB10 "
                                      ),
                            actionButton(inputId = 'ButtonGeneName',
                                     label = 'Submit'
                                     ))),
                    fluidRow(
                      box(title = 'Gene Information',
                          width = 12,
                          uiOutput('ui1'),
                          
                          tags$br(),
                          tags$h4(strong('GO terms')),
                          tags$style(HTML("hr {border-top: 1px solid #000000;}")),
                          tags$hr(),
                          fluidRow(DTOutput("GoTermTableInfobox2")),
                          tags$br(),
                          tags$h4("Definition for selected GO term", style = "font-weight: bold;"),
                          tags$hr(),
                          verbatimTextOutput(outputId = 'goTermDefi2')
                          
                          
                          
                          
                          
                          )),
                    fluidRow(box(width = 12,
                                 fluidRow(tags$h4("circRNA Transcript Map")),
                                 fluidRow(imageOutput("myImage",
                                                      height = '1000px',
                                                      click = "image_click",
                                                      hover = hoverOpts(
                                                        id = "image_hover",
                                                        delay = 0,
                                                        delayType = "throttle"),))))
                    
                    ), # tabItem circSpan
            
            ## tab circ heatmap ====
            tabItem(tabName = 'TabcircHp',
                    h2("circRNA Heatmap"),
                    
                    fluidRow(box(width = 12,
                        tags$head(includeScript("./www/returnClick_circHp.js")),
                        fluidRow(
                          textInput(inputId = "textInput_eneName_heatmap", 
                                    width = '250px',
                                    label = "gene name/gene ID/transcript ID:",
                                    placeholder = "eg: ENSG00000135776, ABCB10 "
                                    )), 
                        fluidRow(actionButton(inputId = 'ButtonGeneNameHeatmap',label = 'Submit')),
                        )),
                    flowR.hp.box,
                    fluidRow(box(width =12,
                      fluidRow(plotlyOutput('plotlyHp1')),
                      fluidRow(plotlyOutput('plotlyHp2')),
                      fluidRow(downloadButton('DL_circ', 'Download data circ(sum)'),
                               downloadButton('DL_circ_cellLines', 'Download data circ(cellLines)')),
                      tags$br(),
                      tags$br()
                      )),
                    ),
            ## tab go search ====
            tabItem(tabName = "goSearch",
                    h2("GO Search", style = "font-weight: bold;"),
                    fluidRow(
                      box(width = 12,
                        flowR1.GoSelectButton,
                        tags$style("#actBGoTerm {margin-left:10px;}"),  
                        fluidRow(actionButton(inputId = 'actBGoTerm', label = 'submit'))
                        )
                    ),
                    fluidRow(
                      box(width = 12,
                        flowR2.table.goTerm,
                        flowR2.verbtext.goTerm.defi,
                      )
                    ),
                    fluidRow(box(width = 12,
                                 fluidRow(h4('Genes for selected GO term', style = "font-weight: bold;")),
                                 flowR3.table.goInfo)),
                    fluidRow(
                      box(width = 12,
                        fluidRow(flowR301.goHP),
                        fluidRow(column(12,
                                        actionButton(inputId = 'GoGeneIDcm', label = 'circRNA map'), 
                                        actionButton(inputId = 'GoGeneIDHp', label = 'Heatmap'), 
                                        downloadButton('DL_goInfo', 'Download data'))),
                        tags$br(),
                        tags$br()
                        ))
                    ),
            ## tab download ====
            tabItem(tabName = 'tabDownload', 
                    h2('Download Data by cell lines'),
                    fluidRow(box(width = 12,
                                 fluidRow(selectInput(inputId = "cellLineSelect", label = "Cell lines selection",
                                                      width = '400px',
                                                      choices = cellLineNames,
                                                      multiple = TRUE)),
                                 fluidRow(radioButtons(inputId = "radioConv", label = "gene/circRNA data",
                                                       choices = c("circRNA level" = "circl",
                                                                   "Gene level" = "genel")
                                 )),
                                 # fluidRow(checkboxInput(inputId = "GOtermIn", label = "GO term included", value = TRUE)),
                                 fluidRow(downloadButton(outputId = 'dlbulk1', 'Download data'))
                                 )
                             ),
                    h2('Download Data by genes'),
                    fluidRow(box(width = 12,
                                 textAreaInput(inputId = "textAi1", 
                                               label = "Gene names", 
                                               value = "", 
                                               width = "300px",
                                               height = '400px'),
                                 tags$br(),
                                 fluidRow(radioButtons(inputId = "radioConv1", label = "gene/circRNA data",
                                                       choices = c("circRNA level" = "circl",
                                                                   "Gene level" = "genel")
                                 )),
                                 fluidRow(downloadButton(outputId = 'dlbulk2', 'Download data')),
                                 tags$br(),
                                 tags$br()
                                 ))
                    
                    ),
            ### tab help ====
            tabItem(tabName = "tabHelp", #16
                    fluidRow(column(12, h1(strong("Help and Documentation"), style = "color: #0073B7;font-size:35px")), align = "left"),
                    box( width = '1000px', 
                         # status = "info", 
                         solidHeader = TRUE, 
                         # title="Help and Documentation",
                         br(),
                         fluidRow(column(12,h4(HTML("<ul><li>circ2GO is a database that provide a complete circRNA expression profile for 60 cell lines, which contain 50 LUAD cell lines, 7 non-LUAD cell lines from lung cancer tissue, 3 cell lines from normal lung tissue. "))), align = "left"),
                         fluidRow(column(12,h4(HTML("<ul><li>The <strong>circRNA Search</strong> options can be used to search the database for single gene's circRNA expression profile. The <strong>circRNA transcript map</strong> shows circRNA location on transcripts visually. The <strong>circRNA heatmap</strong> shows circRNA read counts in each cell lines. Hovering mouse on the heatmap will give the read count value and cell line names.The read count data in heatmap of searched gene could be downloaded in CSV file.</li></ul>"))), align = "left"),
                         fluidRow(column(12,h4(HTML("<ul><li>The <strong>GO Search</strong> option can be used to search the database based on GO terms. With a list of genes searched out, users can get heatmap and visualized circRNA map respectively.</li></ul>"))), align = "left"),
                         fluidRow(column(12,h4(HTML("<ul><li>The dataset could be download bulkly the option <strong>Data Dwnload</strong>. Data for cell lines can be downloaded selectively</li></ul>"))), align = "left"),
                         br(),
                         br()
                         )
                    ) # end tabItem #16    	 
            ),# tabItems
        

        
        footer.all.pages
        )
    
    
    # 1.3) ui ====
    ui <- dashboardPage(
        dashboardHeader(title = img(src = "image/circ2go2.png", height = 60)),
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
                                              end='')
        )
        
        
        ## 140) mainPage ====
        output$logoImage <- renderImage({
            outfile <- './www/image/circ2GO.png'
            list(src = outfile,
                 contentType = 'image/png',
                 width = 400,
                 height = 300,
                 alt = "This is alternate text")
            }, deleteFile = FALSE)
        output$logoImage2 <- renderImage({
          outfile <- './www/image/circ2logo.png'
          list(src = outfile,
               contentType = 'image/png',
               width = 400,
               height = 300,
               alt = "This is alternate text")
        }, deleteFile = FALSE)
        
        
        ## 141) circSpan ====

        ## 1412) get new geneName ====
        observeEvent(input$ButtonGeneName,{
            vals$mpV1 <- isolate(input$textInput.geneName)
            vals$boxInfo <- get_box_info(vals$mpV1)
            
            a <- vals$boxInfo$geneID 
            d <- datasets$df.go
            e <- d[d$Gene.stable.ID == a,c(9:12)]
            e <- e[order(e$GO.term.name),]
            colnames(e) <- gsub('\\.', ' ', colnames(e))
            # edit
            vals$goTable2 <- e
            output$GoTermTableInfobox2 <- renderDT(e[,c(-2)], selection = 'single', rownames = FALSE)
        })
            
        ## 14121) plot circ map ====
        observeEvent(input$ButtonGeneName, {
            wait()
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
        
        ## 14122) box info table click ====
        observeEvent(input$GoTermTableInfobox2_cell_clicked, {
          list.tb1.select <- isolate(input$GoTermTableInfobox2_cell_clicked)
          if (length(list.tb1.select) >0 ){
            go.term.index <- list.tb1.select[[1]]
            output$goTermDefi2 <- renderText(as.character(vals$goTable2[go.term.index, 2]))
          }
        })
        
        ## 1413) output click coordination ====
        observeEvent(input$image_click, {
            x <- isolate(input$image_click$x)
            y <- isolate(input$image_click$y)
            vals$co.default <- c(x,y)
            teID <- get_transID(vals$co.default)
            if (teID == 'empty'){return()}
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
            wait()
            removeModal()
            vals$plotlys <- heatmap_geneID(get_transID(vals$co.default))
            updateTabItems(session, "sideBarID", 'TabcircHp')
            # output$hp.vb.geneID <- renderText(get_transID(vals$co.default))
        })
        
        output$co_xy <- renderText({
            vals$co.default
        })
        output$ID_text <- renderText(get_transID(vals$co.default))
        
        

        
        ## 142) Go Table ====
        observeEvent(input$actBGoTerm, {
            x.s <- grep(isolate(input$autoC1), x.long, value=TRUE)
            x.s1 <- grep(isolate(input$autoC2), x.s, value=TRUE)
            vals$tb.term <- y4[y4$GO.term.name %in% x.s1,]
            colnames(vals$tb.term) <- gsub('\\.', ' ', colnames(vals$tb.term))
            })
        
        output$GoTermTable <- renderDT(vals$tb.term[,c(-3)], selection = 'single', rownames = FALSE)
        
        # output$GoTermTable <- renderDT(iris, selection = 'single', rownames=FALSE)
      
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
                colnames(vals$tb.info) <- gsub('\\.', ' ', colnames(vals$tb.info))
                }
            })
        output$goTermDefi <- renderText(vals$go.term.defi)
        output$GoInfoTable <- renderDT(vals$tb.info, selection = 'single', rownames = FALSE)
        
        ##1423)table GO heatmap ====
        
        observeEvent(input$GoInfoTable_cell_clicked,{
            list.tb2.select <- isolate(input$GoInfoTable_cell_clicked)
            if(length(list.tb2.select) > 0){
                go.info.index <- list.tb2.select[[1]]
                vals$GoInfoTable.IDselected <- vals$tb.info[go.info.index,1]
                }
            })
        
        output$go_vb_geneID_hp <- renderText(vals$GoInfoTable.IDselected)
        
        observeEvent(input$GoGeneIDHp,{
            wait()
            print('GoGeneIDHpButton:')
            vals$boxInfo <- get_box_info(vals$GoInfoTable.IDselected)
            vals$plotlys <- heatmap_geneID(vals$GoInfoTable.IDselected)
            updateTabItems(session, "sideBarID", 'TabcircHp')
            })
        ##1423) table GO circMap ====
        observeEvent(input$GoGeneIDcm,{
          wait()
          vals$boxInfo <- get_box_info(vals$GoInfoTable.IDselected)
          data.geneID <- plot_circSpan_1(datasets, vals$GoInfoTable.IDselected)
          if (nrow(data.geneID$df.geneID.transcripts) == 0){
            return(NULL)
          }
          outfile <- data.geneID$file.output
          png1 <- readPNG(outfile)
          width1 <- dim(png1)[2]
          height1 <- dim(png1)[1]
          output$myImage <- renderImage({
            print(paste('height: ',height1))
            print(width1)
            list(src = outfile,
                 # contentType = 'image/svg+xml',
                 contentType = 'image/png',
                 width = width1,
                 height = height1,
                 alt = "This is alternate text")
          }, deleteFile = TRUE)
          updateTabItems(session, "sideBarID", 'circSpan')
        })
        

        
        ## 1424) download table ====
        output$DL_goInfo <- downloadHandler(
          filename = function() { 
            paste("dataset-", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            write.csv(vals$tb.info, file)
          })
        
        ## 143) heatmap geneID ====
        ## 1430) infobox/goTable ====
        observeEvent(input$ButtonGeneNameHeatmap, {
            wait()
            vals$mpV1 <- isolate(input$textInput_eneName_heatmap)
            vals$boxInfo <- get_box_info(vals$mpV1)
            a <- vals$boxInfo$geneID 
            d <- datasets$df.go
            e <- d[d$Gene.stable.ID == a,c(9,10,11,12)]
            e <- e[order(e$GO.term.name),]
            colnames(e) <- gsub('\\.', ' ', colnames(e))
            vals$goTable1 <- e
            output$GoTermTableInfobox1 <- renderDT(e[,c(1,3,4)], selection = 'single', rownames = FALSE)
            
            vals$plotlys <- heatmap_geneID(vals$mpV1)
        })
        
        ## 1431) goTable defination ====
        # GoTermTable_cell_clicked
        
        observeEvent(input$GoTermTableInfobox1_cell_clicked, {
          list.tb1.select <- isolate(input$GoTermTableInfobox1_cell_clicked)
          if (length(list.tb1.select) >0 ){
            go.term.index <- list.tb1.select[[1]]
            output$goTermDefi1 <- renderText(as.character(vals$goTable1[go.term.index, 2]))}
        })
        
        ## 1432) plots ====
        output$plotlyHp1 <- renderPlotly(vals$plotlys$p1)
        output$plotlyHp2 <- renderPlotly(vals$plotlys$p2)
        ## 1433) downloader ====
        output$DL_circ <- downloadHandler(
          filename = function() { 
            paste("dataset_", isolate(input$textInput_eneName_heatmap), ".csv", sep="")
          },
          content = function(file) {
            df <- readRDS('./0_data/df.geneID.circ.rdata')
            write.csv(df, file)
          })
        
        output$DL_circ_cellLines <- downloadHandler(
          filename = function() { 
            paste("dataset_", isolate(input$textInput_eneName_heatmap), ".csv", sep="")
          },
          content = function(file) {
            df <- readRDS('./0_data/df.geneID.circ.cellLines.rdata')
            write.csv(df, file)
          })
        ## 144) bulk downlaod ====
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
            write.csv(e2[1:100,], file)
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
    } # end of server
    ## return ui/server ====
    list1 <- list(ui = ui, 
                  server = server)
}