## 0.0) common values: ====
df.hml <- list(type = 'hp')

## 0) function: ====
inline = function (x) {
  tags$div(style="display:inline-block;", x)
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
                  chrID=paste(chrID, ':', sep = ''),
                  strand=paste('(',strand,')Strand(GRCh38)', sep = ''),
                  start=paste(start,'-',sep = ''),
                  end=paste(end,',',sep=''))
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
    x.short <- c(x.short,'  ')
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
            menuItem(text = "circ2GO Main Page", tabName = "mainPage", selected = TRUE),
            menuItem(text = "circRNA Search", icon = icon("bars"),startExpanded = TRUE,
                     menuSubItem(text = 'circRNA Transcript Map',
                                 tabName = 'circSpan'),
                     menuSubItem(text = 'circRNA Heatmap',
                                 tabName = 'TabcircHp')
                     ),
            menuItem(text = "GO Search", tabName = 'goSearch', icon = icon("bars")),
            menuItem(text = "circRNA-miRNA Search", tabName = 'tabMicro', icon = icon("bars")),
            menuItem(text = 'Data Download', tabName = 'tabDownload', icon = icon('bars')),
            menuItem(text = 'Help', tabName = 'tabHelp'),
            menuItem(text = 'About', tabName = 'tabAbout')
            
        )
    )
    ## 1.2.2) body ====
    ## 1.2.21) flowR for GO panel ====
    
    selip <- function(x, label){
      selectizeInput(inputId = x,
                     label = label,
                     choices = x.short,
                     selected = '  ',
                     size = 10,
                     width = '250px',
                     multiple = FALSE, # allow for multiple inputs
                     options = list(create = FALSE) 
      )}
    
    flowR1.GoSelectButton <- fluidRow(
      inline(selip('autoC1', 'Key Word 1')),
      tags$style('#and1 {padding-top:0px; font-weight: bold;}'),
      inline(HTML('<p id = "and1">&nbspAND&nbsp<br><span style = "font-size:10px">&nbsp</span></p>')),
      inline(selip('autoC2', 'Key Word 2'))
    )
    flowR1.GoSelectButton1 <- fluidRow(
      column(3, 
             selectizeInput(inputId = 'autoC1',
                            label = 'Key word 1',
                            choices = x.short,
                            selected = '  ',
                            size = 10,
                            multiple = FALSE, # allow for multiple inputs
                            options = list(create = FALSE),
                            width = '250px')
             ),
      column(1,
             tags$style('#and1 {margin-top:30px; font-weight: bold;}'),
             HTML('<p id = "and1">&nbspAND&nbsp</p>')
             ),
      column(3,
             selectizeInput(inputId = 'autoC2',
                              label = 'Key word 2',
                              choices = x.short,
                              selected = '  ',
                              size = 10,
                              multiple = FALSE, # allow for multiple inputs
                              options = list(create = FALSE),
                              width = '250px')
             )
      ) # flowR1.GoSelectButton

    
    
    flowR3.0.goHP <- fluidRow(
       actionButton(inputId = 'GoGeneIDHp', label = 'Heatmap')
    )
    
    flowR301.goHP <- fluidRow(
        tags$h4("circMap/Heatmap/Download for selected gene:", style = "font-weight: bold;padding-left:10px"),
        verbatimTextOutput(outputId = 'go_vb_geneID_hp')
    )
    ## 1.2.211) flowR for heatmap panel ====
    flowR.hp.vbtext <- fluidRow(
        tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
        tags$hr(style="border-color: black;")
    )
    
    flowR.hp.box <- fluidRow(
        tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
        box(title = 'Gene Information',
            width = 12,
            uiOutput('uiInfoboxHp'),
            tags$br()
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
                                          href = "http://www.dkfz.de",
                                          target = "_blank"
                                        )		
                                 )	   	   
                               )
                             ))
    
    ## 1.2.22) body itself ====
    body <- dashboardBody(
       #  tags$script(HTML('
       #  $(document).ready(function() {
       #    $("header").find("nav").append(\'<span style="font-size: 20px;padding-top:0px; background: inherit;" class="myClass"><strong><font color="white"> circ2go Database  </font></strong></span>\');
       #  })
       # ')),
        tags$head(tags$style(HTML('.form-group, .selectize-control {margin: 0px}
                                               .box-body {padding-left: 30px;padding-right: 30px;}'))),
        tabItems(
            ## tab main page ====
            tabItem(tabName = "mainPage", #1
                    fluidPage(
                      fluidRow(column(12,
                                      # color: #0073B7
                                      tags$style(HTML("#logoImage {text-align: center;}")),
                                      imageOutput('logoImage'),
                                      fluidRow(column(12, h1(strong("Welcome to circ2GO"), style = "color: black;font-size:40px")), align = "left"),
                                      fluidRow(column(12,h4(HTML("
                                      Thousands of circular RNAs (circRNAs) are expressed in eukaryotic cells. Here, we provide a dataset of circRNA expression for 148811 circRNAs based on 3.5 billion sequencing reads of rRNA-depleted RNA derived from 60 human lung cell lines including 50 lung adenocarcinoma (LUAD) cell lines. Additionally, the circRNA expression data is merged with orthogonal data about gene and transcript information, functional annotation by gene ontology and microRNA binding site prediction.
                                      <br><br>In circ2GO, you can explore the circRNA profiles, visualize circRNAs expressed from each gene and differentially spliced transcript, view the gene ontology annotation and microRNA binding sites for each circRNA and download the respective data. Most importantly, reverse searches allow the search for circRNAs derived from genes with certain biological processes, molecular functions or cellular components assigned by gene ontology or possessing specific microRNA binding sites to discover circRNAs in the areas of your interest.
                                                                 "))), align = "left"),
                                      br(),
                                      br(),
                                      fluidRow(column(4, offset = 2,  
                                                      fluidRow(
                                                        valueBox(
                                                          tags$p(HTML("<h5>Genes</h5>")),
                                                          subtitle = div(
                                                            fluidRow(
                                                              column(6, h3("12251", align = "left"))
                                                            ),	
                                                            style = "text-align: right;"
                                                          ), 
                                                          color = "green",
                                                          width = 10
                                                        )
                                                      )),
                                               column(4, 
                                                      fluidRow(
                                                        valueBox(
                                                          tags$p(HTML("<h5>circRNAs</h5>")),
                                                          subtitle = div(
                                                            fluidRow(
                                                              column(6, h3("148811", align = "left"))
                                                            ),	
                                                            style = "text-align: right;"
                                                          ), 
                                                          color = "green",
                                                          width = 10
                                                        )
                                                      ))
                                               # column(4, 
                                               #        fluidRow(
                                               #          valueBox(
                                               #            tags$p(HTML("<h5>miRNAs</h5>")),
                                               #            subtitle = div(
                                               #              fluidRow(
                                               #                column(6, h3("897", align = "left"))
                                               #              ),	
                                               #              style = "text-align: right;"
                                               #            ), 
                                               #            color = "green",
                                               #            width = 10
                                               #          )
                                               #        ))
                                      ), # fluidRow with valueBoxes
                                      br(),br(),br(),br(),
                                      fluidRow(column(12, h4(HTML("<strong>References</strong>"), style = "color: #000000"))),
                                      fluidRow(column(12,h4(HTML("The Circular RNA Landscape of Non-Small Cell Lung Cancer Cells. <a href='https://www.mdpi.com/2072-6694/12/5/1091' target ='_blank'>Nele Van Der Steen <i>et al.</i> Cancers 2020</a>.")))),
                                      br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br()
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
                         # fluidRow(column(12,h4(HTML("<strong>Related publications</strong>")))),
                         br(),
                         br(),
                         # fluidRow(column(12,h4(HTML("<strong>Related database</strong>")))),
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
                            includeScript("./www/returnClick_textInputGeneName.js"),
                            textInput(inputId = "textInputGeneName",
                                      width = '250px',
                                      # width = "100%",
                                      label = "gene name/gene ID/transcript ID:",
                                      # value = "ENSG00000135776",
                                      placeholder = "eg: ENSG00000135776, ABCB10 "
                                      ),
                            br(),
                            fluidRow(actionButton(inputId = 'ButtonGeneName',
                                                  label = 'Submit', width = '80px')
                                     )
                            )),
                    fluidRow(
                      box(title = 'Gene Information',
                          width = 12,
                          uiOutput('uiCircInfobox'),
                          br(),
                          uiOutput('uiCircGotermtb'),
                          uiOutput('uiCircGotermDef')
                          )),
                    uiOutput('uiCircMap')
                    
                    ), # tabItem circSpan
            
            ## tab circ heatmap ====
            tabItem(tabName = 'TabcircHp',
                    h2("circRNA Heatmap", style = "font-weight: bold;"),
                    
                    fluidRow(box(width = 12,
                        fluidRow(
                          includeScript("./www/returnClick_textInput_geneName_heatmap.js"),
                          textInput(inputId = "textInput_geneName_heatmap", 
                                    width = '250px',
                                    label = "gene name/gene ID/transcript ID:",
                                    placeholder = "eg: ENSG00000135776, ABCB10 "
                                    )), 
                        br(),
                        fluidRow(actionButton(inputId = 'ButtonGeneNameHeatmap',label = 'Submit', width = '80px')),
                        )),
                    flowR.hp.box,
                    uiOutput('uiPlotly'),
                    fluidRow(box(width =12,
                                 fluidRow(downloadButton('DL_circ', 'Download Sum'),
                                          HTML('&nbsp&nbsp&nbsp'),
                                          downloadButton('DL_circ_cellLines', 'Download Cell Lines'))
                                 )),
                    br(),br()
                    ),
            ## tab go search ====
            tabItem(tabName = "goSearch",
                    h2("GO Search", style = "font-weight: bold;"),
                    fluidRow(box(width =12,
                                 fluidRow(h4('Basic search', style = "font-weight: bold;")),
                                 # tags$style(HTML("hr {border-top: 1px solid #000000;}")),
                                 # tags$hr(),
                                 includeScript("./www/returnClick_goSearchBasic.js"),
                                 textInput(inputId = 'goSearchBasic', label = ' ', value = '',
                                           width = '250px',
                                           placeholder = 'e.g.,protein binding, GO:000122, 122'),
                                 HTML('<i>Searching genes by entering a complete GO term</i><br><br>'),
                                 fluidRow(actionButton(inputId = 'actBGoTermBasicSubmit', label = 'Submit', width = '80px'),
                                          HTML('&nbsp&nbsp&nbsp'),
                                          actionButton(inputId = 'actBGoTermBasicClear', label = 'Reset', width = '80px')
                                          )
                                 )),
                    fluidRow(
                      box(width = 12,
                          fluidRow(h4('Advanced search', style = "font-weight: bold;")),
                          # tags$style(HTML("hr {border-top: 1px solid #000000;}")),
                          # tags$hr(),
                          flowR1.GoSelectButton, 
                          HTML('<i><strong>Step 1:</strong> Searching for a GO term by entering 2 key words<br><strong>Step 2:</strong> Searching genes by selecting a GO term from table generated below</i><br><br>'),
                          tags$style("#actBGoTermSubmit {margin-left:10px;}"),  
                          fluidRow(actionButton(inputId = 'actBGoTermSubmit', label = 'Submit', width = '80px'),
                                   HTML('&nbsp&nbsp&nbsp'),
                                   actionButton(inputId = 'actBGoTermClear', label = 'Reset', width = '80px')
                                   ),
                          tags$br(),
                          fluidRow(uiOutput('uitermtb')),
                          fluidRow(uiOutput('uiDefi'))
                        
                    )),
                    fluidRow(uiOutput('uiInfotb')),
                    fluidRow(box(width = 12,
                                 fluidRow(flowR301.goHP),
                                 fluidRow(column(12,
                                                 actionButton(inputId = 'GoGeneIDcm', label = 'circRNA map'), 
                                                 HTML('&nbsp&nbsp&nbsp'),
                                                 actionButton(inputId = 'GoGeneIDHp', label = 'Heatmap'), 
                                                 HTML('&nbsp&nbsp&nbsp')
                                 )))),
                    br(),br()
                    ), # end of tab item
            ## tab microRNA search ====
            tabItem(tabName = 'tabMicro',
                    h2("circRNA-miRNA Search", style = "font-weight: bold;"),
                    fluidRow(box(width = 12,
                                 # fluidRow(h4('Search', style = "font-weight: bold;")),
                                 # tags$style(HTML("hr {border-top: 1px solid #000000;}")),
                                 # tags$hr(),
                                 includeScript("./www/returnClick_miSearch.js"),
                                 textInput(inputId = 'miSearch', label = ' ', value = '',
                                           width = '350px',
                                           placeholder = 'e.g.,hsa-let-7a-2-3p,CEPT1,ENSG00000031698...'),
                                 br(),
                                 fluidRow(actionButton(inputId = 'miSearchSubmit', label = 'submit', width = '80px'),
                                          HTML('&nbsp&nbsp&nbsp'),
                                          actionButton(inputId = 'miSearchClear', label = 'Reset', width = '80px'),
                                          # edit
                                          br(),br(),
                                          tags$style("#nofound {margin-left:20px;}"),
                                          uiOutput('nofound')
                                 )
                                 )),
                    fluidRow(uiOutput('uiMiTb'),
                             box(width =12,
                                 tags$h4("circMap/Heatmap/Download for selected gene:", style = "font-weight: bold;"),
                                 # tags$style("#miSelectId {margin-left:10px;}"),
                                 verbatimTextOutput(outputId = 'miSelectId'),
                                 fluidRow(actionButton(inputId = 'miIdCm', label = 'circRNA map'), 
                                          HTML('&nbsp&nbsp&nbsp'),
                                          actionButton(inputId = 'miIdHp', label = 'Heatmap')
                                          ))
                             ),
                    br(),br(),br(),br(),br(),br(),br(),br()
                    ),
            ## 12221) tab download ====
            tabItem(tabName = 'tabDownload', 
                    ## 122211) circRNA data ====
                    h2('Download circRNA Data by cell lines'),
                    fluidRow(box(width = 12,
                                 fluidRow(includeScript("./www/returnClick_cellLineSelect.js"),
                                          selectInput(inputId = "cellLineSelect", label = "Cell lines selection",
                                                      width = '400px',
                                                      choices = cellLineNames,
                                                      multiple = TRUE)),
                                 fluidRow(radioButtons(inputId = "radioConv", label = "gene/circRNA data",
                                                       choices = c("circRNA level" = "circl",
                                                                   "Gene level" = "genel")
                                 )),
                                 # fluidRow(checkboxInput(inputId = "GOtermIn", label = "GO term included", value = TRUE)),
                                 fluidRow(downloadButton(outputId = 'dlbulk1', 'Download data')),
                                 br(),br()
                                 )
                             ),
                    h2('Download circRNA Data by genes'),
                    fluidRow(box(width = 12,
                                 textAreaInput(inputId = "textAi1", 
                                               label = "Gene names", 
                                               placeholder = 'IDs splited by comma/Semi-colon/tabs/spaces are acceptable.',
                                               value = "", 
                                               width = "300px",
                                               height = '400px'),
                                 tags$br(),
                                 fluidRow(radioButtons(inputId = "radioConv1", label = "gene/circRNA data",
                                                       choices = c("circRNA level" = "circl",
                                                                   "Gene level" = "genel")
                                 )),
                                 fluidRow(downloadButton(outputId = 'dlbulk2', 'Download data')),
                                 br(),br()
                                 )),
                    ## 122212) miRNA data ====
                    h2('Download miRNA Data'),
                    fluidRow(box(width = 12,
                                 textAreaInput(inputId = "miMultText", 
                                               label = "Gene names/circRNA ID/miBase ID/miRNA family ID", 
                                               placeholder = 'IDs splited by comma/Semi-colon/tabs/spaces are acceptable.',
                                               value = "", 
                                               width = "300px",
                                               height = '400px'),
                                 tags$br(),
                                 fluidRow(downloadButton(outputId = 'miDlIds', 'Download data by IDs'),
                                          HTML('&nbsp&nbsp&nbsp'),
                                          downloadButton(outputId = 'miDlAll', 'Download all data')
                                          ),
                                 tags$br(),
                                 tags$br()
                    )),
                    br(),br(),br()
                    
                    ),
            ### tab help ====
            tabItem(tabName = "tabHelp", #16
                    fluidRow(column(12, h1(strong("Help and Documentation"), style = "color: #0073B7;font-size:35px")), align = "left"),
                    box( width = '1000px', 
                         # status = "info", 
                         solidHeader = TRUE, 
                         # title="Help and Documentation",
                         br(),
                         fluidRow(column(12,h4(HTML("<ul><li><strong>circ2GO</strong> is a database that provides a comprehensive circRNA expression profiles for 60 human cell lines comprising 50 lung adenocarcinoma (LUAD), seven non-LUAD non-small cell lung cancer (NSCLC) and three non-transformed lung cell lines. Importantly, the circ2GO data was derived from 3.5 billion RNA-seq reads from rRNA-depleted RNA and hence not biased due to polyA enrichment.
                                                    </li></ul>"))), align = "left"),
                         fluidRow(column(12,h4(HTML("<ul><li>The <strong>circRNA Search</strong> option can be used to search the database for circRNA expression profiles of a single gene. The <strong>circRNA Transcript Map</strong> visualizes the location and expression of all circRNAs on all transcripts of the respective gene. The <strong>circRNA Heatmap</strong> depicts the circRNA read counts in each cell line. Hovering with the mouse over the heatmap gives the read count value (normalized to library size) and cell line name. The read count data shown in the heatmap can be downloaded as a CSV file.
                                                    </li></ul>"))), align = "left"),
                         fluidRow(column(12,h4(HTML("<ul><li>The <strong>GO Search</strong> option allows the search of the database for the gene ontology terms assigned to the gene from which the circRNA is derived. Thereby, the user can search for circRNAs potentially linked to the <strong>biological process, molecular function or cellular component</strong> of interest. For the resulting list of circRNAs, a transcript map and heatmap can be obtained.
                                                    </li></ul>"))), align = "left"),
                         fluidRow(column(12,h4(HTML("<ul><li>The <strong>circRNA-miRNA search </strong> option provide a database of  circRNA-miRNA pairs predicted by <a href = ‘http://www.targetscan.org/vert_72/’>TargetScan </a> and <a href=’https://genomebiology.biomedcentral.com/articles/10.1186/gb-2003-5-1-r1’>miRanda</a>. It covers 897 high-confidence human miRNAs from <a href='http://www.mirbase.org/'>miRbase</a> and 25166 high abundant circRNAs from our circRNA database. Users could get miRNAs as potential binding partners of a circRNA, and vice versa. The counts of circRNA-miRNA  binding sites predicted byTargetSscan and miRanda  are also included in results. Users can use circRNA ID/gene ID/gene name/MiRbase ID/miRNA family ID when searching, fuzzy search is supported.</li></ul>"))), align = "left"),
                         fluidRow(column(12,h4(HTML("<ul><li>With the <strong>Data Download</strong> option, batch searches for cell lines or genes can be carried out and the resulting datasets can be downloaded as CSV files.</li></ul>"))), align = "left"),
                         br(),
                         br()
                         )
                    ) # end tabItem #16    	 
            ),# tabItems
        

        
        footer.all.pages
        )
    
    
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
        source('./1_transcript_coordinate_scale/library.r')
        source('./1_transcript_coordinate_scale/1_get_circRNA_transID_data.r')
        source('./1_transcript_coordinate_scale/2_data_for_geneID.r')
        source('./2_plot_in_frame/1_plot_in_frame.R')
        
      
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
    ## return ui/server ====
    list1 <- list(ui = ui, 
                  server = server)
}