
## 0) function: ====



get_transID <- function(vals){
    ## fetch df.th : ====
    # df.th <- readRDS(file = '~/yan150/report_work_weekly/week2020_16_online_panel/2_shiny/3_rstudio_shiny_project/2_plot_in_frame/df.th.rdata')
    df.th <- readRDS(file = paste(dir.base, '2_plot_in_frame/df.th.rdata', sep = '/'))
    co.start <- readRDS(file = paste(dir.base, '2_plot_in_frame/co.start.rdata', sep = '/'))
    p.height.px <- readRDS(file = paste(dir.base, './2_plot_in_frame/p.height.px.rdata', sep = '/'))
    padding <- 30
    row.px <- (p.height.px - padding*2)/(nrow(df.th)+1)
    x <- vals$co.default[1]
    y <- vals$co.default[2]
    row.index <- (y - padding)%/%row.px
    ID <- df.th[row.index,]$ID
    
    if (co.start > x & x> padding){
        return(ID)
    }else{
        return('teID')
    }
}

## 1) func: get_shiny_modules ====
get_shiny_modules <- function(){
    ## 1.1) common values ====
    outfile <- paste(dir.base, '2_plot_in_frame/test.png', sep = '/')
    source('./4_go_search/1_go_search.r')
    df.go <- datasets$df.go
    x <- df.go$GO.term.name
    x.short <- unique(split_vector(x, sep = ' '))
    x.long <- unique(x)
    
    y <- df.go[,c(8:12)]
    y2 <- y[y$GO.term.name != '',]
    y4 <- y2[!duplicated(y2$GO.term.name),]
    

    ## 1.2) modules for panel ====
    ## 1.2.1) sidebar ====
    sidebar <- dashboardSidebar(
        sidebarSearchForm(textId = "searchText", buttonId = "searchButton",
                          label = "Search..."),
        sidebarMenu(
            menuItem(text = "Main Page", tabName = "mainPage", selected = TRUE),
            menuItem(text = "circRNA Span", tabName = "circSpan", icon = icon("dashboard")),
            menuItem(text = "GO Search", tabName = 'goSearch', icon = icon("th"), 
                     badgeLabel = "new", badgeColor = "green"),
            menuItem(text = "link", icon = icon("file-code-o"), 
                     href = "https://github.com/rstudio/shinydashboard/")
        )
    )
    ## 1.2.2) body ====
    ## 1.2.21) flowR for GO panel ====
    flowR1.GoSelectButton <- fluidRow(
        column(3,
               selectizeInput(inputId = 'autoC1',
                              label = 'Search',
                              choices = x.short,
                              selected = '',
                              size = 10,
                              multiple = FALSE, # allow for multiple inputs
                              options = list(create = FALSE)) # if TRUE, allows newly created inputs
        ),
        
        
        column(3,
               selectizeInput(inputId = 'autoC2',
                              label = 'Search',
                              choices = x.short,
                              selected = '',
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
        tags$h4("Selected Go term defination"),
        verbatimTextOutput(outputId = 'goTermDefi'
                           # placeholder = 'Go term defination'
                           ),
        tags$hr(style="border-color: black;")
    )
    
    flowR3.table.goInfo <- fluidRow(
        dataTableOutput("GoInfoTable")
    )
    
    
    ## 1.2.22) body itself ====
    body <- dashboardBody(
        tabItems(
            tabItem(tabName = "mainPage",
                    h2("first tab should be shown.")
            ),
            tabItem(tabName = "circSpan",
                    h2("Dashboard tab content"),
                    tags$p("p creates a paragraph of text."),
                    
                    fluidRow(
                        textInput(inputId = "textInput.geneName", 
                                  label = "Gene name:",
                                  value = "ENSG00000135776.4", 
                                  width = '250px'
                        ),
                        actionButton(inputId = 'Button.geneName',
                                     label = 'Go'
                                     )
                    ),
                    
                    helpText("When you click the button above, you should see",
                             "the output below update to reflect the value you",
                             "entered at the top:"),
                    textOutput("textOutput.geneName"),
                    tags$hr(style="border-color: black;"),
                    textOutput("co_xy"),
                    textOutput("ID_text"),
                    imageOutput("myImage",
                                click = "image_click",
                                hover = hoverOpts(
                                    id = "image_hover",
                                    delay = 0,
                                    delayType = "throttle"
                                ),)
            ),
            tabItem(tabName = "goSearch",
                    h2("Widgets tab content"),
                    flowR1.GoSelectButton,
                    flowR2.table.goTerm,
                    flowR2.verbtext.goTerm.defi,
                    flowR3.table.goInfo
            )
        )
    )
    
    
    # 1.3) ui ====
    ui <- dashboardPage(
        dashboardHeader(title = "Simple tabs"),
        sidebar,
        body
    )
    ## 1.4) server ====
    server <- function(input, output, session) {
        ## 141) circSpan ====
        ## 1411)default values ====
        vals <- reactiveValues(co.default = c(0,0), geneName = 'ENSG00000135776.4')
        
        ## 1412) get new geneName ====
        observeEvent(input$Button.geneName, {
            vals$geneName <- input$textInput.geneName
            data.geneID<- plot_circSpan_1(datasets, vals$geneName)
            output$textOutput.geneName <- renderPrint({vals$geneName})
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
            }, deleteFile = FALSE)
        })
        
        ## 1413) output click coordination ====
        observeEvent(input$image_click, {
            x <- input$image_click$x
            y <- input$image_click$y
            vals$co.default <- c(x,y)

        })
        output$co_xy <- renderText({
            vals$co.default
        })
        output$ID_text <- renderText(get_transID(vals))
        
        
        ## 1414) plot circSpan ====
        output$image_clickinfo <- renderPrint({
            cat("Click:\n")
            str(input$image_click)
        })
        
        ## 142) Go Table ====
        ##1421) flow2: table go terms ====    
        tbl <- eventReactive(input$actB.goTerm, {
            x.s <- grep(input$autoC1, x.long, value=TRUE)
            x.s1 <- grep(input$autoC2, x.s, value=TRUE)
            y5 <- y4[y4$GO.term.name %in% x.s1,]
            return(y5)
        })
        output$GoTermTable <- renderDataTable({
            tbl()[,c(-3)]
        }, selection = 'single')
        
        
        ##1422) flow3: table GO info ====
        tbl2 <- eventReactive(input$GoTermTable_cell_clicked, {
            req(length(input$GoTermTable_cell_clicked) > 0)
            go.term.index <- input$GoTermTable_cell_clicked[[1]]
            go.term.name <- tbl()[go.term.index,2]
            go.term.defi <- tbl()[go.term.index,3]
            
            x1 <- df.go[df.go$GO.term.name == go.term.name,]
            x11 <- unique(x1$Gene.stable.ID)
            x2 <- datasets$df.circ.conv.lite
            x3 <- x2[x2$geneID %in% x11,]
            list1 <- list(go.term.defi = go.term.defi,
                          df.circ = x3
            )
            return(list1)
        })
        output$goTermDefi <- renderText(tbl2()$go.term.defi)
        
        output$GoInfoTable <- renderDataTable({
            tbl2()$df.circ
        }, selection = 'single')
        

    }
    ## return ui/server ====
    list1 <- list(ui = ui, 
                  server = server)
}


