dir.base <- '/home/lyuya/shiny.server/circ2go'
setwd(dir.base)

## 1) repare big datasets ====

## 2) running shiny ====
source('./3_shiny_frame/3.1_shiny_frame.r')
us <- get_shiny_modules()
shinyApp(us$ui, us$server)
