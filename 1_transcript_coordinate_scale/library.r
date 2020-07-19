# lib.loc <- '~/R/x86_64-pc-linux-gnu-library/4.0'
# lib.loc <- '~/R/x86_64-redhat-linux-gnu-library/3.6'
# firstPath <- '/home/lyuya/R/x86_64-redhat-linux-gnu-library/3.6'
# .libPaths() <- c(firstPath, .libPaths)

lib.loc <- .libPaths()[1]
library(shiny, lib.loc = lib.loc)
library(DT)
library(stringr)
library(shinydashboard)
# library(plotly)
library(plotly, lib.loc = lib.loc)
# library(png)
library(png, lib.loc = lib.loc)

library(viridisLite, lib.loc = lib.loc)
library(farver, lib.loc = lib.loc)
library(heatmaply, lib.loc = lib.loc)
# install.packages('viridis')
# library(dqshiny)
# library(rhandsontable)