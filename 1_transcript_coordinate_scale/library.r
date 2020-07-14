# lib.loc <- '~/R/x86_64-pc-linux-gnu-library/4.0'
# lib.loc <- '~/R/x86_64-redhat-linux-gnu-library/3.6'
lib.loc <- .libPaths()[1]
library(shiny)
library(DT)
library(stringr)
library(shinydashboard)
# library(plotly)
library(plotly, lib.loc = lib.loc)
# library(png)
library(png, lib.loc = lib.loc)

library(viridisLite)
library(heatmaply)
# install.packages('viridis')
# library(dqshiny)
# library(rhandsontable)