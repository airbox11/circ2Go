library(shiny)
library(DT)
library(stringr)
library(shinydashboard)
library(plotly)
library(png)
# library(dqshiny)
# library(rhandsontable)
dir.base <- '~/yan150/report_work_weekly/week2020_16_online_panel/2_shiny/3_rstudio_shiny_project'
setwd(dir.base)

## 1) repare big datasets ====
source('./1_transcript_coordinate_scale/1_get_circRNA_transID_data.r')
source('./1_transcript_coordinate_scale/2_data_for_geneID.r')
datasets <- load_dataset()
source('./2_plot_in_frame/1_plot_in_frame.R')

## 2) running shiny ====
source('./3_shiny_frame/3.1_shiny_frame.r')
us <- get_shiny_modules()
shinyApp(us$ui, us$server)
# rm(list = ls())
