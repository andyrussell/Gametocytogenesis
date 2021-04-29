#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Sorted gametocyte data"),
  
  
  fluidRow(column(width=9,plotOutput("topplot",click="topplot_click")),column(width=3, textInput("gene","Search"),
                                                                              checkboxInput("slow","Include slows",value=T), checkboxInput("highlight","Circle early responders of interest",value=F))
    ),
  h3(textOutput("gene")),textOutput("geneproduct"),
    # Show a plot of the generated distribution
  fluidRow(column(width=6,
                 # h1(textOutput("gene_name"),
                 
                     plotOutput("basicplot")),column(width=6,
                                                     # h1(textOutput("gene_name"),
                                                     plotOutput("mergedeffectplot")))
  )
      
      
  )

