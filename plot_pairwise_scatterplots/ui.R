#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

shinyUI(fluidPage(
  titlePanel(h1("Compare E-MAP scores")),
  
  sidebarLayout( position = "left",
    
    sidebarPanel( h5("choose genes or mutants to compare"),
        helpText("Choose yeast genes (e.g. MOG1 or YJR074W) or Gsp1 mutants (e.g. T34E) to plot."),
        textInput("gene1", label = h5("first gene/ORF/mutant...")),
        textInput("gene2", label = h5("second gene/ORF/mutant...")),
        helpText("Choose which GO slims categories to higlight."),
        checkboxGroupInput("checkGO", label = h5("GO slims"),
            choices = list(
                "ribosome" = "ribosome",                         
                "transcription and mRNA processing" = "transcription and mRNA processing",
                "Golgi and ER" = "Golgi and ER",                  
                "peroxisome" = "peroxisome",                 
                "vacuole" = "vacuole",                        
                "mitochondrion" = "mitochondrion",                    
                "chromatin" = "chromatin",                   
                "cytoskeleton" = "cytoskeleton",                   
                "cell cycle" = "cell cycle",                  
                "budding" = "budding",                          
                "lipids" = "lipids",                      
                "nuclear transport" = "nuclear transport",             
                "metabolic" = "metabolic"),
            selected = 0
        )
    ),
    
    mainPanel(
      plotOutput("scatterplot")
    )
  )
))




