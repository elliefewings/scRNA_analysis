#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load libraries
library(shiny)
library(ggplot2)
library(gridExtra)
library(Seurat)
library(stringr)
library(shinyBS)
library(magrittr)

# Change data load maximum
options(shiny.maxRequestSize = 500*1024^2)

# Define UI for application
ui <- shinyUI(fluidPage(
                        # Set style for vertical alignment
                        tags$head(tags$style(
                          HTML('
                               #vert { 
                               display: flex;
                               align-items: center;
                               margin-top: 50px;
                               }
                               .tooltip .tooltip-inner {
                               max-width: 100%;
                               }
                               '))),
                        
                        # Create Input
                        fileInput("file", 'Choose Rdata to upload',
                                  accept = c('.Rdata')),
                        actionButton("load", "Generate report"),
                        bsTooltip("file", title="Select Rdata generated from pipeline, 's02_cluster_identity.Rdata'", trigger="hover", placement = "bottom"),
                        
                        # Application title
                        titlePanel(textOutput("title")),
                        fluidRow(tags$hr(style="border-color: black;")),
                        
                        # Summary data
                        fluidRow(column(3, align="left", h3(tags$b("Script options"))), column(7, align="center", h3(tags$b("Heatmap of Differential Features"))), 
                                 column(2, downloadButton("download", "Download Differential Features")
                                 )),
                        fluidRow(column(3, 
                                        wellPanel(
                                          uiOutput("input"),
                                          uiOutput("inlong"),
                                          uiOutput("output"),
                                          uiOutput("outlong"),
                                          uiOutput("npc"),
                                          uiOutput("res")
                                        )),
                                 column(9, align="right", id="vert", plotOutput("heatmap", width="100%", height="900px"))),
                        
                        # PCA plots
                        fluidRow(tags$hr(style="border-color: black;")),
                        fluidRow(column(3, align="center", h3(tags$b("Cell Type Markers"))), column(7, align="center", h3(tags$b("UMAP Clustering"))), column(2, textInput("gene", "Show feature:", ""), actionButton("search", "Search"))),
                        fluidRow(column(3, align="center", offset=1, id="vert", wellPanel(tableOutput("markers"))), column(8, id="vert", align="center", plotOutput("umap", width="80%", height="700px"))),
                        fluidRow(tags$hr(style="border-color: black;"))
                        
                        ))

# Define server logic
server <- shinyServer(function(input, output, session) {
  
  # Load data
  load_Rdata <- function(){
    if(is.null(input$file)){return(NULL)} 
    rdata <- isolate({input$file})
    load(rdata$datapath, envir = .GlobalEnv)
  }
  
  # Create event when report load button is activated
  observeEvent(input$load,{
    load_Rdata()
    
    #Set text outputs
    output$title <- renderText({paste("Clustering and Cell Identity Report:", sample)})
    
    inshort <- ifelse(nchar(opt2$input) > 50, 
                      substr(opt2$input, nchar(opt2$input)-50, nchar(opt2$input)),
                      opt2$input) %>% sub(".*?/", "/", .)
    
    output$input <- renderText({HTML(paste("<b>","Input File:", "</b>", "...", inshort))})
    output$inlong <- renderUI({
      bsTooltip("input", title=opt2$input, trigger="hover", placement = "bottom")
    })
    
    #Set long and short versions of output for tooltip
    outshort <- ifelse(nchar(opt2$output) > 50, 
                       substr(opt2$output, nchar(opt2$output)-50, nchar(opt2$output)),
                       opt2$output) %>% sub(".*?/", "/", .)
    output$output <- renderText({HTML(paste("<b>", "Output Directory:", "</b>", "...", outshort))})
    output$outlong <- renderUI({
      bsTooltip("output", title=opt2$output, trigger="hover", placement = "bottom")
    })
    output$npc <- renderText({HTML(paste("<b>", "Number of Principle Components:", "</b>", opt2$npc))})
    output$res <- renderText({HTML(paste("<b>", "Resolution for Clustering:", "</b>", opt2$res))})
    
    #Create action for download button
    output$download <- downloadHandler(
      filename = function() {
        paste(sample, "_cluster_differential_features.csv", sep = "")
      },
      content = function(file) {
        write.csv(all.markers, file, row.names = FALSE)
      }
    )
    
    #Set heatmap
    output$heatmap <- renderPlot(heat)
    
    #Set markers table
    #colnames(ids) <- c("")
    
    #Set PCA and feature plot event
    if (opt2$markers != "") { 
      
      #Set markers table
      output$markers <- renderTable(ids)
      output$umap <- renderPlot(pca2 + theme(text = element_text(size=20)))
      
      featurePlotter <- eventReactive(input$search, {
        if (input$gene == "") {
          plot.g <- pca2 + theme(text = element_text(size=20))
          
        }
        if (input$gene != "") {
          genen <- str_to_title(input$gene) %>% str_replace("MT-", mt.patt)
          plot.g <- FeaturePlot(data, features = genen) + theme(text = element_text(size=20))
        }
        plot.g
      })
      
      observeEvent(input$search, {
        output$umap <- renderPlot({
          plot.g <- featurePlotter()
          plot.g
        })
      })
    } else {
      
      #Set markers table
      output$markers <- renderTable(
        validate(
        need(exists("ids"), "No markers supplied")
      ))
      
      output$umap <- renderPlot(pca + theme(text = element_text(size=20)))
      
      featurePlotter <- eventReactive(input$search, {
        if (input$gene == "") {
          plot.g <- pca + theme(text = element_text(size=20))
          
        }
        if (input$gene != "") {
          genen <- str_to_title(input$gene) %>% str_replace("MT-", mt.patt)
          plot.g <- FeaturePlot(data, features = genen) + theme(text = element_text(size=20))
        }
        plot.g
      })
      
      observeEvent(input$search, {
        output$umap <- renderPlot({
          plot.g <- featurePlotter()
          plot.g
        })
      })
      
      
    }
    
    
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)

