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

# Load extra functions
source("./source.R")
options(shiny.maxRequestSize = 200*1024^2)

# Define UI for application
ui <- shinyUI(fluidPage(theme= "button.css",
                        
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
                               ')),
                          tags$link(rel = "stylesheet", type = "text/css", href = "button.css")),
                        
                        # Create Input
                        fileInput("file", 'Choose Rdata to upload',
                                  accept = c('.Rdata')),
                        actionButton("load", "Generate report"),
                        bsTooltip("file", title="Select Rdata generated from pipeline, 's01_qc_processing.Rdata'", trigger="hover", placement = "bottom"),
                        
                        # Application title
                        titlePanel(textOutput("title")),
                        fluidRow(tags$hr(style="border-color: black;")),
                        
                        # Summary data
                        fluidRow(column(4, align="left", h3(tags$b("Script options"))), column(6, align="center", h3(tags$b("Quality Control Plots"))), 
                                 column(2, switchButton(inputId="filter", label = "Show Filtered Data", value = FALSE, col = "GB", type = "TF")
                                 )),
                        fluidRow(column(4, 
                                        wellPanel(
                                          uiOutput("input"),
                                          uiOutput("inlong"),
                                          uiOutput("output"),
                                          uiOutput("outlong"),
                                          uiOutput("mincells"),
                                          uiOutput("minfeatures"),
                                          uiOutput("maxfeatures"),
                                          uiOutput("maxpercentmt")
                                        ),
                                        wellPanel(align="center", tableOutput("sum"))),
                                 column(8, align="right", id="vert", plotOutput("initQC", width="100%", height="600px"))),
                        
                        # PCA plots
                        fluidRow(tags$hr(style="border-color: black;")),
                        fluidRow(column(5, align="center", h3(tags$b("Selection of Principle Components"))), column(5, align="center", h3(tags$b("Plotting Principle Components"))), column(2, textInput("gene", "Show feature:", ""), actionButton("search", "Search"))),
                        fluidRow(column(6, align="right", id="vert", plotOutput("npcs", width="80%", height="400px")), column(6, id="vert", align="left", plotOutput("pca", width="80%", height="400px"))),
                        fluidRow(tags$hr(style="border-color: black;")),
                        
                        # Demultiplex plots
                        fluidRow(column(5, align="center", h3(tags$b(textOutput("doubtitle")))), column(7, align="center", h3(tags$b(textOutput("hashtitle"))))),
                        fluidRow(column(5, align="right", id="vert", plotOutput("doublets", width="80%", height="700px")), column(7, id="vert", align="right", plotOutput("hashtags", width="80%", height="700px")))
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
    output$title <- renderText({paste("QC and Data Processing Report:", sample)})
    
    inshort <- ifelse(nchar(opt$input) > 50, 
                      substr(opt$input, nchar(opt$input)-50, nchar(opt$input)),
                      opt$input) %>% sub(".*?/", "/", .)
    
    output$input <- renderText({HTML(paste("<b>","Input Directory:", "</b>", "...", inshort))})
    output$inlong <- renderUI({
      bsTooltip("input", title=opt$input, trigger="hover", placement = "bottom")
    })
    
    #Set long and short versions of output for tooltip
    outshort <- ifelse(nchar(opt$output) > 50, 
                       substr(opt$output, nchar(opt$output)-50, nchar(opt$output)),
                       opt$output) %>% sub(".*?/", "/", .)
    output$output <- renderText({HTML(paste("<b>", "Output Directory:", "</b>", "...", outshort))})
    output$outlong <- renderUI({
      bsTooltip("output", title=opt$output, trigger="hover", placement = "bottom")
    })
    output$mincells <- renderText({HTML(paste("<b>", "Minimum Cells Filter:", "</b>", opt$mincells))})
    output$minfeatures <- renderText({HTML(paste("<b>", "Minimum Features Filter:", "</b>", opt$minfeatures))})
    output$maxfeatures <- renderText({HTML(paste("<b>", "Maximum Features Filter:", "</b>", opt$maxfeatures))})
    output$maxpercentmt <- renderText({HTML(paste("<b>", "Maximum Percentage MT Filter:", "</b>", opt$maxpercentmt))})
    
    #Set QC plots 
    plotselection <- reactive({ 
      if (input$filter == FALSE) {
        qc1.1 <- qc1[[1]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
        qc1.2 <- qc1[[2]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
        qc1.3 <- qc1[[3]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
        qc2 <- qc2 + theme(legend.position = "none") + labs(title = "Percentage mitochondrial features vs number of molecules detected")
        qc3 <- qc3 + theme(legend.position = "none") + labs(title = "Number of unique features vs number of molecules detected")
        plot <- grid.arrange(arrangeGrob(qc1.1, qc1.2, qc1.3, ncol=3), arrangeGrob(qc2, qc3, ncol=2), heights=c(2.5/4, 1.5/4), ncol=1)
      } else {
        qc1.1 <- qc1.f[[1]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
        qc1.2 <- qc1.f[[2]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
        qc1.3 <- qc1.f[[3]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
        qc2 <- qc2.f + theme(legend.position = "none") + labs(title = "Percentage mitochondrial features vs number of molecules detected")
        qc3 <- qc3.f + theme(legend.position = "none") + labs(title = "Number of unique features vs number of molecules detected")
        plot <- grid.arrange(arrangeGrob(qc1.1, qc1.2, qc1.3, ncol=3), arrangeGrob(qc2, qc3, ncol=2), heights=c(2.5/4, 1.5/4), ncol=1)
      }
      plot
    })
    
    output$initQC <- renderPlot({plotselection()})
    
    #Set pre/post filtering table
    row.names(data.meta.summ) <- c(
      "<b>Number of cells</b>",
      "<b>Median nCount_RNA</b>",
      "<b>Minimum nCount_RNA</b>",
      "<b>Maximum nCount_RNA</b>",
      "<b>Median nFeature_RNA</b>",
      "<b>Minimum nFeature_RNA</b>",
      "<b>Maximum nFeature_RNA</b>",
      "<b>Median percent.mt</b>",
      "<b>Minimum percent.mt</b>",
      "<b>Maximum percent.mt</b>"
    )
    colnames(data.meta.summ)[1:2] <- c("Pre-filtering", "Post-filtering")
    output$sum <- renderTable(data.meta.summ, spacing = "l", rownames=TRUE, digits=0, hover=TRUE, sanitize.text.function=function(x){x})
    
    #Set npcs plot
    output$npcs <- renderPlot(npcs$plot + theme(text = element_text(size=20)))
    
    #Set PCA and feature plot event
    output$pca <- renderPlot(pca + theme(text = element_text(size=20), legend.position = "none"))
    
    featurePlotter <- eventReactive(input$search, {
      if (input$gene == "") {
        plot4 <- pca + theme(text = element_text(size=20), legend.position = "none")
        
      }
      if (input$gene != "") {
        genen <- str_to_title(input$gene) %>% str_replace("MT-", mt.patt)
        plot4 <- FeaturePlot(data, features = genen) + theme(text = element_text(size=20), legend.position = "none")
      }
      plot4
    })
    
    observeEvent(input$search, {
      output$pca <- renderPlot({
        plot4 <- featurePlotter()
        plot4
      })
    })
    
    #Plot hashtag data
    
    if (!is.null(opt$hashtag)) {
    output$doubtitle <- renderText("Number of Doublets Identified")
    output$hashtitle <- renderText("Expression Counts Over Hashtags")
    output$doublets <- renderPlot(doublet)
    output$hashtags <- renderPlot(ridge)
    } 
    
    
    
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)

