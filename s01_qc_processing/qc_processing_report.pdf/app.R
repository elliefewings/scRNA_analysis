#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#load(rdata)

# Define UI for application
ui <- shinyUI(fluidPage(
  
  #Set style for vertical alignment
  tags$head(tags$style(
    HTML('
         #vert { 
            display: flex;
            align-items: center;
            margin-top: 50px;
         }
         '))),
  
   # Application title
   titlePanel(htmlOutput("title")),
   fluidRow(HTML("<h3>For interactive report, load RData output here: <a href='https://saezlab.shinyapps.io/qc_processing_report'> </a></h3>")),
   fluidRow(tags$hr(style="border-color: black;")),
   #Summary data
   fluidRow(column(4, align="left", h3(tags$b("Script options"))), column(8, align="center", h3(tags$b("Quality Control Plots")))),
   fluidRow(column(4, 
                   wellPanel(
                     uiOutput("input"),
                     uiOutput("output"),
                     uiOutput("mincells"),
                     uiOutput("minfeatures"),
                     uiOutput("maxfeatures"),
                     uiOutput("maxpercentmt")
                     ),
                   wellPanel(align="center", tableOutput("sum"))),
            column(8, align="right", id="vert", plotOutput("initQC", width="100%", height="600px"))),
  
   #PCA plots
   fluidRow(tags$hr(style="border-color: black;")),
   fluidRow(column(6, align="center", h3(tags$b("Selection of Principle Components"))), column(6, align="center", h3(tags$b("Plotting Principle Components")))),
   fluidRow(column(6, align="center", id="vert", plotOutput("npcs", width="80%", height="400px")), column(6, id="vert", align="center", plotOutput("pca", width="80%", height="400px"))),
   fluidRow(tags$hr(style="border-color: black;")),
  
    # Clustree plot
    fluidRow(column(12, align="center",  h3(tags$b("Selection of K from Clustree")))),
    fluidRow(column(12, align="center", id="vert", plotOutput("clust", width="60%", height="800px"))),
    fluidRow(tags$hr(style="border-color: black;")),
  
    #Demultiplex plots
    fluidRow(column(5, align="center", h3(tags$b(textOutput("doubtitle")))), column(7, align="center", h3(tags$b(textOutput("hashtitle"))))),
    fluidRow(column(5, align="right", id="vert", plotOutput("doublets", width="80%", height="700px")), column(7, id="vert", align="right", plotOutput("hashtags", width="80%", height="700px")))
      
   ))

# Define server logic
server <- shinyServer(function(input, output, session) {
   
  #Set text outputs
  #Set title based on data quality
  headtitle <- NULL
  
  #Set title and colour if > 1000 cells
  headtitle[nsinglets > 1000] <- paste("QC Report:", sample)
  
  #Set title and colour if between 500 and 1000 cells
  headtitle[nsinglets > 500 & nsinglets <= 1000] <- paste("QC Report: ", sample, ' <font style=color:orange !important >(WARNING: Fewer than 1000 singlets)</font>', sep="")
  
  #Set title and colour if < 500 cells
  headtitle[nsinglets <= 500] <- paste("QC Report: ", sample, ' <font color="red">(WARNING: Fewer than 500 singlets)</font>', sep="")
  
  output$title <- renderText({HTML(headtitle)})

  inshort <- ifelse(nchar(opt$input) > 50, 
                    substr(opt$input, nchar(opt$input)-50, nchar(opt$input)),
                    opt$input) %>% sub(".*?/", "/", .)
  
  output$input <- renderText({HTML(paste("<b>","Input Directory:", "</b>", "...", inshort))})
  
  outshort <- ifelse(nchar(opt$output) > 50, 
                     substr(opt$output, nchar(opt$output)-50, nchar(opt$output)),
                     opt$output) %>% sub(".*?/", "/", .)
  output$output <- renderText({HTML(paste("<b>", "Output Directory:", "</b>", "...", outshort))})
  output$mincells <- renderText({HTML(paste("<b>", "Minimum Cells Filter:", "</b>", opt$mincells))})
  output$minfeatures <- renderText({HTML(paste("<b>", "Minimum Features Filter:", "</b>", opt$minfeatures))})
  output$maxfeatures <- renderText({HTML(paste("<b>", "Maximum Features Filter:", "</b>", opt$maxfeatures))})
  output$maxpercentmt <- renderText({HTML(paste("<b>", "Maximum Percentage MT Filter:", "</b>", opt$maxpercentmt))})
  
  #Set QC plots 
  output$initQC <- renderPlot({
  qc1.1 <- qc1[[1]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  qc1.2 <- qc1[[2]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  qc1.3 <- qc1[[3]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  qc2 <- qc2 + theme(legend.position = "none", title = element_text(size=10)) + labs(title = "Percentage mitochondrial features vs number of molecules detected")
  qc3 <- qc3 + theme(legend.position = "none", title = element_text(size=10)) + labs(title = "Number of unique features vs number of molecules detected")
  plot <- grid.arrange(arrangeGrob(qc1.1, qc1.2, qc1.3, ncol=3), arrangeGrob(qc2, qc3, ncol=2), heights=c(2.5/4, 1.5/4), ncol=1)
  plot})
  
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
  
  # Plot clustree data
  output$clust <- renderPlot(clust)
  #Plot hashtag data
  if (!is.null(opt$hashtag)) {
    output$doubtitle <- renderText("Number of Doublets Identified")
    output$hashtitle <- renderText("Expression Counts Over Hashtags")
    output$doublets <- renderPlot(doublet)
    output$hashtags <- renderPlot(ridge)
  } 
  
})

