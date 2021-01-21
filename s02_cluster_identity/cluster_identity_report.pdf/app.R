#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

load(rdata)

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
         '))),
  
  # Application title
  titlePanel(textOutput("title")),
  fluidRow(HTML("<h3>For interactive report, load RData output here: <a href='https://saezlab.shinyapps.io/cluster_identity_report'> </a></h3>")),
  fluidRow(tags$hr(style="border-color: black;")),
  
  # Summary data
  fluidRow(column(3, align="left", h3(tags$b("Script options"))), column(9, align="center", h3(tags$b("Heatmap of Differential Features")))
           ),
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
  fluidRow(column(3, align="center", h3(tags$b(textOutput("headmarkers")))), column(8, align="center", h3(tags$b("UMAP Clustering")))),
  fluidRow(column(3, align="center", offset=1, id="vert", wellPanel(tableOutput("markers"))), column(8, id="vert", align="center", plotOutput("umap", width="80%", height="700px"))),
  fluidRow(tags$hr(style="border-color: black;"))
    ))

# Define server logic
server <- shinyServer(function(input, output, session) {
  
    #Set text outputs
    output$title <- renderText({paste("Clustering and Cell Identity Report:", sample)})
    
    inshort <- ifelse(nchar(opt2$input) > 50, 
                      substr(opt2$input, nchar(opt2$input)-50, nchar(opt2$input)),
                      opt2$input) %>% sub(".*?/", "/", .)
    
    output$input <- renderText({HTML(paste("<b>","Input File:", "</b>", "...", inshort))})

    #Set long and short versions of output for tooltip
    outshort <- ifelse(nchar(opt2$output) > 50, 
                       substr(opt2$output, nchar(opt2$output)-50, nchar(opt2$output)),
                       opt2$output) %>% sub(".*?/", "/", .)
    output$output <- renderText({HTML(paste("<b>", "Output Directory:", "</b>", "...", outshort))})
    output$npc <- renderText({HTML(paste("<b>", "Number of Principle Components:", "</b>", opt2$npc))})
    output$res <- renderText({HTML(paste("<b>", "Resolution for Clustering:", "</b>", opt2$res))})
    
    #Set heatmap
    output$heatmap <- renderPlot(heat)

    #Render id table label
    output$headmarkers <- renderText("Cell Type Markers")
    
    #Create shortened ID table if generated from Panglao
    if (opt2$markers == "") {
      ids <- ids[1:4]
      output$headmarkers <- renderText("Cell Type Markers (predicted by PanglaoDB)")
    }
    
    output$markers <- renderTable(ids)
    
    #Render labelled clusters       
    output$umap <- renderPlot(pca2 + theme(text = element_text(size=20)))
    
  })


