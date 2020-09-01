# Find inflection points for selection of principle components
get_npcs <- function(seurat_object, create_plot = T){
  library(ggplot2)
  
  std_data = Stdev(object = seurat_object, reduction = "pca")
  ndims = length(x = std_data)
  elbow_data = data.frame(dims = 1:ndims, stdev = std_data[1:ndims])
  
  reference = elbow_data[,2]
  difference_vect = c(elbow_data[2:nrow(elbow_data),2],0)
  difference = which((reference - difference_vect) < 0.05)
  
  difference_consec = c(difference[2:length(difference)],0) - difference
  names(difference_consec) = difference
  
  npcs = as.numeric(names(which(difference_consec ==1)[1]))
  
  if(create_plot){
    
    plt = ggplot(elbow_data, aes(x = dims, y = stdev)) +
      geom_point() + 
      geom_vline(xintercept=npcs) +
      geom_text(aes(npcs+5, max(stdev), label=paste("Inferred NPCs:", npcs), vjust=0), size=6) +
      xlab("Number of principle components") +
      ylab("Standard deviation of principle components") +
      theme(text = element_text(size=17), axis.text=element_text(size=12))
    
  }
  
  out <- list(npcs=npcs, plot=plt)
  return(out)
  
}

#Sourced from Rochette Sébastien. (2016, Apr. 19). "A nice RShiny On/Off switch button". Retrieved from https://statnmap.com/2016-04-19-a-nice-rshiny-onoff-switch-button/.

switchButton <- function(inputId, label, value=FALSE, col = "GB", type="TF") {
  
  # color class
  if (col != "RG" & col != "GB") {
    stop("Please choose a color between \"RG\" (Red-Green) 
         and \"GB\" (Grey-Blue).")
  }
  if (!type %in% c("OO", "TF", "YN")){
    warning("No known text type (\"OO\", \"TF\" or \"YN\") have been specified, 
            button will be empty of text") 
  }
  if(col == "RG"){colclass <- "RedGreen"}
  if(col == "GB"){colclass <- "GreyBlue"}
  if(type == "OO"){colclass <- paste(colclass,"OnOff")}
  if(type == "TF"){colclass <- paste(colclass,"TrueFalse")}
  if(type == "YN"){colclass <- paste(colclass,"YesNo")}
  
  # No javascript button - total CSS3
  # As there is no javascript, the "checked" value implies to
  # duplicate code for giving the possibility to choose default value
  
  if(value){
    tagList(
      tags$div(class = "form-group shiny-input-container",
               tags$div(class = colclass,
                        tags$label(label, class = "control-label"),
                        tags$div(class = "onoffswitch",
                                 tags$input(type = "checkbox", name = "onoffswitch", class = "onoffswitch-checkbox",
                                            id = inputId, checked = ""
                                 ),
                                 tags$label(class = "onoffswitch-label", `for` = inputId,
                                            tags$span(class = "onoffswitch-inner"),
                                            tags$span(class = "onoffswitch-switch")
                                 )
                        )
               )
      )
    )
  } else {
    tagList(
      tags$div(class = "form-group shiny-input-container",
               tags$div(class = colclass,
                        tags$label(label, class = "control-label"),
                        tags$div(class = "onoffswitch",
                                 tags$input(type = "checkbox", name = "onoffswitch", class = "onoffswitch-checkbox",
                                            id = inputId
                                 ),
                                 tags$label(class = "onoffswitch-label", `for` = inputId,
                                            tags$span(class = "onoffswitch-inner"),
                                            tags$span(class = "onoffswitch-switch")
                                 )
                        )
               )
      )
    ) 
  }
}

app <- function(data) {
  
  load(data)
  
  shinyApp(options = list(height = 6000, width=2000), ui = shinyUI(fluidPage(includeCSS(paste(script.dir, "/qc_processing_report/www/button.css", sep="")),
                               
                               #Set style for vertical alignment
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
                               
                               # Application title
                               titlePanel(paste("QC and Data Processing Report:", sample)),
                               fluidRow(tags$hr(style="border-color: black;")),
                               #Summary data
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
                                               #wellPanel(
                                               #downloadButton("downloadData", "Download")),
                                               wellPanel(align="center", tableOutput("sum"))),
                                        column(8, align="right", id="vert", plotOutput("initQC", width="100%", height="600px"))),
                               
                               #PCA plots
                               fluidRow(tags$hr(style="border-color: black;")),
                               fluidRow(column(5, align="center", h3(tags$b("Selection of Principle Components"))), column(5, align="center", h3(tags$b("Plotting Principle Components"))), column(2, textInput("gene", "Show feature:", ""), actionButton("search", "Search"))),
                               fluidRow(column(6, align="right", id="vert", plotOutput("npcs", width="80%", height="400px")), column(6, id="vert", align="left", plotOutput("pca", width="80%", height="400px"))),
                               fluidRow(tags$hr(style="border-color: black;"))
                               
                               )),

# Define server logic required to draw a histogram
  server = shinyServer(function(input, output, session) {
  
  #Set input box
  #Set long and short versions of input for tooltip
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
  colnames(data.meta.summ) <- c("Pre-filtering", "Post-filtering")
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
  
})

)}

