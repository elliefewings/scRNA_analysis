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
      geom_text(aes(npcs+10, max(stdev), label=paste("Inferred NPCs:", npcs), vjust=0), size=6) +
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
