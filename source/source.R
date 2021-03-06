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


##Assign identity function
assign.identity <- function(seurat_object, markers){
  
  #Find significant marker genes for each cluster
  all.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  all.markers <- all.markers[!grepl(",", all.markers$cluster),]
  
  #Create one row for each cell.type label and marker gene
  identity <- strsplit(markers$Markers, ",")
  markers <- data.frame(Cell.Type = rep(markers$Cell.Type, sapply(identity, length)), Markers = unlist(identity))
  
  #Make sure both gene names are in sentence case
  markers$Markers <- str_to_sentence(markers$Markers)
  all.markers$gene <- str_to_sentence(all.markers$gene)
  
  #Merge identity on if available
  add <- merge(markers, all.markers, by.x="Markers", by.y="gene", all.y=TRUE)
  
  proportionID <- add %>% group_by(cluster, Cell.Type) %>% mutate(count=length(Cell.Type)) %>% group_by(cluster) %>% mutate(prop=count/length(cluster))
  
  #Subset identity and cluster numbers
  ids <- proportionID %>% select("Cell.Type", "cluster", "prop") %>% unique()
  
  #Sort by proportion of genes in marker list
  ids <- ids[order(ids$cluster, ids$prop),]
  
  #Add formatted column
  ids$name <- paste(ids$Cell.Type, "(", round(ids$prop, 4), ")", sep="")
  
  #Add label column
  ids <- ids %>% group_by(cluster) %>% summarise(label = toString(rev(name))) %>% ungroup()
  
  #Remove "NA," values 
  ids$label <- gsub("NA*\\(.*?\\),", "", ids$label)
  
  #Convert "NA" character to NA
  ids$label[grepl("NA", ids$label)] <- NA
  
  ids <- separate(ids, col = label, into = c("label1", "label2", "label3"), sep=",")
  
  #Reorder
  ids <- ids[order(ids$cluster),]
  
  #Return
  return(ids)
}
  
# Assign cell identities using Panglao
assign.pldb <- function(seurat_object, pldb){
  
  #Find significant marker genes for each cluster
  all.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  
  #Reformat
  markers <- data.frame(Cell.Type = pldb$cell.type, Markers = pldb$official.gene.symbol)
  
  #Make sure both gene names are in sentence case
  markers$Markers <- str_to_sentence(markers$Markers)
  all.markers$gene <- str_to_sentence(all.markers$gene)
  
  #Merge identity on if available
  add <- merge(markers, all.markers, by.x="Markers", by.y="gene", all.y=TRUE)
  
  proportionID <- add %>% group_by(cluster, Cell.Type) %>% mutate(count=length(Cell.Type)) %>% group_by(cluster) %>% mutate(prop=count/length(cluster))
  
  #Subset identity and cluster numbers
  ids <- proportionID %>% select("Cell.Type", "cluster", "prop") %>% unique()
  
  #Sort by proportion of genes in marker list
  ids <- ids[order(ids$cluster, ids$prop),]
  
  #Gather identities if multiple per cluster
  ids <- ids %>% group_by(cluster) %>% mutate(label=paste(rev(paste(Cell.Type, "(", round(prop, 4), ")", sep=""))[1], 
                                                          rev(paste(Cell.Type, "(", round(prop, 4), ")", sep=""))[2], 
                                                          rev(paste(Cell.Type, "(", round(prop, 4), ")", sep=""))[3],
                                                          rev(paste(Cell.Type, "(", round(prop, 4), ")", sep=""))[4], sep=",")) %>% select(-Cell.Type, -prop) %>% unique()
  
  #Remove "NA," values 
  ids$label <- gsub("NA*\\(.*?\\),", "", ids$label)
  
  #Convert "NA" character to NA
  ids$label[grepl("NA\\(", ids$label)] <- NA
  
  ids <- separate(ids, col = label, into = c("label1", "label2", "label3"), sep=",")
  
  #Reorder
  ids <- ids[order(ids$cluster),]
  
  #Return
  return(ids)
}
