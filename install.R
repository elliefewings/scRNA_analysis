
libs <- c("Seurat", "dplyr", "GetoptLong", "optparse", "magrittr", "stringr", "ggplot2", "webshot", "shiny", "gridExtra", "RColorBrewer", "clustree", "tidyr")

for (i in libs) {
  if (! suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))) { 
    install.packages(i, repos = "https://ftp.fau.de/cran/")
    if (! suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))) {
      stop(paste("Unable to install package: ", i, ". Please install manually and restart.", sep=""))
    }
  }
}

if (is.null(webshot:::find_phantom())) {
  webshot::install_phantomjs()
}
