repo <- "http://cran.us.r-project.org"
install.packages('devtools', repos = repo)
require(devtools)
install_version("ggplot2", version = "1.0.1", repos = repo)
install_version("RSQLite", version = "1.0.0", repos = repo)
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5") # too bad I cannot install specified version for rhdf5
