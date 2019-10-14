# Function to reads ina gene expression dataset in GCT format and converts it into an R data frame
GSEA.Gct2Frame <- function(filename = "NULL") { 
  ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
   ds <- ds[-1]
   return(ds)
}
