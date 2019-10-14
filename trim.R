# function to trim the gene set data
trim<- function (x) gsub("^\\s+|\\s+$", "", x)