#Analysis of differences between different T-cells

#Safepoint
#save.image(file="")

#Setup
library( ggplot2 )
library( Matrix )
library( locfit )
library( purrr )
library( furrr )
library( tidyverse )
source( "spmvar.R" )

#Load data
counts1 <- Seurat::Read10X( "~/Desktop/T Cells/rLN1/")
counts2 <- Seurat::Read10X( "~/Desktop/T Cells/rLN2/")
counts3 <- Seurat::Read10X( "~/Desktop/T Cells/rLN3/")
#cellinfo <- read.delim("~/sds/sd17l002/u/isabel/rawMatrix", stringsAsFactors=FALSE )

#Step1.5 Maybe figure out how to make one counts matrix from the other matrices
#Step2 Filter T-cell
#Step3 Get significant genes
#Step4 Plot histograms marking different cell conditions