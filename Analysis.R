#Analysis of differences between T-cells

#Safepoint
#save.image(file="")

#Setup
library( ggplot2 )
library( Matrix )
library( locfit )
library( purrr )
library( furrr )
library( tidyverse )
library(scales)
source( "spmvar.R" )

#Load data

DLBCL1 <- Seurat::Read10X( "~/Desktop/T Cells/DLBCL1")
DLBCL2 <- Seurat::Read10X( "~/Desktop/T Cells/DLBCL2")
DLBCL3 <- Seurat::Read10X( "~/Desktop/T Cells/DLBCL3")
FL1 <- Seurat::Read10X( "~/Desktop/T Cells/FL1")
FL2 <- Seurat::Read10X( "~/Desktop/T Cells/FL2")
FL3 <- Seurat::Read10X( "~/Desktop/T Cells/FL3")
FL4 <- Seurat::Read10X( "~/Desktop/T Cells/FL4")
rLN1<- Seurat::Read10X( "~/Desktop/T Cells/rLN1")
rLN2 <- Seurat::Read10X( "~/Desktop/T Cells/rLN2")
rLN3 <- Seurat::Read10X( "~/Desktop/T Cells/rLN3")

# Make list of sample names
samplenames <- c(DLBCL1, DLBCL2, DLBCL3, FL1, FL2, FL3, FL4, rLN1, rLN2,  rLN3)
names(samplenames) <- names(data)

# Function to calculate PCA, UMAP, colSums for a sample
calc_umap_etc <- function( samplenames ) {
  
  cnts <- samplenames
  
  #Select informative genes and do PCA
  frac <- t(t(cnts) / colSums(cnts))
  gene_means <- rowMeans( frac )
  gene_vars <- rowVars_spm( frac )
  poisson_vmr <- mean( 1 / colSums( cnts ) )
  
  informative_genes <-  names(which(gene_vars / gene_means  >  1.5 * poisson_vmr ))
  pca <- irlba::prcomp_irlba( t(log1p(frac[informative_genes,]/poisson_vmr)), n = 20)$x 
  
  #Do UMAP
  umap <- uwot::umap( pca, n_neighbors = 30, min_dist = .3, metric = "cosine" )
  colnames(umap) <- c( "UMAP1", "UMAP2" )
  
  # Make data frame, add two more useful columns
  ans <- cbind( as.data.frame(pca), as.data.frame(umap) )
  
  ans$colsums <- colSums(cnts)
  
  ans  
}  


data <- list()
data$DLBCL1 <- calc_umap_etc(DLBCL1)
data$DLBCL2 <- calc_umap_etc(DLBCL2)
data$DLBCL3 <- calc_umap_etc(DLBCL3)
data$FL1 <- calc_umap_etc(FL1)
data$FL2 <- calc_umap_etc(FL2)
data$FL3 <- calc_umap_etc(FL3)
data$FL4 <- calc_umap_etc(FL4)
data$rLN1 <- calc_umap_etc(rLN1)
data$rLN2 <- calc_umap_etc(rLN2)
data$rLN3 <- calc_umap_etc(rLN3)


#Function that adds raw and smoothed gene expression to data

sample <- DLBCL1
samplename <- "DLBCL1"
gene <- "CD3E"
add_gene <- function( gene ) {
    for( samplename in names(data) ) {
      cat( samplename, " " )
      data[[ samplename ]][[ paste0( "raw_", gene ) ]] <<- samplenames[[samplename]][ gene, ]
      data[[ samplename ]][[ paste0( "smooth_", gene ) ]] <<- 
        suppressWarnings(predict( locfit.raw( 
          as.matrix( data[[ samplename ]][ ,paste0( "PC", 1:15 ) ] ), 
          data[[ samplename ]][[ paste0( "raw_", gene ) ]],
          base = log( data[[ samplename ]]$colsums ),
          alpha = .1, deg = 1, family = "poisson", ev = dat() ) ) ) 
  }
    
}

#Smooth interesting marker genes
add_gene("CD3E")
add_gene("GZMK")
add_gene("PDCD1")
add_gene("LAG3")
add_gene("HAVCR2") #alias for TIM3
add_gene("GATA3")
add_gene("MT1A")
add_gene("MT2A")

#Sorting into categories
#Creating a list of interesting genes
genelist <- c("CD3E", "GZMK", "PDCD1", "LAG3", "HAVCR2", "GATA3", "MT1A", "MT2A")
#Using the Multimode package to predict the lines in the histogram to filter pos/neg
gene <- "CD3E"

for (gene in genelist[1:2]) {
  lines <- sapply( names(data), function (s) {
    lines <- data[[s]][[ paste0("smooth_", gene)]]
    lines <- lines[lines<.9]
    multimode::locmodes(lines^.15, 2 )$location
  })
  
}



#modes_positions is a tibble including the positions of the three lines for each sample
modes_positions <- lines %>% t %>% as_tibble( rownames="sample" ) %>%
  gather( mode, pos, V1:V3 )



#Use predicted lines by locfit to filter for SATB2 positive cells, including all cells
#beyond the second peak line and 90% of cells between valley and second peak line

#Creating new list
location <- list()
for (gene in genelist) {
  for( s in names(data) ) {
    a <- data[[s]][[paste0("smooth_", gene)]]
    a <- a[ a < .9 ]
    location[[s]][[paste0("locmodes", gene)]] <- multimode::locmodes( a ^.15, 2 )$location^(1/.15)
    location[[s]][[paste0("thresh", gene)]] <- location[[s]][[paste0("locmodes", gene)]][2]
  }
}



hist(data$DLBCL2$smooth_CD3E[data$DLBCL2$smooth_CD3E< .9])
abline(v=location$DLBCL2$locmodesCD3E, col="yellow")
