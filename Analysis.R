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
source("~/Desktop/sc_methods_dev/src/functions_universal.R")
source("~/Desktop/sc_methods_dev/src/functions_specific.R")


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
tFL1 <- Seurat::Read10X( "~/Desktop/T Cells/tFL1")
tFL2 <- Seurat::Read10X( "~/Desktop/T Cells/tFL2")

# Make list of sample names
samplenames <- c(DLBCL1, DLBCL2, DLBCL3, FL1, FL2, FL3, FL4, rLN1, rLN2,  rLN3, tFL1, tFL2)
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
data$tFL1 <- calc_umap_etc(tFL1)
data$tFL2 <- calc_umap_etc(tFL2)


#Function that adds raw and smoothed gene expression to data

sample <- tFL1
samplename <- "tFL1"
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
add_gene("CD4")
add_gene("CD8B")
add_gene("CD79B")

#Sorting into categories
#Creating a list of interesting genes
genelist <- c("CD3E", "GZMK", "CD4", "CD8B", "PDCD1", "LAG3", "HAVCR2", "GATA3", "MT1A", "MT2A" )

#Use predicted lines by locfit to filter for SATB2 positive cells, including all cells
#beyond the second peak line and 90% of cells between valley and second peak line

#Creating new list
location <- list()
for (gene in genelist) {
  for( s in names(data)) {
    a <- data[[s]][[paste0("smooth_", gene)]]
    a <- a[ a < .9 ]
    location[[s]][[paste0("locmodes", gene)]] <- multimode::locmodes( a ^.15, 2 )$location^(1/.15)
    location[[s]][[paste0("thresh", gene)]] <- location[[s]][[paste0("locmodes", gene)]][2]
  }
}

#Trying different powertransformation #Trying more modes
location.2 <- list()
for (gene in genelist) {
  for( s in names(data) ) {
    a <- data[[s]][[paste0("smooth_", gene)]]
    a <- a[ a < .9 ]
    location.2[[s]][[paste0("locmodes", gene)]] <- multimode::locmodes( a ^.2, 3 )$location^(1/.2)
    location.2[[s]][[paste0("thresh", gene)]] <- location.2[[s]][[paste0("locmodes", gene)]][4]
  }
}

#Doing multimode only on CD3E pos cells
s <- samplename
location.3 <- list()
for (gene in genelist) {
  for( s in names(data)[11:12] ) {
    a <- data[[s]][[paste0("smooth_", gene)]]
    a <- a[CD3Epos[[s]]]
    a <- a[ a < .9 ]
    location.3[[s]][[paste0("locmodes", gene)]] <- multimode::locmodes( a ^.15, 2 )$location^(1/.15)
    location.3[[s]][[paste0("thresh", gene)]] <- location.3[[s]][[paste0("locmodes", gene)]][2]
  }
}

#Produce histograms
for (sample in names(data)[2]){
  for (gene in genelist[2]) {
    hist(data[[sample]][[paste0("smooth_", gene)]][data[[sample]][[paste0("smooth_", gene) ]] < .9]^.15, main=c(gene, sample, "mit B Zellen"), breaks=100)
    abline(v=location[[sample]][[paste0("locmodes", gene)]]^.15, col="yellow")
  }
}

#only FL4
for (sample in names(data)[2]){
  for (gene in genelist[2]) {
    da <- data[[sample]][[paste0("smooth_", gene)]][ CD3Epos[[ sample ]] ]
    da <- da[da < .9] 
    da <- da^.15 %>%
    hist(main=c(gene, sample), breaks=100)
    abline(v=location.3[[sample]][[paste0("locmodes", gene)]]^.15, col="yellow")
  }
}

location$rLN3$threshCD4man <- 0.17^(1/.15)
location$rLN2$threshCD4man <- 0.4^(1/.05)
location$rLN1$threshCD4man <- 0.6^(1/.05)
location$FL4$threshCD4man <- 0.55^(1/.05)
location$FL3$threshCD4man <- 0.55^(1/.05)
location$FL2$threshCD4man <- 0.575^(1/.05)
location$FL1$threshCD4man <- 0.575^(1/.05)
location$DLBCL3$threshCD4man <- 0.85^(1/.01)
location$DLBCL2$threshCD4man <- 0.85^(1/.01)
location$DLBCL1$threshCD4man <- 0.7^(1/.01)

location$FL4$threshCD8Bman <- 0.25^(1/.15)


#Do Heatmap for rLN1
# as_tibble( data$rLN1 )%>%
#  select(smooth_CD3E, smooth_PDCD1, smooth_LAG3, smooth_GATA3, smooth_HAVCR2, smooth_GZMK) %>%
#   filter(smooth_CD3E > location$FL1$threshCD3E) %>%
#   as.matrix()%>%
#   heatmap(., scale="none")

numerise_smooth <- function( gene ) {
  #Will give error when using for samples that are not FL4, no manual CD8 cutoff
  sampleDF <- data[[ sample ]][ data[[ sample ]][[ "smooth_CD3E" ]] > location[[ sample ]][[ "threshCD3E" ]], ]
  if (gene == "CD4" || gene == "CD8B") {
    return( as.numeric(sampleDF[, paste0("smooth_", gene)] > location[[sample]][[paste0("thresh", gene, "man")]]) )
  }
  else {
    return( as.numeric(sampleDF[, paste0("smooth_", gene)] > location[[sample]][[paste0("thresh", gene)]]) )
  }
}


for (sample in names(data)) {
  heatmap( sapply(genelist, numerise_smooth), scale ="none", main=sample)
}


#Filter for CD3E and GZMK pos and B cell marker CD79B
filter_gene_pos <- function(sample, gene) {
  Genepos<- data[[sample]][[paste0("smooth_", gene)]] > location[[sample]][[paste0("thresh", gene)]]
}


CD3Epos <- sapply(names(data), filter_gene_pos, "CD3E")
GZMKpos <- sapply(names(data), filter_gene_pos, "GZMK")

filter_gene_neg <- function(sample, gene) {
  Geneneg<- data[[sample]][[paste0("smooth_", gene)]] < location[[sample]][[paste0("thresh", gene)]]
}

CD79Bneg <- sapply(names(data), filter_gene_neg, "CD79B")

#Look at T cell division (CD3E pos/neg) in UMAP
lapply(names(data), function(i){
  x <- data[[i]]
  data.frame(
    sample = i,
    UMAP1 = x$UMAP1,
    UMAP2 = x$UMAP2,
    CD3Epos = CD3Epos[[i]],
    GZMKpos = GZMKpos[[i]],
    stringsAsFactors = FALSE)
}) %>% 
  bind_rows() %>%
  mutate(class=case_when(
                        CD3Epos & GZMKpos ~ "T_tox",
                        CD3Epos & !GZMKpos ~ "T_CD4",
                        TRUE ~ "other")) %>%
ggplot(aes( UMAP1, UMAP2,
            col= class))+
  geom_point( size=0.2)+
  facet_wrap(~sample)



#

