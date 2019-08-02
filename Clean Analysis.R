#Analysis of differences between T-cells

#Safepoint
#save.image(file="")

##Setup
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
meta <- read.csv("~/Desktop/T Cells/sample_sheet.csv")

#Load data

samplelist <- list( DLBCL1 <- Seurat::Read10X( "~/Desktop/T Cells/DLBCL1"),
  DLBCL2 <- Seurat::Read10X( "~/Desktop/T Cells/DLBCL2"),
  DLBCL3 <- Seurat::Read10X( "~/Desktop/T Cells/DLBCL3"),
  FL1 <- Seurat::Read10X( "~/Desktop/T Cells/FL1"),
  FL2 <- Seurat::Read10X( "~/Desktop/T Cells/FL2"),
  FL3 <- Seurat::Read10X( "~/Desktop/T Cells/FL3"),
  FL4 <- Seurat::Read10X( "~/Desktop/T Cells/FL4"),
  rLN1<- Seurat::Read10X( "~/Desktop/T Cells/rLN1"),
  rLN2 <- Seurat::Read10X( "~/Desktop/T Cells/rLN2"),
  rLN3 <- Seurat::Read10X( "~/Desktop/T Cells/rLN3"),
  tFL1 <- Seurat::Read10X( "~/Desktop/T Cells/tFL1"),
  tFL2 <- Seurat::Read10X( "~/Desktop/T Cells/tFL2") )

names(samplelist) <- c("DLBCL1", "DLBCL2", "DLBCL3", "FL1", "FL2", "FL3", 
                       "FL4", "rLN1", "rLN2",  "rLN3", "tFL1", "tFL2")
genelist <- c("CD3E", "GZMK", "CD4", "CD8B", "PDCD1", "LAG3", 
              "HAVCR2", "GATA3", "MT1A", "MT2A", "CD79B", "PLAC8", 
              "PDCD1", "IL2RA" )



##PCA, UMAPs, colSums


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



##Smoothing

#Function that adds raw and smoothed gene expression to data
add_gene <- function( gene ) {
  for( samplename in names( data ) ) {
    cat( samplename, " " )
    data[[ samplename ]][[ paste0( "raw_", gene ) ]] <<- samplelist[[ samplename ]][ gene, ]
    data[[ samplename ]][[ paste0( "smooth_", gene ) ]] <<- 
      suppressWarnings(predict( locfit.raw( 
        as.matrix( data[[ samplename ]][ ,paste0( "PC", 1:15 ) ] ), 
        data[[ samplename ]][[ paste0( "raw_", gene ) ]],
        base = log( data[[ samplename ]]$colsums ),
        alpha = .1, deg = 1, family = "poisson", ev = dat() ) ) ) 
  }
}

#Smooth interesting marker genes
add_gene( "CD3E" )
add_gene( "GZMK" )
add_gene( "PDCD1" )
add_gene( "LAG3" )
add_gene( "HAVCR2" ) #alias for TIM3
add_gene( "GATA3" )
add_gene( "MT1A" )
add_gene( "MT2A" )
add_gene( "CD4" )
add_gene( "CD8B" )
add_gene( "CD79B" )
add_gene( "CD70" )
add_gene( "MAST1" )
add_gene( "ELFN1-AS1" )
add_gene( "PLAC8" ) 
add_gene( "PDCD1" )
add_gene( "IL2RA" )

#Convert data into tibble
data2 <- data %>%
  bind_rows( .id="sample" ) %>%
  as_tibble 
#Add cellnames
cellnames2 <- lapply(samplelist, colnames)
cellnames <- do.call(c, (cellnames2))
data2 <- add_column(data2, cellnames)


##Locmodes

#New List with three locmodes
location <- list()
for (gene in genelist) {
  for( s in names(data)) {
    a <- data[[s]][[paste0("smooth_", gene)]]
    a <- a[ a < .9 ]
    location[[s]][[paste0("locmodes", gene)]] <- multimode::locmodes( a ^.15, 2 )$location^(1/.15)
  }
}
#Setting manual threshholds for CD4, because locmode is not very good with CD4
location$rLN3$locmodesCD4[2] <-     0.17^( 1/.15 )
location$rLN2$locmodesCD4[2] <-     0.4^( 1/.05 )
location$rLN1$locmodesCD4[2] <-     0.6^( 1/.05 )
location$FL4$locmodesCD4[2] <-      0.55^( 1/.05 ) 
location$FL3$locmodesCD4[2] <-      0.55^( 1/.05 )
location$FL2$locmodesCD4[2] <-      0.575^( 1/.05 )
location$FL1$locmodesCD4[2] <-      0.575^( 1/.05 )
location$DLBCL3$locmodesCD4[2] <-   0.85^( 1/.01 )
location$DLBCL2$locmodesCD4[2] <-   0.85^( 1/.01 )
location$DLBCL1$locmodesCD4[2] <-   0.7^( 1/.01 )
location$FL4$locmodesCD8B[2] <-     0.25^( 1/.15 )
location$DLBCL3$locmodesCD79B[2] <- 0.01672978^( 1/.5 )
#Make list into tibble
location.tibble <- location %>%
  bind_rows( .id="sample" ) %>%
  as_tibble
#Creating a named lookup vector with the thresholds
create_thresh_vector <- function( gene) {
  ans <- c(location$DLBCL1[[paste0("locmodes", gene)]][2],
                  location$DLBCL2[[paste0("locmodes", gene)]][2],
                  location$DLBCL3[[paste0("locmodes", gene)]][2],
                  location$FL1[[paste0("locmodes", gene)]][2],
                  location$FL2[[paste0("locmodes", gene)]][2],
                  location$FL3[[paste0("locmodes", gene)]][2],
                  location$FL4[[paste0("locmodes", gene)]][2],
                  location$rLN1[[paste0("locmodes", gene)]][2],
                  location$rLN2[[paste0("locmodes", gene)]][2],
                  location$rLN3[[paste0("locmodes", gene)]][2],
                  location$tFL1[[paste0("locmodes", gene)]][2],
                  location$tFL2[[paste0("locmodes", gene)]][2])
  names(ans) <- names(samplelist)
  ans
}

CD3Ethresh <- create_thresh_vector("CD3E")
GZMKthresh <- create_thresh_vector("GZMK")
CD79Bthresh <- create_thresh_vector("CD79B")
PLAC8thresh <- create_thresh_vector("PLAC8")
PDCD1thresh <- create_thresh_vector("PDCD1")
IL2RAthresh <- create_thresh_vector("IL2RA")
  
data2 <- add_column(data2, CD3Epos = data2$smooth_CD3E > CD3Ethresh[[sample]],
                    GZMKpos = data2$smooth_GZMK > GZMKthresh[[sample]],
                    CD79Bpos = data2$smooth_CD79B > CD79Bthresh[[sample]],
                    PLAC8pos = data2$smooth_PLAC8 > PLAC8thresh[[sample]],
                    PDCD1pos = data2$smooth_PDCD1 > PDCD1thresh[[sample]],
                    IL2RApos = data2$smooth_IL2RA > IL2RAthresh[[sample]])

##Histograms

#Produce histograms


filter(data2,smooth_PLAC8 < .9 & CD3Epos & !GZMKpos & !CD79Bpos) %>%
 # filter(data2,smooth_PLAC8 < .9) %>%
    ggplot +
      geom_histogram(aes(smooth_PLAC8^.15), bins=50) +
      facet_wrap(~sample, scales="free_y")
    #geom_abline(v=location == smpl & paste0("locmodes", gene)]]^.15, col="yellow")
    #abline(v=select(filter(location.tibble, sample== smpl), location.tibble[[, paste0("locmodes", gene)]])^.15, col="yellow")
  


#Histogram of subset
for (sample in names(data)[7]){
  for (gene in genelist[2]) {
    da <- data[[sample]][[paste0("smooth_", gene)]][ CD3Epos[[ sample ]] ]
    da <- da[da < .9] 
    da <- da^.15 %>%
      hist(main=c(gene, sample), breaks=100)
    abline(v=location.3[[sample]][[paste0("locmodes", gene)]]^.15, col="yellow")
  }
}


#UMAP
data2%>%
  ggplot()+
  geom_point(aes(UMAP1, UMAP2, col=raw_PLAC8/colsums), size=.1)+
  scale_color_gradientn(colours=rev(rje::cubeHelix(10))[2:10], trans=power_trans(1/2))+
  facet_wrap(~sample)



##Plot DEG "MAST1" in CD4 T cells showing different populations
gene_overlap <- purrr::reduce(lapply(samplelist, function(i) {
  rownames(i)
}), intersect)

plot_hist_gene_celltypes <- function(gene) {
    ans <- dplyr::filter(data2, data2$CD3Epos & !data2$GZMKpos & !data2$CD79Bpos )%>%
    mutate(state=case_when(
    str_starts(ans$sample, "F" ) ~ "indolent",
    str_starts(ans$sample, "r" ) ~ "control",
    TRUE ~ "aggressive" ), 
    class=case_when(
    ans$PLAC8pos  ~ "helper"  ,
    ans$PDCD1pos ~ "follicular helper",
    ans$IL2RApos ~ "regulatory"
  )) %>%
    filter( ans[[paste0("smooth_", gene)]] < .9 )
  ggplot(ans)+
    geom_density(aes(ans[[paste0("smooth_", gene)]]^.15, col=class, group=class))+
    facet_wrap(~sample)
}

plot_hist_gene_celltypes("MAST1")
