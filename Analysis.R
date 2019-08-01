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
meta <- read.csv("~/Desktop/T Cells/sample_sheet.csv")


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

# Make list of raw counts
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
add_gene("CD70")
add_gene("MAST1")
add_gene("ELFN1-AS1")
add_gene("PLAC8")
add_gene("PDCD1")
add_gene("IL2RA")




#Creating a list of interesting genes
genelist <- c("CD3E", "GZMK", "CD4", "CD8B", "PDCD1", "LAG3", 
              "HAVCR2", "GATA3", "MT1A", "MT2A", "CD79B", "PLAC8", 
              "PDCD1", "IL2RA" )

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

# #Trying different powertransformation 
# #Trying more modes
# location.2 <- list()
# for (gene in genelist) {
#   for( s in names(data) ) {
#     a <- data[[s]][[paste0("smooth_", gene)]]
#     a <- a[ a < .9 ]
#     location.2[[s]][[paste0("locmodes", gene)]] <- multimode::locmodes( a ^.2, 3 )$location^(1/.2)
#     location.2[[s]][[paste0("thresh", gene)]] <- location.2[[s]][[paste0("locmodes", gene)]][4]
#   }
# }

# #Doing multimode only on CD3E pos cells
# s <- samplename
# location.3 <- list()
# for (gene in genelist) {
#   for( s in names(data)[11:12] ) {
#     a <- data[[s]][[paste0("smooth_", gene)]]
#     a <- a[CD3Epos[[s]]]
#     a <- a[ a < .9 ]
#     location.3[[s]][[paste0("locmodes", gene)]] <- multimode::locmodes( a ^.15, 2 )$location^(1/.15)
#     location.3[[s]][[paste0("thresh", gene)]] <- location.3[[s]][[paste0("locmodes", gene)]][2]
#   }
# }

#Produce histograms
for (sample in names(data)[7]){
  for (gene in genelist[4]) {
    hist(data[[sample]][[paste0("smooth_", gene)]][data[[sample]][[paste0("smooth_", gene) ]] < .9]^.15, main=c(gene, sample, "with B Cells"), breaks=100)
    abline(v=location[[sample]][[paste0("locmodes", gene)]]^.15, col="yellow")
  }
}

#only FL4
for (sample in names(data)[7]){
  for (gene in genelist[2]) {
    da <- data[[sample]][[paste0("smooth_", gene)]][ CD3Epos[[ sample ]] ]
    da <- da[da < .9] 
    da <- da^.15 %>%
    hist(main=c(gene, sample), breaks=100)
    abline(v=location.3[[sample]][[paste0("locmodes", gene)]]^.15, col="yellow")
  }
}

#Setting manual threshholds for CD4, because locmode is not very good with CD4
location$rLN3$threshCD4 <- 0.17^(1/.15)
location$rLN2$threshCD4 <- 0.4^(1/.05)
location$rLN1$threshCD4 <- 0.6^(1/.05)
location$FL4$threshCD4 <- 0.55^(1/.05)
location$FL3$threshCD4 <- 0.55^(1/.05)
location$FL2$threshCD4 <- 0.575^(1/.05)
location$FL1$threshCD4 <- 0.575^(1/.05)
location$DLBCL3$threshCD4 <- 0.85^(1/.01)
location$DLBCL2$threshCD4 <- 0.85^(1/.01)
location$DLBCL1$threshCD4 <- 0.7^(1/.01)

location$FL4$threshCD8Bman <- 0.25^(1/.15)




#Heatmaps



#Do Heatmap for rLN1
# as_tibble( data$rLN1 )%>%
#  select(smooth_CD3E, smooth_PDCD1, smooth_LAG3, smooth_GATA3, smooth_HAVCR2, smooth_GZMK) %>%
#   filter(smooth_CD3E > location$FL1$threshCD3E) %>%
#   as.matrix()%>%
#   heatmap(., scale="none")

#Makes a numeric vector from smoothed gene expression that lies beyond threshhold
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

#Produces heatmaps
for (sample in names(data)) {
  heatmap( sapply(genelist, numerise_smooth), scale ="none", main=sample)
}



#New colorscale for heatmap to improve distinction
x <- matrix()
gene <- "CD3E"
for (sample in names(data)) {
  for (gene in genelist[1:4]) {
    x[[sample]][[ gene]] <- data[[sample]][[paste0("smooth_", gene )]]^.15 %>%
     scale( center=location[[ sample ]][[ paste0( "thresh", gene ) ]], 
             scale=location[[ sample ]][[ paste0( "locmodes", gene ) ]][3] - 
               location[[ sample ]][[ paste0( "locmodes", gene ) ]][2] )
}
  }
 %>%
  do.call(rbind, .)#%>%
  matrix()#%>%
  bind_rows() #%>%

scale_manual <- function(sample, gene) {
  x <- data[[ sample ]][[ paste0( "smooth_", gene ) ]]
  xi <- location[[ sample ]][[ paste0( "thresh", gene ) ]]
  h <- location[[ sample ]][[ paste0( "locmodes", gene ) ]][ 3 ] - 
    location[[ sample ]][[ paste0( "thresh", gene ) ]]
  y= 1/ (1+ exp((x+xi)/h))
}

scale2_manual <- function(sample, gene, h) {
  x <- data[[ sample ]][[ paste0( "smooth_", gene ) ]]
  xi <- location[[ sample ]][[ paste0( "thresh", gene ) ]]
  y= 1/ (1+ exp((x+xi)/h))
  return(y)
}

#New column with cellid
genes_scaled <- cbind(  CD3E = unlist(sapply( names( data ), scale_manual, "CD3E" ) ), 
                        CD4 = unlist(sapply( names( data ), scale_manual, "CD4" ) ),
                        CD8B = unlist(sapply( names( data ), scale_manual, "CD8B" ) ),
                        GZMK = unlist(sapply( names( data ), scale_manual, "GZMK" ) ) )
cellid <- c(1:length(data$DLBCL1$smooth_CD3E), 1:length(data$DLBCL2$smooth_CD3E), 
            1:length(data$DLBCL3$smooth_CD3E), 1:length(data$FL1$smooth_CD3E), 
            1:length(data$FL2$smooth_CD3E), 1:length(data$FL3$smooth_CD3E),
            1:length(data$FL4$smooth_CD3E), 1:length(data$rLN1$smooth_CD3E),
            1:length(data$rLN2$smooth_CD3E), 1:length(data$rLN3$smooth_CD3E),
            1:length(data$tFL1$smooth_CD3E), 1:length(data$tFL2$smooth_CD3E))
genes_scaled <- cbind(genes_scaled, cellid=cellid)
genes_scaled <- gather(as.data.frame(genes_scaled), key=gene, value=intensity, -cellid)
 
#For just one sample, DLBCL2
genes_scaled <- cbind(  CD3E = scale_manual( "DLBCL2", "CD3E" ), 
                        CD4 = scale2_manual( "DLBCL2", "CD4", 2 ) ,
                        CD8B = scale_manual( "DLBCL2", "CD8B" ) ,
                        GZMK = scale_manual( "DLBCL2", "GZMK" ) )

pheatmap(genes_scaled[order(genes_scaled[, "CD4"]), ], 
         cluster_rows = FALSE)


cellid <- c(1:length(data$DLBCL2$smooth_CD3E))
genes_scaled <- cbind(genes_scaled, cellid=cellid) 

genes_scaled <- genes_scaled %>%
as.data.frame() %>%
  arrange(CD4) %>%
  rownames_to_column("id")%>%
as.matrix() 
  pheatmap(genes_scaled)
  
  
sample <- "DLBCL2"
gene <- "CD4"
  
x <- data[[ sample ]][[ paste0( "smooth_", gene ) ]]
#x <- x[x<.9]
x <- x^.1
#x <- seq(-5, 5, .1)
plot(x, 1/(1+exp(-(x-.25)/.025)))
hist(x)  
  
gene <- "CD8B"
z <- data[[ sample ]][[ paste0( "smooth_", gene ) ]]
#z <- z[z<.9]
z <- z^.1
plot(z, 1/(1+exp(-(z-.25)/.025)))
  #ggplot( aes(gene, id)) +
  #geom_tile(aes(fill=intensity))

z <- 1/(1+exp(-(z-.25)/.025))
x <- 1/(1+exp(-(x-.25)/.025))

gene <- "CD3E"
w <- data[[ sample ]][[ paste0( "smooth_", gene ) ]]
#z <- z[z<.9]
w <- w^.1
w <- 1/(1+exp(-(w-.25)/.025))
plot(w, 1/(1+exp(-(w-.35)/.025)))

testmatrix <- cbind(CD8=z, CD4=x, CD3E=w)
testmatrix <- testmatrix[CD3Epos$DLBCL2]
pheatmap(testmatrix)




heatmap(genes_scaled[4123:9461, ], main ="DLBCL2")

DLBCL2smooth <- as_tibble(data$DLBCL2) %>%
  select(smooth_CD3E, smooth_CD8B, smooth_GZMK, smooth_CD4)

ggplot(DLBCL2smooth) + geom_tile(aes(fill = genes_scaled[4123:9461, ]))



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
#Manually setting threshold to fix DLBCL3 sample, threshold found by Multimode wiht power transformation .5
CD79Bneg$DLBCL3 <-  data$DLBCL3$smooth_CD79B < 0.01672978^(1/.5)

#Look at T cell division (CD3E pos/neg) in UMAP
lapply(names(data), function(i){
  x <- data[[i]]
  data.frame(
    sample = i,
    UMAP1 = x$UMAP1,
    UMAP2 = x$UMAP2,
    CD3Epos = CD3Epos[[i]],
    GZMKpos = GZMKpos[[i]],
    CD79Bneg = CD79Bneg[[i]],
    CD4 = x$smooth_CD4,
    CD8 = x$smooth_CD8B,
    GZMK= x$smooth_GZMK,
    stringsAsFactors = FALSE)
}) %>% 
  bind_rows() %>%
  mutate(class=case_when(
                        CD3Epos & GZMKpos & CD79Bneg ~ "T_tox",
                        CD3Epos & !GZMKpos & CD79Bneg ~ "T_CD4",
                        TRUE ~ "other")) %>%
ggplot(aes( UMAP1, UMAP2,
            col= sigmoid(GZMK^.15,  0.3)))+
  geom_point( size=0.2)+
  facet_wrap(~sample)+
  scale_color_gradient2(midpoint=0.5)

lapply(names(data), function(i){
  x <- data[[i]]
  data.frame(
    sample = i,
    UMAP1 = x$UMAP1,
    UMAP2 = x$UMAP2,
    CD3Epos = CD3Epos[[i]],
    GZMKpos = GZMKpos[[i]],
    CD79Bneg = CD79Bneg[[i]],
    CD4 = x$smooth_CD4,
    CD8 = x$smooth_CD8B,
    stringsAsFactors = FALSE)
}) %>% 
  bind_rows() %>%
  mutate(class=case_when(
    CD3Epos & GZMKpos & CD79Bneg ~ "T_tox",
    CD3Epos & !GZMKpos & CD79Bneg ~ "T_CD4",
    TRUE ~ "other")) %>%
  ggplot(aes( UMAP1, UMAP2,
              col= class))+
  geom_point( size=0.2)+
  facet_wrap(~sample)


data2 <- data %>%
  bind_rows( .id="sample" ) %>%
  as_tibble 

cellnames <- lapply(samplenames, colnames)
  cellnames <- do.call(c, (cellnames))
  data2 <- add_column(data2, cellnames)


data2%>%
  mutate(smooth_MAST1 = if_else(smooth_MAST1>.9, NA_real_,smooth_MAST1 ))%>%
ggplot +
  geom_point(aes( 
    x = UMAP1, y = UMAP2, 
    col = smooth_MAST1^.15 ) , size=.2) +
  facet_wrap( ~ sample ) +
  scale_color_gradientn( colours=rev(rje::cubeHelix(10))[2:10] )


data %>%
  bind_rows( .id="sample" ) %>%
  as_tibble %>%
ggplot +
  geom_point( aes( 
    x = raw_GZMK^.15,
    y = raw_CD70^.15 ), size=.1, alpha=.5 ) +
  facet_wrap( ~ sample )


##Cytotoxic T-Cells
#Without DLBCL1, no T-cells
#Create Pseudobulk object
#1. Filter for CD3E pos, GZMK pos, BCellmarker neg
#2. Sum up raw UMIs for each gene (rowSums) found in data$sample$raw_gene -> only subset of all genes,
#need to work with raw counts
#3. Create Matrix from vectors

gene_overlap <- reduce(lapply(samplenames, function(i) {
  rownames(i)
}), intersect)

Pseudobulk <- lapply(names(data),  function(samplename) {
  Matrix::rowSums( samplenames[[ samplename ]] [ gene_overlap , CD3Epos[[ samplename ]] 
                                                 & GZMKpos[[ samplename ]] & CD79Bneg[[ samplename ]] ]) 
  
}) %>%
 do.call(cbind, .)
  colnames(Pseudobulk) <- names(data)
                  
#Create colData
 Condition <-  c("aggressive", "aggressive", "aggressive", 
    "indolent", "indolent", "indolent", "indolent",
    "control", "control", "control", "aggressive", "aggressive")
coldata <- data.frame(Condition = Condition[2:12 ], row.names = names(data)[2:12 ], Sex = meta$Sex[2:12])

#Create DESeq object
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = Pseudobulk[, 2:12],
                              colData = coldata,
                              design = ~ Sex + Condition)
dds$Condition <- relevel(dds$Condition, "control")
#Error because of 0 in DLBCL3, because locmodes failed to set sensible threshhold, no cells remained
#now threshhold with ^0.5 power transformation to replace for CD79B
#Pre-filtering to get ridd of NAs
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <-  results(dds, contrast=c("Condition", "aggressive", "indolent"))

#Visualising results
plotMA(res, ylim=c(-2,2), main="Cytotoxic T Cells")

## CD4 Cells

#Create Pseudobulk object
#1. Filter for CD3E pos, GZMK pos, BCellmarker neg
#2. Sum up raw UMIs for each gene (rowSums) found in data$sample$raw_gene -> only subset of all genes,
#need to work with raw counts
#3. Create Matrix from vectors

Pseudobulk2 <- lapply(names(data),  function(samplename) {
  Matrix::rowSums( samplenames[[ samplename ]] [ gene_overlap , !CD3Epos[[ samplename ]] 
                                                 & !GZMKpos[[ samplename ]] & !CD79Bneg[[ samplename ]] ]) 
  
}) %>%
  do.call(cbind, .)
colnames(Pseudobulk2) <- names(data)

#Create colData
Condition <-  c("aggressive", "aggressive", "aggressive", 
                "indolent", "indolent", "indolent", "indolent",
                "control", "control", "control", "aggressive", "aggressive")
Sex <- meta$Sex[2:12]
coldata2 <- data.frame(Condition, Sex, row.names = names(data))

#Create DESeq object
library("DESeq2")
dds2 <- DESeqDataSetFromMatrix(countData = Pseudobulk2,
                              colData = coldata2,
                              design = ~ Condition+ Sex)
dds2$Condition <- relevel(dds2$Condition, "control")
#Error because of 0 in DLBCL3, because locmodes failed to set sensible threshhold, no cells remained
#now threshhold with ^0.5 power transformation to replace

#Pre-filtering to get rid of NAs
keep2 <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep2,]
dds2 <- DESeq(dds2)
res2 <-  results(dds2, contrast=c("Condition", "aggressive", "indolent"))

#Visualising results
plotMA(res2, ylim=c(-2,2), main="CD4 Cells")



data.frame(UMAP1=data$FL4$UMAP1, UMAP2=data$FL4$UMAP2, 
           CD8 = data$FL4$raw_CD8B, colSums= data$FL4$colsums, CD4=data$FL4$raw_CD4, GZMK=data$FL4$raw_GZMK)%>%
ggplot(aes(UMAP1, UMAP2, col=(CD8/colSums)^.5)) +
geom_point() +
  scale_color_gradientn(colours = rev(rje::cubeHelix(10))[2:10])


#Histograms for CD4 T cells of MAST1
#Filter for Cells of interest
gene <- "ELFN1-AS1"
plot_hist_gene <- function(gene) {
CD4cellsMAST1 <- lapply(names(data),  function(samplename) {
   CD4cells <-samplenames[[ samplename ]] [ gene_overlap , !CD3Epos[[ samplename ]] 
                                                 & !GZMKpos[[ samplename ]] & !CD79Bneg[[ samplename ]] ]
   CD4cellsMAST1 <- CD4cells[ gene, ]
  })%>%
  unlist()%>%
   tibble::enframe(name = "cellnames")%>%
  left_join( data2)
  CD4cellsMAST1 <- filter(CD4cellsMAST1, CD4cellsMAST1[[paste0("smooth_", gene)]] < .9)
  
  CD4cellsMAST1 <- mutate(CD4cellsMAST1, state=case_when(
    str_starts(CD4cellsMAST1$sample, "F") ~ "indolent",
    str_starts(CD4cellsMAST1$sample, "r") ~ "control",
    TRUE ~ "aggressive"
  ))
  
  ggplot(CD4cellsMAST1)+
  geom_density(aes(CD4cellsMAST1[[paste0("smooth_", gene)]]^.15, col=state))+
  facet_wrap(~sample)
    }

plot_hist_gene("MAST1")
plot_hist_gene("ELFN1-AS1")


#Compare whether DESeq2 result and CD4 cell markers match
CD4markers <- c("IL7R","PLAC8",
                "KLF2",
                "CCL5",
                "NKG7",
                "CCL4",
                "PDCD1",
                "TOX",
                "TOX2",
                "CD200",
                "CXCR5",
                "ICOS",
                "IL2RA",
                "FOXP3")
k <- res2[order(res2$padj), ]

CD4markers %in% rownames(k)[k$padj < .15] #only PLAC8 above .14 padj

#Clustering of CD4 T cells
#TH <- PLAC8pos, !PDCD1pos, !IL2RApos
#TFH <- !PLAC8pos, PDCD1pos, !IL2RApos
#TReg <- !PLAC8pos, !PDCD1pos, IL2RApos

PLAC8pos <- sapply(names(data), filter_gene_pos, "PLAC8")
PDCD1pos <- sapply(names(data), filter_gene_pos, "PDCD1")
IL2RApos <- sapply(names(data), filter_gene_pos, "IL2RA")


plot_hist_gene_celltypes <- function(gene) {
  CD4cellsMAST1 <- lapply(names(data),  function(samplename) {
    CD4cells <-samplenames[[ samplename ]] [ gene_overlap , !CD3Epos[[ samplename ]] 
                                             & !GZMKpos[[ samplename ]] & !CD79Bneg[[ samplename ]] ]
    CD4cellsMAST1 <- CD4cells[ gene, ]
  })%>%
    unlist()%>%
    tibble::enframe(name = "cellnames")%>%
    left_join( data2)
  CD4cellsMAST1 <- filter(CD4cellsMAST1, CD4cellsMAST1[[paste0("smooth_", gene)]] < .9)
  
  CD4cellsMAST1 <- mutate(CD4cellsMAST1, state=case_when(
    str_starts(CD4cellsMAST1$sample, "F") ~ "indolent",
    str_starts(CD4cellsMAST1$sample, "r") ~ "control",
    TRUE ~ "aggressive"
  ))
  CD4cellsMAST1 <- mutate(CD4cellsMAST1, state=case_when(
    CD4cellsMAST1$smooth_PLAC8 > location ~ "indolent",
    str_starts(CD4cellsMAST1$sample, "r") ~ "control",
    TRUE ~ "aggressive"
  ))
  ggplot(CD4cellsMAST1)+
    geom_density(aes(CD4cellsMAST1[[paste0("smooth_", gene)]]^.15, col=state))+
    facet_wrap(~sample)
}
  
paste0(gene, "pos") <- data[[sample]][[paste0("smooth_", gene)]] > location[[sample]][[paste0("thresh", gene)]]
  
