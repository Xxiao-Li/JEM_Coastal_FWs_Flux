
rm(list = ls())
graphics.off()
library(fluxweb)
library(igraph)

setwd("E:/00-food web functioning and stability/Data_manuscript/SG")
info <- read.csv("SG_fluxweb.csv" ,fileEncoding = "UTF-8-BOM")
info

diet <- read.csv("SG_diet.csv", header = TRUE, fileEncoding = "UTF-8-BOM")
diet <- as.matrix(diet)
diet[is.na(diet)] <- 0
diet <- apply(diet[,-1],2,as.numeric)
rownames(diet) <- colnames(diet) 
head(diet)
diet

#-----------------------------all functions---------------------------
plotfw=function(net, col=NULL, lab=NULL, size=NULL,
         nylevel=7, maxsize=10, labcex=0.01,
         ynum=6, ylab= "Trophic Level", ...){
  n <- vcount(net)
  if (!is.null(col)){
    V(net)$color <- col
  } else{
    V(net)$color <- rainbow(vcount(net))
  }
  if (!is.null(lab)){
    V(net)$name <- lab
  }
  if (!is.null(size)){
    V(net)$size <- size
  } else {
    V(net)$size <- maxsize/2
  }
  
  tl <- trophiclevels(net)
  dgpred <- tl
  
  bks <- c(0.9, seq(1.9, max(tl), length.out = nylevel))
  ynod <- cut(tl, breaks = bks, include.lowest = TRUE, 
              labels = 1:(length(bks)-1))
  
  maxx <- max(table(ynod))
  xnod <- rep(0,n)
  for (i in 1:nylevel){
    l <- sum(ynod==i)
    
    ltr <- (l/maxx)**(1/2)*maxx
    if (l>1) {
      xnod[ynod==i] <- seq(-ltr,ltr,length.out = l)
    } else {
      xnod[ynod==i] <- 0
    }
  }
  
  coo <- cbind(xnod,tl)
  
  #y axis with 1 and continuous axis from 2 to max TL.
  yax <- c(1, seq(2, max(tl), length.out = ynum-1))
  labax <- round(yax, 1)
  #rescale xax between -1 and 1
  laby <- (yax-min(yax))/(max(yax)-min(yax))*2 - 1
  
  plot(net, layout= coo, vertex.label.color="black", 
       vertex.label.cex=labcex, ...)
  axis(2, at = laby, labels = labax)
  mtext(ylab, side = 2, line=2)
  res <- data.frame(coo, "size"= V(net)$size, "color"= V(net)$color)
  names(res) <- c("x", "y", "size", "color")
  row.names(res) <- V(net)$name
  invisible(res)
}

trophiclevels = function(net){
  #Get adjacency matrix 
  mat <- get.adjacency(net, sparse=F)
  
  #Detect basal node
  basal <- rownames(subset(mat, apply(mat, 2, sum)==0) & apply(mat, 1, sum)!= 0)
  #Compute the shortest path to basal node
  paths_prey <- suppressWarnings(shortest.paths(graph = net, v= V(net), to = V(net)[basal], 
                                                mode = "in", weights = NULL, algorithm = "unweighted"))
  
  paths_prey[is.infinite(paths_prey)] <- NA
  shortest_paths <- suppressWarnings(as.matrix(apply(paths_prey, 1, min, na.rm=TRUE)))
  #for species with no prey apart of them
  shortest_paths[is.infinite(shortest_paths)] <- NA
  
  # Shortest TL
  sTL <- 1 + shortest_paths
  
  # Compute average prey trophic level
  # inspired from cheddar package calculation
  W <- t(mat)
  rs <- rowSums(W)
  W <- W/matrix(rs, ncol = ncol(W), nrow = nrow(W))
  W[0 == rs, ] <- 0
  I <- diag(ncol(W))
  tl0<-rowSums(I - W)
  result <- tryCatch(solve(I - W), error = function(e) e)
  if ("error" %in% class(result)) {
    avtl <- rep(NA, ncol(pm))
    names(avtl) <- colnames(pm)
  }
  else {
    avtl <- rowSums(result)
  }
  
  # Short-weighted TL is the average between 
  # Shortest TL and average prey TL
  SWTL <- (sTL + avtl)/2
  
  return(SWTL)
}

#-----------------------------------------------------------------------

diet_igraph <- graph.adjacency(diet,mode="directed",weighted=TRUE)
vcount(diet_igraph)
ecount(diet_igraph)

# basal specie
V(diet_igraph)$name[degree(diet_igraph, mode="in")==0]
# top predators 
V(diet_igraph)$name[degree(diet_igraph, mode="out")==0]
netmatrix <- get.adjacency(diet_igraph, sparse=F)
heatmap(netmatrix, Rowv=NA, Colv=NA, scale="none")

# ------2. Topological indicators ------------------------------------
# Species richness
unweighttopoly<-c()
S <- vcount(diet_igraph)
unweighttopoly$S <- S
# Connectance
C <- ecount(diet_igraph)/(S*(S-1))
unweighttopoly$C <- C
# Generality
# Identify predator nodes, i.e. taxa with at least one prey
pred <- degree(diet_igraph, mode="in")>0
# Compute mean generality of the food web, i.e. mean number of prey per predators
G <- sum(degree(diet_igraph, mode="in")[pred])/sum(pred)
G
unweighttopoly$G <- G
# Vulnerability
# Identify prey nodes, i.e. taxa with at least one predator
prey <- degree(diet_igraph, mode="out")>0
# Compute the mean vulnerability, i.e. mean number of predators per prey
V <- sum(degree(diet_igraph, mode="out")[prey])/sum(prey)
V
unweighttopoly$V <- V
# Mean shortest path
# shortest path length between all pair of nodes
sp <- shortest.paths(diet_igraph)
# Mean shortest path length between different species
# [upper.tri(sp)] remove the diagonal.
# Diagonal is by default set to 0 (path to itself)
# which artificially lower the mean shortest path.
ShortPath <- mean(sp[upper.tri(sp)]) 
ShortPath
unweighttopoly$ShortPath <- ShortPath
# Trophic level
# Compute the trophic level for each node  
tlnodes <- trophiclevels(diet_igraph)
unweighttopoly$tlnodes <- tlnodes
# Calculate the average trophic level of the food web
TL <- mean(tlnodes)
TL
unweighttopoly$TL <- TL
# Omnivory
# Link the trophic level to the interactions
webtl <- netmatrix*as.vector(tlnodes)
# Remove the trophic level when no interactions
webtl[webtl==0] <- NA
#Compute the standard of the trophic levels of prey 
omninodes <- apply(webtl,2,sd, na.rm=TRUE)
# Average the standard deviation over all taxa (with more than 2 preys)  
Omni <- mean(omninodes, na.rm=TRUE)
Omni
unweighttopoly$Omni <- Omni

write.csv(unweighttopoly, "unweighttopoly.csv")

#-Node weighted metrics------------------------------------------------
biomasses <- info$meanB
#Calculate percentage biomass per functional group
percB <- tapply(biomasses, info$trogroup, sum)/sum(biomasses)*100
barplot(as.matrix(percB), col=colFG, ylab="%")
#Visual representation parameter
Vscale <- 25 #multiplying factor
Vmin <- 4 #minimum size of the node
#scale the size of the node to the mean biomass
nodmax <- max(biomasses)
sizeB <- (biomasses/nodmax)*Vscale +Vmin
#Plot the food web
plotfw(diet_igraph, col=info$colfg, size=sizeB,
       edge.width=0.3, edge.arrow.size=0.3)
# Node-weighted Connectance
node_weighted <- c()
nwC <- sum(degree(diet_igraph)*biomasses)/(2*sum(biomasses)*(vcount(diet_igraph)-1))
nwC
node_weighted$nwC <- nwC
# Node-weighted generality
# Identify predators, i.e. taxa with at least one prey
pred <- degree(diet_igraph, mode="in")>0
# Compute weigthed in-degree average among predators
nwG <- sum((degree(diet_igraph, mode="in")*biomasses)[pred])/(sum(biomasses[pred]))
nwG
node_weighted$nwG <- nwG
# Node-weighted vulnerability
# Identify prey, i.e. taxa with at least one predator
prey <- degree(diet_igraph, mode="out")>0
# Compute weigthed out-degree average among prey
nwV <- sum((degree(diet_igraph, mode="out")*biomasses)[prey])/(sum(biomasses[prey]))
nwV
node_weighted$nwV <- nwV
#Node-weighted Trophic level
# tlnodes <- trophiclevels(net)
# Compute a weighted average of trophic levels
nwTL <- sum(tlnodes*biomasses)/sum(biomasses)
nwTL
node_weighted$nwTL <- nwTL

write.csv(node_weighted, "node_weighted.csv")

##---------------------------------------------------------------------

levels(info$types)
# Color the functional groups
colFG<- c("red","orange", "khaki", "blue", "green", "cyan")
# Assign the colors to each node (taxon)
info$colfg <- colFG[as.numeric(info$types)]
plotfw(diet_igraph, col=info$colfg, 
       edge.width=0.3, edge.arrow.size=0.3)


#--------------------calculating fluxes--------------------------------
##metabolic rate
nb.species <- length(info$trogroup)
nb.species
meta.rate <- rep(NA,nb.species)
meta.rate <- 0.71 * info$bodymass ^ (-0.25)
meta.rate[meta.rate == Inf] <- 0
meta.rate
losses <- meta.rate


##efficiencies
efficiencies = rep (NA, nb.species)
efficiencies[info$loss.type == "detritus"] <- 0.158
efficiencies[info$loss.type == "plant"] <- 0.545
efficiencies[info$loss.type == "animals"] <- 0.906
efficiencies
##calculating fluxes

biomasses <- info$meanB
mat.fluxes <- fluxing (diet,biomasses,losses,efficiencies, bioms.prefs = FALSE,bioms.losses = TRUE,ef.level = "prey")
mat.fluxes
write.csv()

basals <- colSums(mat.fluxes) == 0
names(basals)
plants <- basals
plants[which(names(basals) == "waterdetritus"  |names(basals) == "sedimentdetritus"  ) ] = FALSE
herbivory <- sum(rowSums(mat.fluxes[plants,]))
herbivory
carnivory <- sum(rowSums(mat.fluxes[!basals,]))
carnivory
detritivory <- sum(mat.fluxes[(names(basals) == "waterdetritus"  | names(basals) == "sedimentdetritus"),])
detritivory
total <- sum(mat.fluxes)
#_____________________________________________________________________

wightedtopology = function(fluxes, loop=FALSE){
  res<-c()
  # The flux matrix
  W.net<-as.matrix(fluxes) #fluxmatrix from fluxweb
  
  res$trophiclevel <- tlnodes
  res$metaboliclosses <- losses * biomasses
  res$assimilation <- efficiencies
  
  ### Taxon-specific Shannon indices of inflows
  # sum of k species inflows --> colsums
  sum.in<-apply(W.net, 2, sum)
  res$ingestion <- sum.in

  # Diversity of k species inflows
  # columns divided by the total col sum
  H.in.mat<- t(t(W.net)/sum.in)*t(log(t(W.net)/sum.in)) 
  H.in.mat[!is.finite(H.in.mat)] <- 0 #converts NaN to 0's
  H.in<- apply(H.in.mat, 2, sum)*-1
  
  # Effective number of prey or resources = N(R,k) 
  # The reciprocal of H(R,k) --> N (R,k) is the equivalent number of prey for species k
  N.res<-ifelse(sum.in==0, H.in, exp(H.in))
 
    ### Taxon-specific Shannon indices of outflows
  # sum of k speies outflows --> rowsums
  sum.out<-apply(W.net, 1, sum)
  
  res$outflow <- sum.out
  res$production <- sum.out + res$metaboliclosses
  res$excretion <- sum.in - sum.out - res$metaboliclosses
  res$PB <- res$production / biomasses
  res$QB <- res$ingestion / biomasses  
  res$efficiencies <- res$production / res$ingestion
  
  # Diversity of k species outflows
  # rows divided by the total row sum
  H.out.mat<- (W.net/sum.out)*log(W.net/sum.out) 
  H.out.mat[!is.finite(H.out.mat)] <- 0 #converts NaN to 0's
  H.out<- apply(H.out.mat, 1, sum)*-1
  
  # Effective number of predators or consumers = N(C,k) 
  # The reciprocal of H(C,k) --> N (C,k) is the equivalent number of predators for species k
  N.con<-ifelse(sum.out==0, H.out, exp(H.out))
  
  ### Quantitative Weighted connectance
  no.species<-ncol(W.net)
  
  # The weighted link density (LDw) is:
  # In the weighted version the effective number of predators for species i is weighted by i's 
  # contribution to the total outflow the same is the case for the inflows
  tot.mat<- sum(W.net)
  # LD.w <- (sum((sum.in/tot.mat)*N.res) + sum((sum.out/tot.mat)*N.con))/2
  # equivalent to next formula, but next one is closer to manuscript
  LD <- 1/(2*tot.mat)*(sum(sum.in*N.res) + sum(sum.out*N.con))
  
  #Weighted connectance
  res$lwC<- LD/ifelse(loop, no.species, no.species-1)
  
  # positional.index
  pos.ind<- sum.in*N.res/(sum.in*N.res+sum.out*N.con) #postional index
  basal.sp<- pos.ind[pos.ind==0] #basal species = 0
  top.sp<- pos.ind[pos.ind==1] #defintion according to Bersier et al. 2002 top species = [0.99, 1]
  
  con.sp<-length(pos.ind)-length(basal.sp)# all consumer taxa except basal
  # weighted quantitative Generality
  res$lwG<-sum(sum.in*N.res/sum(W.net))
  
  res.sp<- length(pos.ind)-length(top.sp)
  # weighted quantitative Vulnerability
  res$lwV<-sum(sum.out*N.con/sum(W.net))
  
  return(res)
}

flux.output <- wightedtopology(mat.fluxes)
write.csv(flux.output, "linkwightedtopology.csv")  #all flux metrics
#------------------------------------------------------------------

growth.rates = rep(NA, dim(mat.fluxes)[1])
growth.rates[colSums(mat.fluxes) == 0] = 0.1

stability.value(mat.fluxes,biomasses,losses,efficiencies,growth.rates,bioms.prefs = FALSE, bioms.losses = TRUE, ef.level = "prey",
                full.output = FALSE)

write.csv(out$values,"eigenvalues.csv")
max(Re(out$values))



