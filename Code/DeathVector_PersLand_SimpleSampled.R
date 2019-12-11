### Introduction to Topological Data Analysis with R 
### Peter Bubenik
### October 30, 2019

# The following three commands only need to be run once to install
# the required R packages on your system. Once you've run them once
# please add a "#" symbol to the start of each line.
#install.packages("TDA") # R Topological Data Analysis package
#install.packages("deldir") # R Delaunay and Voronoi package
#install.packages("kernlab") # R Support Vector Machine package
#install.packages("png") # R image processing
#install.packages("imager") # R image processing

# The following commands load the packages functions into memory.
library(TDA) 
library(deldir) 
library("Matrix") # package for sparse matrices
library("kernlab") 
library(png) # For loading images
#library(imager) # For processing images
par(pty="s") # force the plotting region to be square

###############################
########## Functions ##########
###############################

# Euclidean distance between vectors
euclidean.distance <- function(u, v) sqrt(sum((u - v) ^ 2))

# Sample points from an annulus
sample.annulus <- function(num.pts, inner.radius, outer.radius){
  theta <- runif(num.pts) * 2 * pi
  radius <- sqrt(runif(num.pts, inner.radius^2, outer.radius^2))
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  cbind(x,y)
}

# Plot the Voronoi cells and dual and Delaunay complex
plot.delaunay <- function(X){
  DelVor <- deldir(X[,1], X[,2], suppressMsge = TRUE)
  # Voronoi cells:
  plot(DelVor, pch=20, col=c("black","red","blue"), wlines= ("tess"))
  # Voronoi cells and their dual (the Delaunay complex):
  plot(DelVor, pch=20, col=c("black","red","blue"), wlines= ("both"))
  # Delaunay complex:
  plot(DelVor, pch=20, col=c("black","red","blue"), wlines= ("triang"))
}

# Plot the Vietoris-Rips complex up to some filtration value
plot.rips <- function(X,max.filtration){
  plot(X, pch=20, col="blue", asp=1)
  num.pts <- dim(X)[1]
  for(i in 1:num.pts)
    for(j in 1:num.pts)
      if (euclidean.distance(X[i,],X[j,]) < max.filtration)
        lines(rbind(X[i,],X[j,]))
}

# Plot representative cycles for Delaunay complex
plot.delaunay.cycle <- function(X){
  PH.output <- alphaComplexDiag(X, maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                location = TRUE)
  PD <- PH.output[["diagram"]]
  ones <- which(PD[, 1] == 1)
  persistence <- PD[ones,3] - PD[ones,2]
  cycles <- PH.output[["cycleLocation"]][ones[order(persistence)]]
  for (i in 1:length(cycles)){
    plot(X, pch=20, col="blue", asp=1)
    for (j in 1:dim(cycles[[i]])[1])
      lines(cycles[[i]][j,,])
  }
}

### Plot representative cycles for Vietoris-Rips complex
plot.rips.cycle <- function(X){
  PH.output <- ripsDiag(X, maxdimension = 1, maxscale = max.filtration, 
                        library = c("GUDHI", "Dionysus"), location = TRUE)
  PD <- PH.output[["diagram"]]
  ones <- which(PD[, 1] == 1)
  persistence <- PD[ones,3] - PD[ones,2]
  cycles <- PH.output[["cycleLocation"]][ones[order(persistence)]]
  for (i in 1:length(cycles)){
    plot(X, pch=20, col="blue", asp=1)
    for (j in 1:dim(cycles[[i]])[1])
      lines(cycles[[i]][j,,])
  }
}

# Death Vector
death.vector <- function(PD){
  zeroes <- which(PD[, 1] == 0)
  PD.0 <- PD[zeroes,2:3]
  dv <- vector()
  if ((min(PD.0[,"Birth"]) == 0) && (max(PD.0[,"Birth"]) == 0))
    dv <- sort(PD.0[,2], decreasing=TRUE)
  return(dv)
}

# Plot Persistence Landscape
plot.landscape <- function(PL,t.vals){
  plot(t.vals,PL[1,],type="l",ylab="Persistence",xlab="Parameter values",col=1,ylim=c(min(PL),max(PL)))
  for(i in 2:dim(PL)[1])
    lines(t.vals,PL[i,],type="l",col=i)
}

# Matrix of death vectors from a list of persistence diagrams
#death.vector.matrix <- function(PD.list){
#  print(length(PD.list))
#  num.pts <- length(which(PD.list[[1]][,1] == 0))
#  print(num.pts)
#  DVM <- matrix(0L, nrow = length(PD.list), ncol = num.pts - 1)
#  for (c in 1 : length(PD.list)) {
#    DVM[c,] <- death.vector(PD.list[[c]])[-1]
#  }
#  cat(DVM[c,])
#  return(DVM)
#}

#Matrix of death vectors from a list of persistence diagrams
#death.vector.matrix <- function(PD.list){
#  print(length(PD.list))
#  num.pts <- length(which(PD.list[[1]][,1] == 0))
#  print(num.pts)
#  DVM <- matrix(0L, nrow = length(PD.list), ncol = num.pts - 1)
#  print(dim(DVM))
#  print(length(DVM[1,]))
#  for (c in 1 : length(PD.list)) {
#    v = death.vector(PD.list[[c]])[-1]
#    while(length(v) < num.pts-1) {
#      v <- c(v, v[length(v)])
#    }
#    DVM[c,] <- v
#  }
#  return(DVM)
#}

# MODIFIED: 
# Matrix of death vectors from a list of persistence diagrams
death.vector.matrix <- function(PD.list, num.pts){
  DVM <- matrix(0L, nrow = length(PD.list), ncol = num.pts)
  for (c in 1 : length(PD.list)) {
    v = death.vector(PD.list[[c]])[-1]
    while(length(v) < num.pts) {
      cat(v)
      v <- c(v, v[length(v)])
    }
    DVM[c,] <- v
  }
  return(DVM)
}

# Matrix of persistence landscape row vectors from list of persistence landscapes
landscape.matrix.from.list <- function(PL.list){
  n <- length(PL.list)
  m <- ncol(PL.list[[1]])
  max.depth <- integer(n)
  for (i in 1:n)
    max.depth[i] <- nrow(PL.list[[i]])
  K <- max(max.depth)
  PL.matrix <- Matrix(0, nrow = n, ncol = m*K, sparse = TRUE)
  for (i in 1:n)
    for (j in 1:max.depth[i])
      PL.matrix[i,(1+(j-1)*m):(j*m)] <- PL.list[[i]][j,]
  return(PL.matrix)
}

# Convert a vector to a persistence landscape
landscape.from.vector <- function(PL.vector, t.vals){
  m <- length(t.vals)
  K <- length(PL.vector)/m
  PL <- Matrix(0, nrow = K, ncol=m, sparse = TRUE)
  for (i in 1:K){
    PL[i,1:m] <- PL.vector[(1+(i-1)*m):(i*m)]
  }
  return(PL)
}

# Take difference of vectors of potentially different lengths
difference.vectors <-function(vector.1,vector.2){
  length.1 <- length(vector.1)
  length.2 <- length(vector.2)
  difference.vector = numeric(max(length.1,length.2))
  difference.vector = as(difference.vector, "sparseVector")
  difference.vector[1:length.1] = difference.vector[1:length.1] + vector.1
  difference.vector[1:length.2] = difference.vector[1:length.2] - vector.2
}

# Permutation test for two matrices consisting of row vectors
permutation.test <- function(M1 ,M2, num.repeats = 10000){
  # append zeros if necessary so that the matrices have the same number of columns
  num.columns <- max(ncol(M1),ncol(M2))
  M1 <- cbind(M1, Matrix(0,nrow=nrow(M1),ncol=num.columns-ncol(M1)))
  M2 <- cbind(M2, Matrix(0,nrow=nrow(M2),ncol=num.columns-ncol(M2)))
  t.obs <- euclidean.distance(colMeans(M1),colMeans(M2))
  k <- dim(M1)[1]
  M <- rbind(M1,M2)
  n <- dim(M)[1]
  count <- 0
  for (i in 1:num.repeats){
    permutation <- sample(1:n)
    t <- euclidean.distance(colMeans(M[permutation[1:k],]),colMeans(M[permutation[(k+1):n],]))
    if (t >= t.obs)
      count <- count + 1
  }
  return(count/num.repeats)
}

# Import sampled points from .csv file
import.sampled.points <- function(boundary_fname, dark_fname){
  # Import image boundary data from .csv file
  df_boundary <- read.csv(boundary_fname, header=FALSE)
  # Extract x and y coordinates
  boundary_x <- as.vector(df_boundary['V1'])
  boundary_y <- as.vector(df_boundary['V2'])
  # Bind x and y coordinates into single matrix
  boundary = cbind(unlist(boundary_x), unlist(boundary_y))
  # Creat plot to view boundary
  plot(boundary, pch=20, col="blue", asp=1)
  
  # Import image interior data from .csv file
  df_interior <- read.csv(dark_fname, header=FALSE)
  # Extract x and y coordinates
  interior_x <- as.vector(df_interior['V1'])
  interior_y <- as.vector(df_interior['V2'])
  # Bind x and y coordinates into single matrix
  interior = cbind(unlist(interior_x), unlist(interior_y))
  # Create plot to view interior
  plot(interior, pch=20, col="blue", asp=1)
  
  # Bind boundary and interior into single matrix and divide by 150 so all x,y coordinates have max value of 1
  X <- rbind(boundary, interior)/150
  
  # Rename column and row names
  colnames(X) <- c("x","y")
  rownames(X) <- NULL
  
  # Return matrix
  return(X)
}

################################
########## Parameters ##########
################################

num.pts <- 200
outer.radius <- 1
inner.radius.max <- 0.5
inner.radius.min <- 0.25
num.repeats <- 32
t.steps <- 300 #200
min.t <- 0
max.t <- 1
t.vals <- seq(min.t,max.t,(max.t-min.t)/t.steps)
max.filtration <- 1.1 # halt the Rips computation at this value
cost <- 10
num.folds <- 5
num.training <- 24 #21
num.testing <- 8 #20
num.repeats.for.average <- 20
epsilon <- 0.01


############################
########## Images ##########
############################
# Set working directory
setwd('/home/james/Documents/UF-Homework/MTG7396-Topological-Data-Analysis/Project')

# Get all filenames containing image boundary data
boundary_files <- list.files(path="Simple-Sampled-Data-and-Plots", pattern="*_gray_boundary.csv", full.names=TRUE, recursive=FALSE)
boundary_files
# Split into three different classes
bkl_boundary_files = boundary_files[1:41]
bkl_boundary_files
mel_boundary_files = boundary_files[42:81]
mel_boundary_files
nv_boundary_files = boundary_files[82:141]
nv_boundary_files

# Get all filenames containing image dark interior data
dark_interior_files <- list.files(path="Simple-Sampled-Data-and-Plots", pattern="*_gray_dark.csv", full.names=TRUE, recursive=FALSE)
dark_interior_files
# Split into three different classes
bkl_dark_interior_files = dark_interior_files[1:41]
bkl_dark_interior_files
mel_dark_interior_files = dark_interior_files[42:81]
mel_dark_interior_files
nv_dark_interior_files = dark_interior_files[82:141]
nv_dark_interior_files



########################################################
########## COMPUTE DV and PL for SINGLE IMAGE ##########
########################################################

# Import image boundary data from .csv file
X <- import.sampled.points(mel_boundary_files[38], mel_dark_interior_files[38])
X
plot(X, pch=20, col="blue", asp=1)

# Delaunay Complex and its Persistence Diagram
# plot the Voronoi cells and Delaunay complex (3 plots)
plot.delaunay(X)
# compute persistent homology
PH.output <- alphaComplexDiag(X)
PD <- PH.output[["diagram"]]
# filtration values have been squared for some reason - fix by taking square root
PD[,2:3] <- sqrt(PD[,2:3])
# plot the persistence diagram
plot(PD, asp=1, diagLim = c(0,max.t))
legend(0.5*max.t, 0.25*max.t, c("Homology in degree 0","Homology in degree 1"), 
       col = c(1,2), pch = c(19,2), cex = .8, pt.lwd = 2)
# plot the bar code
plot(PD, diagLim = c(0,max.t), barcode=TRUE)

# Vietoris-Rips complex and its Persistence Diagram
# plot the Vietoris-Rips complex up to the maximum filtration value
plot.rips(X,max.filtration)
# Compute persistent homology
PH.output <- ripsDiag(X, maxdimension = 1, maxscale = max.filtration)
# Get persistence diagram from output of persistent homology computation
PD <- PH.output[["diagram"]]
# plot the persistence diagram
plot(PD, asp=1, diagLim = c(0, max.filtration))
legend(0.5*max.filtration, 0.25*max.filtration, c("Homology in degree 0","Homology in degree 1"),
       col = c(1,2), pch = c(19,2), cex = .8, pt.lwd = 2)
# plot the bar code
plot(PD, diagLim = c(0, max.filtration), barcode=TRUE)

# plot representative most-persistent cycle from Delaunay complex
plot.delaunay.cycle(X)

# plot representative most-persistent cycle from Vietoris-Rips complex
plot.rips.cycle(X)

# Use Delaunay complex from now on
PH.output <- alphaComplexDiag(X)
PD <- PH.output[["diagram"]]
PD[,2:3] <- sqrt(PD[,2:3])

# Death Vector
DV <- death.vector(PD)
DVr <- DV[-1]
plot(DVr, type="l", col="blue", ylab="Persistence")

# Persistence Landscapes in dimensions 1
PL <- t(landscape(PD,dimension=1,KK=1:100,tseq=t.vals))
plot.landscape(PL,t.vals)







###########################################################
########## COMPUTE DV and PL for MULTIPLE IMAGES ##########
###########################################################

# Average Death Vector and Average Persistence Landscape
num.repeats = 40
PD.list <- vector("list",num.repeats)
for (c in 1 : num.repeats){
  #X <- sample.annulus(num.pts,inner.radius.max,outer.radius)
  X <- import.sampled.points(nv_boundary_files[c], nv_dark_interior_files[c])
  #PH.output <- alphaComplexDiag(X)
  PH.output <- ripsDiag(X, maxdimension = 1, maxscale = max.filtration)
  PD.list[[c]] <- PH.output[["diagram"]]
  PD.list[[c]][,2:3] <- sqrt(PD.list[[c]][,2:3])
}

DV.matrix <- death.vector.matrix(PD.list, 199)
average.DV <- colMeans(DV.matrix)
plot(average.DV, type="l", col="blue", ylab = "Persistence")

PL.list <- vector("list",num.repeats)
for (c in 1 : num.repeats) {
  PL.list[[c]] <- t(landscape(PD.list[[c]],dimension=1,KK=1:100,tseq=t.vals))
}
PL.matrix <- landscape.matrix.from.list(PL.list)
average.PL.vector <- colMeans(PL.matrix, sparseResult = TRUE)
average.PL <- landscape.from.vector(average.PL.vector,t.vals)
plot.landscape(average.PL,t.vals)

for (c in 1 : num.repeats) {
  plot.landscape(PL.list[[c]],t.vals)
}

PL.matrix

# Second class of samples
PD.list.2 <- vector("list",num.repeats)
for (c in 1 : num.repeats){
  #X <- sample.annulus(num.pts,inner.radius.min,outer.radius)
  X <- import.sampled.points(mel_boundary_files[c], mel_dark_interior_files[c])
  #PH.output <- alphaComplexDiag(X)
  PH.output <- ripsDiag(X, maxdimension = 1, maxscale = max.filtration)
  PD.list.2[[c]] <- PH.output[["diagram"]]
  PD.list.2[[c]][,2:3] <- sqrt(PD.list.2[[c]][,2:3])
}

DV.matrix.2 <- death.vector.matrix(PD.list.2, 199)
average.DV.2 <- colMeans(DV.matrix.2)
plot(average.DV.2, type="l", col="blue", ylab = "Persistence")

PL.list.2 <- vector("list",num.repeats)
for (c in 1 : num.repeats) {
  PL.list.2[[c]] <- t(landscape(PD.list.2[[c]],dimension=1,KK=1:100,tseq=t.vals))
}
PL.matrix.2 <- landscape.matrix.from.list(PL.list.2)
average.PL.vector.2 <- colMeans(PL.matrix.2, sparseResult = TRUE)
average.PL.2 <- landscape.from.vector(average.PL.vector.2,t.vals)
plot.landscape(average.PL.2,t.vals)

for (c in 1 : num.repeats) {
  plot.landscape(PL.list.2[[c]],t.vals)
}

mel_boundary_files


# Third class of samples
PD.list.3 <- vector("list",num.repeats)
for (c in 1 : num.repeats){
  #X <- sample.annulus(num.pts,inner.radius.min,outer.radius)
  X <- import.sampled.points(bkl_boundary_files[c], bkl_dark_interior_files[c])
  #PH.output <- alphaComplexDiag(X)
  PH.output <- ripsDiag(X, maxdimension = 1, maxscale = max.filtration)
  PD.list.3[[c]] <- PH.output[["diagram"]]
  PD.list.3[[c]][,2:3] <- sqrt(PD.list.3[[c]][,2:3])
}

DV.matrix.3 <- death.vector.matrix(PD.list.3, 199)
average.DV.3 <- colMeans(DV.matrix.3)
plot(average.DV.3, type="l", col="blue", ylab = "Persistence")

PL.list.3 <- vector("list",num.repeats)
for (c in 1 : num.repeats) {
  PL.list.3[[c]] <- t(landscape(PD.list.3[[c]],dimension=1,KK=1:100,tseq=t.vals))
}
PL.matrix.3 <- landscape.matrix.from.list(PL.list.3)
average.PL.vector.3 <- colMeans(PL.matrix.3, sparseResult = TRUE)
average.PL.3 <- landscape.from.vector(average.PL.vector.3,t.vals)
plot.landscape(average.PL.3,t.vals)


# The difference in death vectors and persistence landscapes
difference.average.DV <- average.DV - average.DV.2
plot(difference.average.DV, type="l", col="blue", ylab="Persistence")
difference.average.PL.vector <- difference.vectors(average.PL.vector,average.PL.vector.2)
difference.average.PL <- landscape.from.vector(difference.average.PL.vector,t.vals)
plot.landscape(difference.average.PL,t.vals)

# p values for differences in the average landscapes
permutation.test(DV.matrix,DV.matrix.2)
permutation.test(PL.matrix,PL.matrix.2,num.repeats=1000)



# Principal Components Analysis
data.labels <- c(rep(1,nrow(PL.matrix)), rep(2,nrow(PL.matrix.2)))
DV.vectors <- rbind(DV.matrix,DV.matrix.2)
pca.0 <- prcomp(DV.vectors,rank=10)
plot(pca.0,type="l")
plot(pca.0$x[,1:2], col=data.labels, pch=17+(2*data.labels), asp=1)
loading_vectors <- t(pca.0$rotation)
plot(loading_vectors[1,],type="l")
plot(loading_vectors[2,],type="l")
PL.vectors <- rbind(PL.matrix,PL.matrix.2)
pca.1 <- prcomp(PL.vectors,rank=10)
plot(pca.1,type="l")
plot(pca.1$x[,1:2], col=data.labels, pch=17+(2*data.labels), asp=1)
loading_vectors <- t(pca.1$rotation)
plot.landscape(landscape.from.vector(loading_vectors[1,],t.vals),t.vals)
plot.landscape(landscape.from.vector(loading_vectors[2,],t.vals),t.vals)

# View death vectors and persistence landscapes
DV.vectors
PL.vectors

# Save death vectors and persistence landscapes to file
DV.vectors.all <- rbind(DV.matrix,DV.matrix.2,DV.matrix.3)
PL.vectors.all <- rbind(PL.matrix,PL.matrix.2,PL.matrix.3)
write.table(as.matrix(DV.vectors.all),file="Simple-TDA-Data/Simple-DV-Data-Vietoris-Rips.csv",row.names=FALSE, col.names=FALSE)
write.table(as.matrix(PL.vectors.all),file="Simple-TDA-Data/Simple-PL-Data-Vietoris-Rips.csv",row.names=FALSE, col.names=FALSE)

# Support Vector Machine Classification
svm_model <- ksvm(as.matrix(DV.vectors),data.labels,type="C-svc",scaled=c(),kernel="rbfdot",C=cost,cross=num.folds)
print(svm_model)
svm_model <- ksvm(as.matrix(PL.vectors),data.labels,type="C-svc",scaled=c(),kernel="rbfdot",C=cost,cross=num.folds)
print(svm_model)
