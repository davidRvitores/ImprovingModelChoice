

############################################################################################################################
#                         BRIEF EXPLANATION OF THE MAIN ALGORITHMS INVOLVED                                                #
############################################################################################################################


# ----------------------------------------------------------------------------------------------------------------
# EXPLANATION OF THE ALGORITHMS FOR G-CPC
# ----------------------------------------------------------------------------------------------------------------

# GCPC(X,k,M,csh,cvol,niter,tol,nstart,graph)  -  FUNCTION FOR CLUSTERING G-CPC

# Input values:

#       X <- data matrix, containing by rows the observations
#       k <- number of clusters to search
#       M <- number of classes sharing the common parameters within each of them
#       csh, cvol <- constant values for the determinant and shape constraints
#       niter <- maximum number od iterations in the algorithms
#       tol <- tolerance for the convergence in the algorithms
#       nstart <- number of random initializations
#       graph <- T or F. If T, display the solution for each random initializations

# Output values - Returns a list with the following elements of the best solution found between the nstart:

#       par <- list containing the parameters of a model:
#             par$pi <- vector of length k containing the weights
#             par$mu <- matrix containing in the row i the mean of the group i
#             par$Lambda <- list with the k sets of eigenvalues (without imposing prod(Lambda[[i]])=1),
#                         we are including here the shapes and volumes
#             par$ejes <- list with the M common orthogonal matrices
#             par$P <- vector of length k with the group assignation 
#       posterior <- matrix with the posterior probabilities of each observation to arise from each group
#       loglik <- loglikelihood in the best solution
#       df <- degreees of freedom
#       BIC <- value of BIC

# ----------------------------------------------------------------------------------------------------------------

# GCPCclas(X,labels,M,csh,cvol,niter,tol,nstart,graph)  -  FUNCTION FOR DISCRIMINANT ANALYSIS G-CPC

# Input values: all parameters defined as above, and:

#       labels <- vector containing the labels (with levels 1,...,k)

# Output values: same as above

# ----------------------------------------------------------------------------------------------------------------
# EXPLANATION OF THE ALGORITHMS FOR G-PROP - the input and output values are the same as in G-CPC, the functions are:
# ----------------------------------------------------------------------------------------------------------------

# GPROP(X,k,M,csh,cvol,niter,tol,nstart,graph)  -  FUNCTION FOR CLUSTERING G-PROP

# GPROPclas(X,labels,M,csh,cvol,niter,tol,nstart,graph)  -  FUNCTION FOR DISCRIMINANT ANALYSIS G-PROP


# ----------------------------------------------------------------------------------------------------------------
# EXPLANATION OF THE FUNCTION FOR PLOTTING THE RESULTS OF A FITTED MODEL
# ----------------------------------------------------------------------------------------------------------------

# grafElipsesError(X,posterior,par,d,componentes,labels)

# X <- data set
# posterior <- matrix with the posterior probabilities of observation i arising from group j computed with the model
# par <- parameters computed by the model 
# d <- dimension of plot, take the first d components
# componentes <- axis taken for the representation. If componentes=0, canonical axis, and if 
#             componentes=1,...,G the axis par$ejes[[componentes]]
# labels <- original labels of the groups, used to mark the missclassified observations.



# ----------------------------------------------------------------------------------------------------------------
# EXPLANATION OF THE FUNCTIONS FOR CROSS VALIDATION
# ----------------------------------------------------------------------------------------------------------------


# cvGCPC(train,test,labelsTrain,labelsTest,par0,csh,cvol,niter,tol)
# cvPROP(train,test,labelsTrain,labelsTest,par0,csh,cvol,niter,tol)


# Functions to apply the cross validation techniques. Splitting the data set in test and train, these functions 
# compute the missclassification error in the test set with the model computed from the training set. 
# This function is  supposed to be applied after a model has been fitted to the total data set. In order to avoid
# the computational cost of iterating again over a large number of starts, we take as initial solution the 
# parameters par0 estimated in the global model, and compute only one iteration. 

# cvMclust(train,test,labelsTrain,labelsTest,modelo)

# Algorithm to apply the same methodology, with mclust. (with no initial solution trick, but mclust works much more faster)



############################################################################################################################
############################################################################################################################
###                                                                                                                      ###
###                                                 REAL DATA EXAMPLES                                                   ###
###                                                                                                                      ###
############################################################################################################################
############################################################################################################################


# In this script we compute the solution of the real data examples in the paper. The solutions presented in the 
# paper are the best solutions obtained from these algorithms, and they should be obtained taking a large number 
# of random initializations. In order to avoid long run-times, with the number of initializations we are taking
# optimal solutions could not be achieved, but reasonable solutions might also be achieved. 


############################################################################################################################
#                                          IRIS DATASET - CLUSTER ANALYSIS                                                 #      
############################################################################################################################

data(iris)
X <- as.matrix(iris[,1:4])
labels <- as.numeric(factor(as.character(iris[,5]),labels=1:3))

#------------------------------- 2-CPC ------------------------------------

iris2CPC <- GCPC(X,k=3,M=2,csh=100,cvol=100,niter=100,tol=10^(-5),nstart=10,graph=F)
grafElipsesError(X,iris2CPC$posterior,iris2CPC$par,d=4,0,labels)
title(paste("2-CPC CLUSTERING \n"," Loglik =", round(iris2CPC$loglik,4), "\n df = ",
            iris2CPC$df, "\n BIC=",round(iris2CPC$BIC,4)),line=-7,outer=T)
error(labels, apply(iris2CPC$posterior,1,which.max))

#------------------------------- 2-PROP ------------------------------------

iris2PROP <- GPROP(X,k=3,M=2,csh=100,cvol=100,niter=100,tol=10^(-5),nstart=10,graph=F)
grafElipsesError(X,iris2PROP$posterior,iris2PROP$par,d=4,0,labels)
title(paste("2-PROP CLUSTERING \n","Loglik =", round(iris2PROP$loglik,4), "\n df = ",
            iris2PROP$df, "\n BIC=",round(iris2PROP$BIC,4)),line=-7,outer=T)
error(labels, apply(iris2PROP$posterior,1,which.max))

#------------------------------- mclust  ------------------------------------

mod <- Mclust(X,G=3)
summary(mod)

############################################################################################################################
#                                      IRIS DATASET - DISCRIMINANT ANALYSIS                                                #      
############################################################################################################################

#------------------------------- 2-CPC ------------------------------------

iris2CPCb <- GCPCclas(X,labels,M=2,csh=100,cvol=100,niter=100,tol=10^(-6),nstart=5,graph=F)
grafElipsesError(X,iris2CPCb$posterior,iris2CPCb$par,d=4,0,labels)
title(paste("2-CPC DISCRIMINANT ANALYSIS \n"," Loglik =", round(iris2CPCb$loglik,4), "\n df = ",
            iris2CPCb$df, "\n BIC=",round(iris2CPCb$BIC,4)),line=-7,outer=T)
error(labels, apply(iris2CPCb$posterior,1,which.max))

#------------------------------- 2-PROP ------------------------------------

iris2PROPb <- GPROPclas(X,labels,M=2,csh=100,cvol=100,niter=100,tol=10^(-6),nstart=5,graph=F)
grafElipsesError(X,iris2PROPb$posterior,iris2PROPb$par,d=4,0,labels)
title(paste("2-PROP DISCRIMINANT ANALYSIS \n"," Loglik =", round(iris2PROPb$loglik,4), "\n df = ",
            iris2PROPb$df, "\n BIC=",round(iris2PROPb$BIC,4)),line=-7,outer=T)
error(labels, apply(iris2PROPb$posterior,1,which.max))

#------------------------------- mclust ------------------------------------

irisMclust <- MclustDA(X, labels, modelType = "EDDA")
summary(irisMclust)
modelo <- irisMclust$model$`1`$modelName


#----------------------------- Cross validation -------------------------------

N <- nrow(X)
k <- length(unique(labels))
v <- 1:N

# resamplingCV <- 0 # to perform K-fold (take 1<=K<=N) or leave one out (take K=N)
resamplingCV <- 1 # to perform CV(K,p)

K <- 50 # Take a suitable value of K for the method  
# In the paper, K=300 and p=0.8, 0.95 for CV, K=150 for L.O.O
p <- 0.8 # probability of knowing the label of an observation

particion <- list()
i0 <- rep(0,K)
i1 <- rep(0,K)
i2 <- rep(0,K)

cat("Partition: \n")
for(h in 1:K){
  cat(h,"-")
  
  if(resamplingCV==0){
    # Partition for k-fold o LOO:
    particion[[h]] <- sample(v,floor(length(v)/(K+1-h)))
    v <- setdiff(v,particion[[h]])
  }else{
    # Partition created assigning labels with probability p
    particion[[h]] <- v[runif(N)>p]
  }
  
  test <- X[particion[[h]],]
  labelsTest <- labels[particion[[h]]]
  train <- X[-particion[[h]],]
  labelsTrain <- labels[-particion[[h]]]
  
  i0[h] <- cvMclust(train,test,labelsTrain,labelsTest,modelo)
  i1[h] <- cvGCPC(train,test,labelsTrain,labelsTest,par0=iris2CPC$par,csh=100,cvol=100,niter=200,tol=10^(-7))
  i2[h] <- cvGPROP(train,test,labelsTrain,labelsTest,par0=iris2PROP$par,csh=100,cvol=100,niter=200,tol=10^(-7))
}

cat("IRIS Cross Validation results: \n",
    "Mclust: \t Mean error= \t",mean(i0[i0>=0]),"\t sd/sqrt(K)= \t",sd(i0[i0>=0])/sqrt(length(i0[i0>=0])),"\n",
    "2-CPC: \t Mean error= \t",mean(i1),"\t sd/sqrt(K)= \t",sd(i1)/sqrt(K),"\n",
    "2-PROP: \t Mean error= \t",mean(i2),"\t sd/sqrt(K)= \t",sd(i2)/sqrt(K),"\n")



############################################################################################################################
#                                      CRABS DATASET - DISCRIMINANT ANALYSIS                                               #      
############################################################################################################################

library(MASS)
data(crabs)
X <- as.matrix(crabs[4:8])

specSex0 <- as.character(paste(crabs[,1],crabs[,2]))
labels <- as.numeric(factor(specSex0,labels=1:4))

#------------------------------- 2-CPC ------------------------------------

crabs2CPC <- GCPCclas(X, labels, M=2, csh=100000, cvol=100000, niter=100, tol=10^(-5), nstart=5,graph=F)
grafElipsesError(X,crabs2CPC$posterior,crabs2CPC$par,d=5,0,labels)
title(paste("2-CPC DISCRIMINANT ANALYSIS \n"," Loglik =", round(crabs2CPC$loglik,4), "\n df = ",
            crabs2CPC$df, "\n BIC=",round(crabs2CPC$BIC,4)),line=-7,outer=T)
error(labels, apply(crabs2CPC$posterior,1,which.max))

#------------------------------- 2-PROP ------------------------------------

crabs2PROP <- GPROPclas(X, labels, M=2,csh=100000,cvol=100000,niter=100,tol=10^(-5),nstart=5,graph=F)
grafElipsesError(X,crabs2PROP$posterior,crabs2PROP$par,d=5,0,labels)
title(paste("2-PROP DISCRIMINANT ANALYSIS \n"," Loglik =", round(crabs2PROP$loglik,4), "\n df = ",
            crabs2PROP$df, "\n BIC=",round(crabs2PROP$BIC,4)),line=-7,outer=T)
error(labels, apply(crabs2PROP$posterior,1,which.max))

#------------------------------- mclust ------------------------------------

crabsMclust <- MclustDA(X, labels, modelType = "EDDA")
summary(crabsMclust)
modelo <- crabsMclust$model$`1`$modelName

#----------------------------- Cross validation -------------------------------

N <- nrow(X)
k <- length(unique(labels))
v <- 1:N

# resamplingCV <- 0 # to perform K-fold (take 1<=K<=N) or leave one out (take K=N)
resamplingCV <- 1 # to perform CV(K,p)

K <- 10 # Take a suitable value of K for the method  
# In the paper, K=300 and p=0.8, 0.95 for CV, K=200 for L.O.O
p <- 0.8 # probability of knowing the label of an observation

N <- nrow(X)
d <- ncol(X)
k <- length(unique(labels))
v <- 1:N
particion <- list()
c0 <- rep(0,K)
c1 <- rep(0,K)
c2 <- rep(0,K)

cat("Partition: \n")
for(h in 1:K){
  cat(h,"-")
  
  # Partición para hacer k-fold o LOO:
  # particion[[h]] <- sample(v,floor(length(v)/(K+1-h)))
  # v <- setdiff(v,particion[[h]])
  
  # Particion creda mediante particiones aleatorias con probabilidad p
  particion[[h]] <- v[runif(N)>p]
  
  test <- X[particion[[h]],]
  labelsTest <- labels[particion[[h]]]
  train <- X[-particion[[h]],]
  labelsTrain <- labels[-particion[[h]]]
  
  c0[h] <- cvMclust(train,test,labelsTrain,labelsTest,modelo)
  c1[h] <- cvGCPC(train,test,labelsTrain,labelsTest,par0=crabs2CPC$par,csh=100000,cvol=100000,niter=200,tol=10^(-7))
  c2[h] <- cvGPROP(train,test,labelsTrain,labelsTest,par0=crabs2PROP$par,csh=100000,cvol=100000,niter=200,tol=10^(-7))
}

cat("CRABS Cross Validation results: \n",
    "Mclust: \t Mean error= \t",mean(c0[c0>=0]),"\t sd/sqrt(K)= \t",sd(c0[c0>=0])/sqrt(length(c0[c0>=0])),"\n",
    "2-CPC: \t Mean error = \t",mean(c1),"\t sd/sqrt(K)= \t",sd(c1)/sqrt(K),"\n",
    "2-PROP: \t Mean error= \t",mean(c2),"\t sd/sqrt(K)= \t",sd(c2)/sqrt(K),"\n")


############################################################################################################################
#                                  OLIVEOIL DATASET- DISCRIMINANT ANALYSIS                                                 #  
############################################################################################################################

# THIS EXAMPLE NEEDS A LOT OF RANDOM INITIALIZATIONS, EXPENSIVE COMPUTATIONAL TIMES. HERE WE SET 
# NSTART=20, BUT MORE INITIALIZACIONS MIGHT BE NEEDED TO GET THE OPTIMAL SOLUTION.

library(pdfCluster)
data(oliveoil)
X <- as.matrix(oliveoil[,3:10])
labels <- as.numeric(oliveoil[,2])

#------------------------------- 2-CPC ------------------------------------

oliveoil2CPC <- GCPCclas(X,labels, M=2,csh=10000,cvol=10000,niter=100,tol=10^(-5),nstart=20,graph=F)
grafElipsesError(X,oliveoil2CPC$posterior,oliveoil2CPC$par,d=8,0,labels)
title(paste("2-CPC DISCRIMINANT ANALYSIS \n"," Loglik =", round(oliveoil2CPC$loglik,4), "\n df = ",
            oliveoil2CPC$df, "\n BIC=",round(oliveoil2CPC$BIC,4)),line=-7,outer=T)
cat("Region separation = ", c(rep("S",4),rep("I",2),rep("N",3)), "\t (S=south, I=Island of Sardinia, N=Centre-North)\n",
    "2-CPC separation = ", oliveoil2CPC$par$P)
error(labels, apply(oliveoil2CPC$posterior,1,which.max))

#------------------------------- 3-CPC ------------------------------------

oliveoil3CPC <- GCPCclas(X,labels, M=3,csh=10000,cvol=10000,niter=100,tol=10^(-5),nstart=20,graph=F)
grafElipsesError(X,oliveoil3CPC$posterior,oliveoil3CPC$par,d=8,0,labels)
title(paste("3-CPC DISCRIMINANT ANALYSIS \n"," Loglik =", round(oliveoil3CPC$loglik,4), "\n df = ",
            oliveoil3CPC$df, "\n BIC=",round(oliveoil3CPC$BIC,4)),line=-7,outer=T)
cat("Region separation = ", c(rep("S",4),rep("I",2),rep("N",3)), "\t (S=south, I=Island of Sardinia, N=Centre-North)\n",
    "2-CPC separation = ", oliveoil3CPC$par$P)
error(labels, apply(oliveoil3CPC$posterior,1,which.max))

#------------------------------- 3-PROP ------------------------------------

oliveoil3PROP <- GPROPclas(X,labels, M=3,csh=10000,cvol=10000,niter=100,tol=10^(-5),nstart=20,graph=F)
grafElipsesError(X,oliveoil3PROP$posterior,oliveoil3PROP$par,d=8,0,labels)
title(paste("3-PROP DISCRIMINANT ANALYSIS \n"," Loglik =", round(oliveoil3PROP$loglik,4), "\n df = ",
            oliveoil3PROP$df, "\n BIC=",round(oliveoil3PROP$BIC,4)),line=-7,outer=T)
cat("Region separation = ", c(rep("S",4),rep("I",2),rep("N",3)), "\t (S=south, I=Island of Sardinia, N=Centre-North)\n",
    "2-CPC separation = ", oliveoil3PROP$par$P)
error(labels, apply(oliveoil3PROP$posterior,1,which.max))

#------------------------------- mclust ------------------------------------

oliveMclust <- MclustDA(X, labels, modelType = "EDDA")
summary(oliveMclust)
modelo <- oliveMclust$model$`1`$modelName

#----------------------------- Cross validation -------------------------------

N <- nrow(X)
k <- length(unique(labels))
v <- 1:N

# resamplingCV <- 0 # to perform K-fold (take 1<=K<=N) or leave one out (take K=N)
resamplingCV <- 1 # to perform CV(K,p)

K <- 10 # Take a suitable value of K for the method  
# In the paper, K=300 and p=0.8, 0.95 for CV, K=572 for L.O.O
p <- 0.8 # probability of knowing the label of an observation


N <- nrow(X)
d <- ncol(X)
k <- length(unique(labels))
v <- 1:N
particion <- list()
o0 <- rep(0,K)
o1 <- rep(0,K)
o2 <- rep(0,K)
o3 <- rep(0,K)

cat("Partition: \n")
for(h in 1:K){
  cat(h,"-")
  
  # Partición para hacer k-fold o LOO:
  # particion[[h]] <- sample(v,floor(length(v)/(K+1-h)))
  # v <- setdiff(v,particion[[h]])
  
  # Particion creda mediante particiones aleatorias con probabilidad p
  particion[[h]] <- v[runif(N)>p]
  
  test <- X[particion[[h]],]
  labelsTest <- labels[particion[[h]]]
  train <- X[-particion[[h]],]
  labelsTrain <- labels[-particion[[h]]]
  
  o0[h] <- cvMclust(train,test,labelsTrain,labelsTest,modelo)
  o1[h] <- cvGCPC(train,test,labelsTrain,labelsTest,par0=oliveoil2CPC$par,csh=10000,cvol=10000,niter=100,tol=10^(-6))
  o2[h] <- cvGCPC(train,test,labelsTrain,labelsTest,par0=oliveoil3CPC$par,csh=10000,cvol=10000,niter=100,tol=10^(-6))
  o3[h] <- cvGPROP(train,test,labelsTrain,labelsTest,par0=oliveoil3PROP$par,csh=10000,cvol=10000,niter=100,tol=10^(-6))
}

cat("\n OLIVE OIL Cross Validation Results: \n",
    "Mclust: \t Media= \t",mean(o0[c0>=0]),"\t sd/sqrt(K)= \t",sd(o0[o0>=0])/sqrt(length(o0[o0>=0])),"\n",
    "2-CPC: \t Media= \t",mean(o1),"\t sd/sqrt(K)= \t",sd(o1)/sqrt(K),"\n",
    "3-CPC: \t Media= \t",mean(o2),"\t sd/sqrt(K)= \t",sd(o2)/sqrt(K),"\n",
    "3-PROP: \t Media= \t",mean(o3),"\t sd/sqrt(K)= \t",sd(o3)/sqrt(K),"\n")


############################################################################################################################
#                                          CANCER DATASET- CLUSTER ANALYSIS                                                #  
############################################################################################################################

# THIS EXAMPLE NEEDS A LOT OF RANDOM INITIALIZATIONS, EXPENSIVE COMPUTATIONAL TIMES. HERE WE SET 
# NSTART=20, BUT MORE INITIALIZACIONS MIGHT BE NEEDED TO GET THE OPTIMAL SOLUTION.

dataCancer <- read.csv(file="PATH/data.csv",header=T) # change path
lab0 <- read.csv(file="PATH/labels.csv",header=T)     # change path
labels <- as.numeric(factor(lab0[,2]))

cancer <- as.matrix(dataCancer)
col <- ncol(cancer)
N <- nrow(cancer)
cancer <- matrix(as.numeric(cancer[,2:col]),nrow=N)

cuadrados <- apply(cancer^2,2,sum)
nulos <- which(cuadrados<10^(-5))
cancer <- cancer[,-nulos]

# Principal components
cp <- prcomp(cancer,center=T,scale=T)
acumulados <- cumsum(cp$sdev^2)/sum(cp$sdev^2)
acumulados[1:14]
X <- cp$x[,1:14]

#------------------------------- 3-CPC ------------------------------------

cancer3CPC <- GCPC(X,k=5,M=3,csh=1000,cvol=1000,niter=100,tol=10^(-5),nstart=20,graph=F)
grafElipsesError(X,cancer3CPC$posterior,cancer3CPC$par,d=8,0,labels)
title(paste("3-CPC CLUSTERING \n Loglik =", round(cancer3CPC$loglik,4), "\n df = ",cancer3CPC$df, "\n BIC=",round(cancer3CPC$BIC,4)),line=-7,outer=T)
error(labels, apply(cancer3CPC$posterior,1,which.max))

#------------------------------- mclust ------------------------------------

mod <- Mclust(X,G=5)
summary(mod)



