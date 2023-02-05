

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
###                                                 SIMULATION EXAMPLES                                                  ###
###                                                                                                                      ###
############################################################################################################################
############################################################################################################################

# We start loading the data represented in the simulations of the examples. The different objects are explained in the different examples

# save(SigmaMatrixClas,SnMatrixClas XmatrixClas,labelsMatrixClas,muCPC,SigmaCPC,Xcpc,labelsCPC,muPROP,SigmaPROP,Xprop,labelsPROP, 
# file= "simulationData.RData")
load("C:/Users/davii/OneDrive/Desktop/DOCTORADO/TRABAJOS/Parsimonious classification /github/simulationData.RData")


############################################################################################################################
#                       CLASSIFICATION OF COVARIANCE MATRICES (FIGURE 1)                                                   #
############################################################################################################################


# SigmaMatrixClas <- 100 covariance matrices randomly created following the procedure explained in the example
# XmatrixClas <- 200*100 observations randomly generated from the 100 different distributions
# SnMatrixClas <- sample covariance matrices
# labelsMatrixClas <- labels of the observations

M <- 4 # number of classses to search

solCPC <- GCPCclas(XmatrixClas,labelsMatrixClas,M,csh=100,cvol=100,niter=100,tol=10^(-5),nstart=5,graph=F) # Common principal components solution 
grupo <- solCPC$par$P      # Classification in classes obtained
ejes <- solCPC$par$ejes    # Estimation of the axis in each class

par(pty="s")
par(mfrow=c(2,M))
par(mar=rep(3,4))

# Plot of the CPC classification obtained for the sample covariance matrices S1,...,S100
for(i in 1:M){
  ind <- which(grupo==i)
  L <- length(ind)
  plot(ellipse(SnMatrixClas[[ind[1]]],centre=c(0,0)),col=i+1,type="l",xlim=c(-4,4),ylim=c(-4,4),
       xlab="",ylab="",lwd=0.9)
  if(L>1){
    for(l in 1:L){
      lines(ellipse(SnMatrixClas[[ind[l]]],centre=c(0,0),npoints=500),col=i,lwd=0.9)
    }
  }
  for(j in 1:2){
    abline(a=0,b=ejes[[i]][2,j]/ejes[[i]][1,j],col="black",lty=2,lwd=3)
  }
}


solPROP <- GPROPclas(XmatrixClas,labelsMatrixClas,M,csh=100,cvol=100,niter=100,tol=10^(-5),nstart=5,graph=F)  # Proportionality solution
grupo <- solPROP$par$P   # Classification in classes obtained
Sigma <- list()          # list containing the M estimated proportional matrix
for(i in 1:M){
  j <- which(solPROP$par$P==i)[1]
  Sigma[[i]] <- solPROP$par$ejes[[i]]%*%diag(solPROP$par$Lambda[[j]])%*%t(solPROP$par$ejes[[i]])
  Sigma[[i]] <- Sigma[[i]]/det(Sigma[[i]])^(1/2)
}


# Plot of the PROPORTIONALITY classification obtained for the sample covariance matrices S1,...,S100
for(i in 1:M){
  ind <- which(grupo==i)
  L <- length(ind)
  plot(ellipse(SnMatrixClas[[ind[1]]],centre=c(0,0)),col=i+4,type="l",xlim=c(-4,4),ylim=c(-4,4),
       xlab="",ylab="",lwd=0.9)
  if(L>1){
    for(l in 1:L){
      lines(ellipse(SnMatrixClas[[ind[l]]],centre=c(0,0),npoints=500),col=i+4,lwd=0.9)
    }
  }
  for(j in 1:2){
    lines(ellipse(Sigma[[i]],centre=c(0,0),npoints=500),col="black",lwd=2)
  }
}
#par(cex.main=1.5)
#title("4-CPC classification",outer=T,line=-1.5,lwd=5)
#title("4-PROP classification",outer=T,line=-17.6,lwd=5)


############################################################################################################################
#                   2-CPC SIMULATION (FIGURE 2 CLUSTERING, TABLE 4 DISCRIMINANT ANALYSIS)                                  #
############################################################################################################################


# muCPC <- means of the theoretical distributions
# SigmaCPC <- covariance matrices of the theoretical distributions
# Xcpc <- simulated data, 100 observations per group
# labelsCPC <- labels associated to the simulated data

N <- nrow(Xcpc)
k <- length(SigmaCPC)

# Plot of the original groups and distributions

par(pty="s",mfrow=c(1,3),mar=rep(3,4),cex.main=1.3)

plot(as.vector(Xcpc[,1]),as.vector(Xcpc[,2]),col=labelsCPC, xlim=c(min(Xcpc[,1])-1,max(Xcpc[,1]+1)),
     ylim=c(min(Xcpc[,2])-1,max(Xcpc[,2]+1)),
     xlab="",ylab="",pch=".",cex=3,main="Original groups and covariance matrices")

for(i in 1:k){
  lines(ellipse(SigmaCPC[[i]],centre=muCPC[i,],npoints=500),col=i,lwd=2)
}


#################################### CLUSTER ANALYSIS ###################################

#----------------------------- 2-CPC solution ------------------------------------

clusCPC <- GCPC(Xcpc,k=6,M=2,csh=100,cvol=100,niter=200,tol=10^(-8),nstart=5,graph=F)
clusCPC$BIC

# Confusion table
# createSigma <- function to compute the covariance matrices from the orientation, shapes and classes
Sigma <- createSigma(clusCPC$par$Lambda,clusCPC$par$ejes,clusCPC$par$P)
Gaux <- apply(clusCPC$posterior,1,which.max) # Clusters created 
G <- rep(0,N)                                # Clusters created, with different labeling order (closer to the origianl groups)
niv <- rep(0,k)
for(i in 1:k){
  niv[i] <- mfv(labelsCPC[which(Gaux==i)])
  G[which(Gaux==i)] <- mfv(labelsCPC[which(Gaux==i)])
}
table(labelsCPC,G)
sum(labelsCPC!=G)

# Plot
plot(as.vector(Xcpc[,1]),as.vector(Xcpc[,2]),col=G, xlim=c(min(Xcpc[,1])-1,max(Xcpc[,1]+1)),ylim=c(min(Xcpc[,2])-1,max(Xcpc[,2]+1)),
     xlab="",ylab="",pch=".",cex=3,main="2-CPC clustering")
for(i in 1:k){
  lines(ellipse(Sigma[[i]],centre=clusCPC$par$mu[i,],npoints=500),col=niv[i],lwd=2)
  points(clusCPC$par$mu,pch=as.character(clusCPC$par$P),lwd=5,cex=1.7)
}

#------------------------------- mclust solution ------------------------------------

mod <- Mclust(Xcpc,G=6)
summary(mod)

Sigma1 <- list() # list with the covariance matrices estimated by mclust
for(i in 1:k){
  Sigma1[[i]] <- mod$par$var$sigma[,,i]
}

# Confusion table
Gaux <- mod$classification   # Clusters created
G <- rep(0,N)                # Clusters created, with different labeling order (closer to the origianl groups)
niv <- rep(0,k)
for(i in 1:k){
  niv[i] <- mfv(labelsCPC[which(Gaux==i)])
  G[which(Gaux==i)] <- mfv(labelsCPC[which(Gaux==i)])
}
table(labelsCPC,G)
sum(labelsCPC!=G)

# Plot
plot(as.vector(Xcpc[,1]),as.vector(Xcpc[,2]),col=G, xlim=c(min(Xcpc[,1])-1,max(Xcpc[,1]+1)),ylim=c(min(Xcpc[,2])-1,max(Xcpc[,2]+1)),
     xlab="",ylab="",pch=".",cex=3,main="Mclust: VEV model")
for(i in 1:k){
  lines(ellipse(Sigma1[[i]],centre=mod$par$mean[,i],npoints=500),col=niv[i],lwd=2)
  #points(t(mod$par$mean),pch=as.character(solPROP2$par$P),lwd=5,cex=1.7)
}

#################################### DISCRIMINANT ANALYSIS ###################################

#----------------------------- 2-CPC solution ------------------------------------

discCPC <- GCPCclas(Xcpc,labelsCPC,M=2,csh=100,cvol=100,niter=100,tol=10^(-5),nstart=3,graph=F)
discCPC$BIC
# error(labels,classification)  <- function to directly compute the confusion matrix and classification error
error(labelsCPC,apply(discCPC$posterior,1,which.max))

#------------------------------- mclust solution ------------------------------------

mod2 <- MclustDA(Xcpc,labelsCPC, modelType = "EDDA")
summary(mod2)
modelo <- mod2$models$`1`$modelName

Gmclust <- predict.MclustDA(mod2,Xcpc)$classification
error(labelsCPC,Gmclust)

#-------------- Cross-validation analysis for both methods ----------------------------


# resamplingCV <- 0 # to perform K-fold (take 1<=K<=N) or leave one out (take K=N)
resamplingCV <- 1 # to perform CV(K,p)

K <- 300 # Take a suitable value of K for the method 
p <- 0.9 # probability of knowing the label of an observation


v <- 1:N
particion <- list()
m0 <- rep(0,K)
m1 <- rep(0,K)

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
  
  test <- Xcpc[particion[[h]],]
  labelsTest <- labelsCPC[particion[[h]]]
  train <- Xcpc[-particion[[h]],]
  labelsTrain <- labelsCPC[-particion[[h]]]
  
  m0[h] <- cvMclust(train,test,labelsTrain,labelsTest,modelo)
  m1[h] <- cvGCPC(train,test,labelsTrain,labelsTest,par0=discCPC$par,csh=100,cvol=100,niter=100,tol=10^(-6))
}

cat("Results: \n",
    "Mclust: \t Mean error= \t",mean(m0),"\t sd/sqrt(K)= \t",sd(m0)/sqrt(K),"\n",
    "2-CPC: \t Mean error= \t",mean(m1),"\t sd/sqrt(K)= \t",sd(m1)/sqrt(K),"\n")


############################################################################################################################
#                             2-CPC model validation (FIRST ROW IN TABLES 3 AND 6)                                         #
############################################################################################################################

# Now we will compute the proportion of times in which our model 2-CPC improves the BIC value obtained by the best 
# mclust model. This is done by repeating m times the experiments above, for each sample size n within each group

n <- 100 # In the paper, this simulation is done with n=50,100,200 and m=1000 and also N=25 for discriminant analysis
N <- k*n

m <- 10 # In the paper, we are taking m = 1000, it requires long computational times. 
# In addition, we are taking nstart=10, it should be higher to ensure optimality of the solution (also an astute initial
# solution might be considered, and repeat only one start, as in the cross validation techniques)

BestCLUST <- rep(0,m)
BestDA <- rep(0,m)

for(l in 1:m){
  cat(l,"-")
  
  # Data simulation
  for(i in 1:k){
    if(i==1){
      X <- mvrnorm(n,muCPC[i,],SigmaCPC[[i]])
      labels <- rep(i,n)
    }else{
      labels <- c(labels,rep(i,n))
      X <- rbind(X,mvrnorm(n,muCPC[i,],SigmaCPC[[i]]))
    }
  }
  
  # Clustering
  sink("NUL")  # Just to delete the console output, in order to print only the loop repetitions. 
  # If it gives any problem, it must be commented
  # If it gives any error and console stops displaying the output, compile sink() several times
  solCPC <- GCPC(X,k=6,M=2,csh=100,cvol=100,niter=200,tol=10^(-8),nstart=10,graph=F)
  sink()
  bicCPC <- solCPC$BIC
  sink("NUL")
  mod <- Mclust(X,G=6)
  sink()
  bicMCLUST <-  mod$bic
  
  if(bicCPC>=bicMCLUST){
    BestCLUST[l] <- 1
  }
  
  # Discriminant analysis
  sink("NUL")
  solCPC <- GCPCclas(X,labels,M=2,csh=100,cvol=100,niter=100,tol=10^(-5),nstart=10,graph=F)
  sink()
  bicCPC <- solCPC$BIC
  sink("NUL")
  mod2 <- MclustDA(X,labels, modelType = "EDDA")
  sink()
  bicMCLUST <- mod2$bic
  
  if(bicCPC>=bicMCLUST){
    BestDA[l] <- 1
  }
}

cat(" Proportion of times in which clustering 2-CPC improves the BIC value:\n",sum(BestCLUST)/m,"\n",
    "Proportion of times in which discriminant analysis 2-CPC improves the BIC value:\n",sum(BestDA)/m)



############################################################################################################################
#                                  2-PROP SIMULATION (FIGURE 3 AND TABLE 5)                                                #
############################################################################################################################


# muPROP <- means of the theoretical distributions
# SigmaPROP <- covariance matrices of the theoretical distributions
# Xprop <- simulated data, 100 observations per group
# labelsPROP <- labels associated to the simulated data

N <- nrow(Xprop)
k <- length(SigmaPROP)

# Plot of the original groups and distributions

par(pty="s",mfrow=c(1,3),mar=rep(3,4),cex.main=1.3)

plot(as.vector(Xprop[,1]),as.vector(Xprop[,2]),col=labelsPROP, xlim=c(min(Xprop[,1])-1,max(Xprop[,1]+1)),
     ylim=c(min(Xprop[,2])-1,max(Xprop[,2]+1)),
     xlab="",ylab="",pch=".",cex=3,main="Original groups and covariance matrices")

for(i in 1:k){
  lines(ellipse(SigmaPROP[[i]],centre=muPROP[i,],npoints=500),col=i,lwd=2)
}



#################################### CLUSTER ANALYSIS ###################################

#----------------------------- 2-PROP solution ------------------------------------

clusPROP <- GPROP(Xprop,k=6,M=2,csh=100,cvol=100,niter=200,tol=10^(-8),nstart=5,graph=F)
clusPROP$BIC

# Confusion table
Sigma <- createSigma(clusPROP$par$Lambda,clusPROP$par$ejes,clusPROP$par$P)
Gaux <- apply(clusPROP$posterior,1,which.max) # Clusters created 
G <- rep(0,N)                                # Clusters created, with different labeling order (closer to the origianl groups)
niv <- rep(0,k)
for(i in 1:k){
  niv[i] <- mfv(labelsPROP[which(Gaux==i)])
  G[which(Gaux==i)] <- mfv(labelsPROP[which(Gaux==i)])
}
table(labelsPROP,G)
sum(labelsPROP!=G)


# Plot
plot(as.vector(Xprop[,1]),as.vector(Xprop[,2]),col=G, xlim=c(min(Xprop[,1])-1,max(Xprop[,1]+1)),
     ylim=c(min(Xprop[,2])-1,max(Xprop[,2]+1)),
     xlab="",ylab="",pch=".",cex=3,main="2-PROP clustering")
for(i in 1:k){
  lines(ellipse(Sigma[[i]],centre=clusPROP$par$mu[i,],npoints=500),col=niv[i],lwd=2)
  points(clusPROP$par$mu,pch=as.character(clusPROP$par$P),lwd=5,cex=1.7)
}

#------------------------------- mclust solution ------------------------------------

mod <- Mclust(Xprop,G=6)
summary(mod)

Sigma1 <- list() # list with the covariance matrices estimated by mclust
for(i in 1:k){
  Sigma1[[i]] <- mod$par$var$sigma[,,i]
}

# Confusion table
Gaux <-mod$classification       # Clusters created 
G <- rep(0,N)                   # Clusters created, with different labeling order (closer to the origianl groups)
niv <- rep(0,k)
for(i in 1:k){
  niv[i] <- mfv(labelsPROP[which(Gaux==i)])
  G[which(Gaux==i)] <- mfv(labelsPROP[which(Gaux==i)])
}
table(labelsPROP,G)
sum(labelsPROP!=G)

# Plot
plot(as.vector(Xprop[,1]),as.vector(Xprop[,2]),col=G, xlim=c(min(Xprop[,1])-1,max(Xprop[,1]+1)),
     ylim=c(min(Xprop[,2])-1,max(Xprop[,2]+1)),
     xlab="",ylab="",pch=".",cex=3,main="Mclust: VVV model")
for(i in 1:k){
  lines(ellipse(Sigma1[[i]],centre=mod$par$mean[,i],npoints=500),col=niv[i],lwd=2)
  #points(t(mod$par$mean),pch=as.character(solPROP2$par$P),lwd=5,cex=1.7)
}


#################################### DISCRIMINANT ANALYSIS ###################################

#----------------------------- 2-PROP solution ------------------------------------

discPROP <- GPROPclas(Xprop,labelsPROP,M=2,csh=100,cvol=100,niter=200,tol=10^(-5),nstart=5,graph=F)
discPROP$BIC
error(labelsPROP,apply(discPROP$posterior,1,which.max))

#------------------------------- mclust solution ------------------------------------

mod2 <- MclustDA(Xprop,labelsPROP, modelType = "EDDA")
summary(mod2)
modelo <- mod2$models$`1`$modelName

Gmclust <- predict.MclustDA(mod2,Xprop)$classification
error(labelsPROP,Gmclust)

#-------------- Cross-validation analysis for both methods ----------------------------


# resamplingCV <- 0 # to perform K-fold (take 1<=K<=N) or leave one out (take K=N)
resamplingCV <- 1 # to perform CV(K,p)

K <- 300 # Take a suitable value of K for the method 
p <- 0.9 # probability of knowing the label of an observation


v <- 1:N
particion <- list()
m0 <- rep(0,K)
m1 <- rep(0,K)


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
  test <- Xprop[particion[[h]],]
  labelsTest <- labelsPROP[particion[[h]]]
  train <- Xprop[-particion[[h]],]
  labelsTrain <- labelsPROP[-particion[[h]]]
  
  m0[h] <- cvMclust(train,test,labelsTrain,labelsTest,modelo)
  m1[h] <- cvGPROP(train,test,labelsTrain,labelsTest,par0=discPROP$par,csh=100,cvol=100,niter=100,tol=10^(-6))
}

cat("Results: \n",
    "Mclust: \t Mean error= \t",mean(m0),"\t sd/sqrt(K)= \t",sd(m0)/sqrt(K),"\n",
    "2-PROP: \t Mean error= \t",mean(m1),"\t sd/sqrt(K)= \t",sd(m1)/sqrt(K),"\n")




############################################################################################################################
#                               2-PROP model validation (SECOND ROW IN TABLES 3 AND 6)                                     #      
############################################################################################################################

# Now we will compute the proportion of times in which our model 2-PROP improves the BIC value obtained by the best 
# mclust model. This is done by repeating m times the experiment above, for each sample size n within each group

n <- 100 # In the paper, this simulation is done with n=50,100,200 and m=1000 and also N=25 for discriminant analysis
N <- k*n

m <- 10 # In the paper, we are taking m=1000, it requires long computational times. 
# In addition, we are taking nstart=10, it should be higher to ensure optimality of the solution (also an astute initial
# solution might be considered, and repeat only one start, as in the cross validation techniques)

BestCLUST <- rep(0,m)
BestDA <- rep(0,m)


for(l in 1:m){
  cat(l,"-")
  
  # Data simulation
  for(i in 1:k){
    if(i==1){
      X <- mvrnorm(n,muPROP[i,],SigmaPROP[[i]])
      labels <- rep(i,n)
    }else{
      labels <- c(labels,rep(i,n))
      X <- rbind(X,mvrnorm(n,muPROP[i,],SigmaPROP[[i]]))
    }
  }
  
  # Clustering
  sink("NUL")  # Just to delete the console output, in order to print only the loop repetitions. 
  # If it gives any problem, it must be commented
  # If it gives any error and console stops displaying the output, compile sink() several times
  solPROP <- GPROP(X,k=6,M=2,csh=100,cvol=100,niter=200,tol=10^(-8),nstart=10,graph=F)
  sink()
  bicPROP <- solPROP$BIC
  sink("NUL")
  mod <- Mclust(X,G=6)
  sink()
  bicMCLUST <-  mod$bic
  
  if(bicPROP>=bicMCLUST){
    BestCLUST[l] <- 1
  }
  
  # Discriminant analysis
  sink("NUL")
  solPROP <- GPROPclas(X,labels,M=2,csh=100,cvol=100,niter=100,tol=10^(-5),nstart=10,graph=F)
  sink()
  bicPROP <- solPROP$BIC
  sink("NUL")
  mod2 <- MclustDA(X,labels, modelType = "EDDA")
  sink()
  bicMCLUST <- mod2$bic
  
  if(bicPROP>=bicMCLUST){
    BestDA[l] <- 1
  }
}

cat(" Proportion of times in which clustering 2-CPC improves the BIC value:\n",sum(BestCLUST)/m,"\n",
    "Proportion of times in which discriminant analysis 2-CPC improves the BIC value:\n",sum(BestDA)/m)

