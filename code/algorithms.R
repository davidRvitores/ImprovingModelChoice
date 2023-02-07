### RUN THIS CODE TO WORK WITH G-CPC AND G-PROP

################### 1) LIBRARIES REQUIRED  #################################################################################

library(MASS)
library(tclust)
library(mclust)
library(ellipse)
library(mvtnorm)
library(fastmatrix)
library(pracma)
library(modeest)
library(RColorBrewer)

# Color palette used for the plots 

paleta <- brewer.pal(n=12,name="Paired")
paleta[6] <- paleta[10]
palette(paleta)


################## 2) AUXILIAR FUNCIONS ####################################################################################

# Function to create k covariance matrices given:

# Lambda <- list with the k sets of eigenvalues
# ejes <- list with the M common orthogonal matrices
# P <- vector of length k with the group assignation 

createSigma <- function(Lambda,ejes,P){
  k <- length(Lambda)
  Sigma <- list()
  for(i in 1:k){
    j <- P[i]
    Sigma[[i]] <- ejes[[j]]%*%diag(Lambda[[i]])%*%t(ejes[[j]])
  }
  Sigma
}

# ----------------------------------------------------------------------------

# Compute the loglikelihood 

# X <- data matrix
# par <- list containing the parameters of a model:
# par$pi <- vector of length k containing the weights
# par$mu <- matrix containing in the row i the mean of the group i
# par$Lambda <- list with the k sets of eigenvalues (without imposing prod(Lambda[[i]])=1),
#               we are including here the shapes and volumes
# par$ejes <- list with the M common orthogonal matrices
# par$P <- vector of length k with the group assignation 

loglikelihood <- function(X,par){
  N <- nrow(X)
  k <- length(par$Lambda)
  M <- length(par$P)
  
  Sigma <- createSigma(par$Lambda,par$ejes,par$P)
  
  suma <- array(NA,dim=c(N,k))
  for(i in 1:k){
    suma[,i] <- par$pi[i]*dmvnorm(X,mean= as.numeric(par$mu[i,]),sigma = Sigma[[i]], log=F)
  }
  loglik <- 0
  for(j in 1:N){
    loglik <- loglik + log(sum(suma[j,])) 
  }
  loglik
}

# ---------------------------------------------------------------------------

# Function to plot the result of a procedure applied to a data matrix X given by:

# par <- define as above
# posterior <- matrix with the posterior probabilities of observation i arising from group j
# d <- dimension of plot, take the first d components
# componentes <- axis taken for the representation. If componentes=0, canonical axis, and if 
#             componentes=1,...,G the axis par$ejes[[componentes]]

grafElipses <- function(X,posterior,par,d,componentes){
  
  palette(brewer.pal(n=12,name="Paired"))
  
  dim <- ncol(X)
  G <- apply(posterior,1,which.max)
  mu <- par$mu
  k <- ncol(posterior)
  if(componentes==0){
    w <- k+1
    par$ejes[[k+1]] <- diag(dim)
  }else{
    w <- componentes
  }
  
  Sn <- createSigma(par$Lambda,par$ejes,par$P)
  
  A <- matrix(rep(0,d^2),nrow=d)
  position <- 1
  for(j in 1:(d-1)){
    for(i in (j+1):d){
      A[i,j] <- position
      position <- position+1
    }
  }
  for(i in 1:d){
    A[i,i] <- position
    position <- position+1
  }
  
  
  layout(mat = A,
         heights = rep(1,4), # Heights of the two rows
         widths = rep(1,4)) # Widths of the two columns
  par(pty="m")
  par(mar=c(0.5,0.5,0.5,0.5))
  #layout.show(9)
  
  
  for(l in 1:(d-1)){
    for(m in (l+1):d){
      X1 <- as.matrix(X)%*%par$ejes[[w]][,c(l,m)]
      mu1 <- mu%*%par$ejes[[w]][,c(l,m)]
      plot(X1,col=(G),cex=3.5,pch=".",xaxt="n",yaxt="n",
           xlim=c(min(X1[,1])-0.5,max(X1[,1])+0.5),
           ylim=c(min(X1[,2])-0.5,max(X1[,2])+0.5))
      #points(X1[which(grupo!=species),],col="red",cex=1.5,pch="x")
      for(i in 1:k){
        lines(ellipse(t(par$ejes[[w]][,c(l,m)])%*%Sn[[i]]%*%par$ejes[[w]][,c(l,m)]
                      ,centre=mu1[i,], npoints=500,level=0.9),col=(i),lwd=1.4)
        points(mu1,pch=as.character(par$P),lwd=5,cex=1.7)
        
      }
    }
  }
  for(i in 1:d){
    plot(100,100,xlim=c(-1,1),ylim=c(-1,1),
         xaxt="n",yaxt="n")
    text(0,0, paste("Var",i),cex=2)
    #text(0,0, variables[i], cex=2.2)
  }
  
  title("G-CPC/G-PROP PLOT",outer=T,line=-2)
  #title(paste("Axis class=",componentes),outer=T,line=-4)
  
}

# ----------------------------------------------------------------------------------------

# Same function as before, including the parameter labels in order to mark the missclassified observations

grafElipsesError <- function(X,posterior,par,d,componentes,labels){
  
  paleta <- brewer.pal(n=12,name="Paired")
  paleta[6] <- paleta[10]
  palette(paleta)
  
  dim <- ncol(X)
  G <- apply(posterior,1,which.max)
  mu <- par$mu
  k <- ncol(posterior)
  if(componentes==0){
    w <- k+1
    par$ejes[[k+1]] <- diag(dim)
  }else{
    w <- componentes
  }
  
  Sn <- createSigma(par$Lambda,par$ejes,par$P)
  
  A <- matrix(rep(0,d^2),nrow=d)
  position <- 1
  for(j in 1:(d-1)){
    for(i in (j+1):d){
      A[i,j] <- position
      position <- position+1
    }
  }
  for(i in 1:d){
    A[i,i] <- position
    position <- position+1
  }
  
  
  layout(mat = A,
         heights = rep(1,4), # Heights of the two rows
         widths = rep(1,4)) # Widths of the two columns
  par(pty="m")
  par(mar=c(0.5,0.5,0.5,0.5))
  #layout.show(9)
  
  classification <- apply(posterior,1,which.max)
  N <- length(labels)
  k <- length(unique(labels))
  grupo <- rep(0,N)
  niv <- character()
  for(i in 1:k){
    grupo[which(classification==i)] <- mfv(labels[which(classification==i)])
  }
  mal <- which(grupo!=labels)
  
  
  for(l in 1:(d-1)){
    for(m in (l+1):d){
      X1 <- as.matrix(X)%*%par$ejes[[w]][,c(l,m)]
      mu1 <- mu%*%par$ejes[[w]][,c(l,m)]
      plot(X1,col=(G),cex=1,pch=(labels+2),xaxt="n",yaxt="n",
           xlim=c(min(X1[,1])-0.5,max(X1[,1])+0.5),
           ylim=c(min(X1[,2])-0.5,max(X1[,2])+0.5))
      for(i in 1:k){
        lines(ellipse(t(par$ejes[[w]][,c(l,m)])%*%Sn[[i]]%*%par$ejes[[w]][,c(l,m)]
                      ,centre=mu1[i,], npoints=500,level=0.9),col=(i),lwd=1.4)
        points(mu1,pch=as.character(par$P),lwd=5,cex=2.5)
      }
      #points(X1[mal,],col=G[mal],cex=1,pch=(labels[mal]+2))
      points(X1[mal,],col="red",pch=1,cex=3,lwd=1)
    }
  }
  for(i in 1:d){
    plot(100,100,xlim=c(-1,1),ylim=c(-1,1),
         xaxt="n",yaxt="n")
    text(0,0, paste("Var",i),cex=2)
    #text(0,0, variables[i], cex=2.2)
  }
  
  #title("G-CPC/G-PROP PLOT",outer=T,line=-2)
  #title(paste("Axis class=",componentes),outer=T,line=-4)
  
}






# ----------------------------------------------------------------------------

# Classification error fot discriminant analysis, or clustering error in cluster analysis

# labels <- original labels of the data
# classification <- classification in groups created

error <- function(labels,classification){
  N <- length(labels)
  k <- length(unique(labels))
  grupo <- rep(0,N)
  niv <- character()
  for(i in 1:k){
    grupo[which(classification==i)] <- mfv(labels[which(classification==i)])
  }
  cat("Confusion Matrix: \n")
  print(table(labels,grupo))
  cat("\n error = ",sum(labels!=grupo),"/",N,"=",sum(labels!=grupo)/N,"\n")
}



# ----------------------------------------------------------------------------

# Posterior probabilities of the Wishart distribution, compute only the variable part of the density

# S <- sample covariance matrix
# Sigma <- list of k different matrices
# n <- number of observations used to compute S

post <- function(S,Sigma,n){
  M <- length(Sigma) 
  suma <- rep(0,M)
  for(k in 1:M){
    for(j in 1:M){
      suma[k] <- suma[k] + (det(Sigma[[j]])/det(Sigma[[k]])*
                              exp(sum(diag((inv(Sigma[[j]])-inv(Sigma[[k]]))%*%S))))^(-n/2)
    }
    suma[k] <- 1/suma[k]
  }
  return(suma) 
}


# ----------------------------------------------------------------------------

# Restricted Lambda values by the shape restriction 

# Lambda defined as above
# n <- vector with the computes values n_i
# c <- constant value >=1

LambdaRestr <- function(Lambda, n, c){
  zero.tol <- 10^(-8)
  ni.ini <- n
  K <- length(Lambda)
  p <- length(Lambda[[1]])
  
  autovalues <- array(NA,dim=c(p,K))
  for(i in 1:K){
    autovalues[,i] <- Lambda[[i]]
  }
  
  d <- t(autovalues)
  
  p = nrow (autovalues)
  K = ncol (autovalues)
  
  n <- sum(ni.ini)
  
  nis = matrix(data=ni.ini,nrow=K,ncol=p)
  
  
  d_=sort(c(d,d/c))
  dim=length(d_)
  d_1=d_
  d_1[dim+1]=d_[dim]*2
  d_2=c(0,d_)
  ed=(d_1+d_2)/2
  dim=dim+1;
  
  
  if ((max(d[nis>0]) <= zero.tol)){
    return (matrix (0, nrow = p, ncol = K))
  }
  
  
  if (max(d[nis>0])/min(d[nis>0])<=c){
    d[nis==0]=mean(d[nis>0])
    return (Lambda)
  }
  
  t <- s <- r <- array(0,c(K,dim))
  sol <- sal <- array(0,c(dim))
  
  for (mp_ in 1:dim)
  {
    for (i in 1:K)
    {  
      r[i,mp_]=sum((d[i,]<ed[mp_]))+sum((d[i,]>ed[mp_]*c))
      s[i,mp_]=sum(d[i,]*(d[i,]<ed[mp_]))
      t[i,mp_]=sum(d[i,]*(d[i,]>ed[mp_]*c))
    }
    
    sol[mp_]=sum(ni.ini/n*(s[,mp_]+t[,mp_]/c))/(sum(ni.ini/n*(r[,mp_])))
    
    e = sol[mp_]*(d<sol[mp_])+d*(d>=sol[mp_])*(d<=c*sol[mp_])+(c*sol[mp_])*(d>c*sol[mp_])
    o=-1/2*nis/n*(log(e)+d/e)
    
    sal[mp_]=sum(o)
  }
  
  eo=which.max(c(sal))
  m=sol[eo]
  
  
  L <- (m*(d<m)+d*(d>=m)*(d<=c*m)+(c*m)*(d>c*m))
  for(i in 1:K){
    Lambda[[i]] <- L[i,]
  }
  return(Lambda)
}


#  --------------------------------------------------------------------------------

# Best orthogonal matrix using the procedure of Browne and McNicholas

# Sn <- list with the k values of S_i computed
# n <- vector with the computes values n_i
# Lambda <- actual values of the eigenvalues
# B <- initilization orthogonal matrix
# tol <- tolerance to stop the iteration

ejesMCN <- function(Sn, n,Lambda,B,tol){
  
  k <- length(Sn)
  d <- nrow(Sn[[1]])
  
  invA <- list()
  for(i in 1:k){
    invA[[i]] <- n[i]*diag(Lambda[[i]]^(-1))
  }
  
  w <- rep(0,k)
  for(i in 1:k){
    w[i] <- eigen(Sn[[i]])$values[1]
  }
  
  iter <- 0
  diff <- 1
  
  while(iter<500 && diff > tol){
    B0 <- B
    iter <- iter+1
    
    Ft <- array(0,dim=c(d,d))
    for(i in 1:k){
      Ft <- Ft + invA[[i]]%*%t(B)%*%Sn[[i]]-w[i]*invA[[i]]%*%t(B)
    }
    decomposition <- svd(Ft)
    B <- decomposition$v%*%t(decomposition$u)
    # print(B)
    
    diff <- sum(abs(B%*%t(B0)))-sum(diag(abs(B%*%t(B0))))
  }
  B
}

### -----------------------------------------------------------------------------

# Function to apply the volume restriction
# LambdaR <- list with the best shape constrained eigenvalues (imposing prod(LambdaR[[i]])=1)
# Lambda <- list with the best eigenvalues without any restriction 

volRestrCPC <- function(LambdaR,Lambda,n,c){
  
  k <- length(Lambda)
  d <- length(Lambda[[1]])
  N <- sum(n)
  
  # optimal unrestricted volume
  
  vopt <- rep(0,k)
  for(j in 1:k){
    for(l in 1:d){
      vopt[j] <- vopt[j] + Lambda[[j]][l]/LambdaR[[j]][l]
    }
    vopt[j] <- 1/d * vopt[j]
  }
  j <- 3
  
  # possible values for the maximum mopt
  
  if(max(vopt)/min(vopt)<=c){
    return(vopt)
  }else{
    e <- sort(c(vopt,vopt/c))
    f <- c(e-10^(-16),e[2*k]+1)
    
    m <- rep(0,2*k+1)
    for(i in 1:(2*k+1)){
      menores <- as.numeric(vopt<f[i])
      mayores <- as.numeric(vopt>c*f[i])
      num <- 0
      den <- d*sum(n*(mayores+menores))
      for(j in 1:k){
        aux <- sum(Lambda[[j]]/LambdaR[[j]])
        num <- num + n[j]*aux*(menores[j]+mayores[j]/c)
      }
      m[i] <- num/den
    }
    
    sol <- rep(0,2*k+1)
    for(i in 1:(2*k+1)){
      vol <- rep(0,k)
      for(l in 1:k){
        vol[l] <- min(max(vopt[l],m[i]),c*m[i])
      }
      for(j in 1:k){
        for(l in 1:d){
          sol[i] <- sol[i] +n[j]*(log(vol[j]*LambdaR[[j]][l])+
                                    Lambda[[j]][l]/(vol[j]*LambdaR[[j]][l]))
        }
      }
    }
    J <- which.min(sol)
    for(l in 1:k){
      vol[l] <- min(max(vopt[l],m[J]),c*m[J])
    }
    vol
  }
}



####################################### 3) G-CPC ALGORITHM AND INTERNAL STEPS ####################################################

# In all the algorithms, the parameters are defined as above

# ----------------------------------------------------------------------------------------

# E- step in the EM algorithm

pasoE1 <- function(X,par){
  k <- length(par$Lambda)
  N <- nrow(X)
  Sigma <- createSigma(par$Lambda,par$ejes,par$P)
  
  logPosterior <- matrix(rep(0,N*k),nrow=N)
  result <- matrix(rep(0,N*k),nrow=N)
  for(i in 1:k){
    logPosterior[,i] <- log(par$pi[i])+dmvnorm(X,as.numeric(par$mu[i,]),Sigma[[i]], log=T)
  }
  posterior <- exp(logPosterior)
  for(i in 1:N){
    if(sum(posterior[i,])==0){
      result[i,which.max(logPosterior[i,])] <- 1
    }else{
      result[i,] <- abs(posterior[i,]/sum(posterior[i,]))
    }
  }
  result
}

# ----------------------------------------------------------------------------

# Compute the values of pi, mu, n, Sn in the M step of the EM algorithm

pasoM1 <- function(X,posterior){
  
  k <- ncol(posterior)
  N <- nrow(posterior)
  d <- ncol(X)
  
  Sn <- list()
  mu <- array(NA,dim=c(k,d))
  n <- rep(0,k)
  pi <- rep(0,k)
  
  for(i in 1:k){
    n[i] <- sum(posterior[,i])
    pi[i] <- n[i]/sum(N)
    Sn[[i]] <- matrix(rep(0,d^2),nrow=d)
    A <- cov.wt(X,posterior[,i],method="ML")
    mu[i,] <- A$center
    Sn[[i]] <- A$cov
  }
  res <- list()
  res$mu <- mu
  res$Sn <- Sn
  res$n <- n
  res$pi <- pi
  return(res)
}


# -------------------------------------------------------------------------------

# Compute the best P, it is the step P-V in the theory

pasoE2 <- function(Sn,n,par,csh){
  
  k <- length(Sn)
  M <- length(par$ejes)
  # distribucion posterior de las variables latentes
  posterior2 <- matrix(rep(0,M*k),nrow=k)
  for(i in 1:k){
    Sigma <- list()
    # Matrices de covarianzas estimadas
    for(j in 1:M){
      L <- list()
      L[[1]] <-diag(t(par$ejes[[j]])%*%Sn[[i]]%*%par$ejes[[j]])
      Lj <- diag(LambdaRestr(L,1,csh)[[1]])
      Sigma[[j]] <- par$ejes[[j]]%*%Lj%*%t(par$ejes[[j]])
      #Sigma[[j]] <- par$ejes[[j]]%*%diag(diag(t(par$ejes[[j]])%*%Sn[[i]]%*%par$ejes[[j]]))%*%t(par$ejes[[j]])
    }
    posterior2[i,] <- post(Sn[[i]],Sigma,n[[i]])
    par$P[i] <-  which.max(posterior2[i,])
  }
  for(j in 1:M){
    if(sum(par$P==j)==0){
      hecho <- 0
      while(hecho==0){
        l <- sample(1:k,1)
        if(sum(par$P==par$P[l])>1){
          par$P[l] <- j
          hecho <- 1
        }
      }
    }
  }
  par
}

# ------------------------------------------------------------------------------

# Compute the best Lambda (with the determinant and shape constraints) and axis given 
# the actual values of Sn, n.

# niter <- number of iterations
# tol <- tolerance to end the iterations
# csh, cvol <- constants for the constraints


pasoM2 <- function(Sn,n,par,csh,cvol,niter,tol){
  
  k <- length(Sn)
  M <- length(par$ejes)
  d <- nrow(Sn[[1]])
  N <- sum(n)
  
  P <- par$P
  B <- par$ejes
  LambdaR <- par$Lambda
  for(h in 1:k){
    LambdaR[[h]] <- LambdaR[[h]]/prod(LambdaR[[h]])^(1/d)
  }
  vol <- rep(0,k)
  
  loglik <- 10^8
  diff <- 1
  i <- 1
  
  
  while(i<niter && diff>tol){
    #cat(diff,"\n")
    B0 <- B
    LambdaR0 <- LambdaR
    vol0 <- vol
    
    Lambda <- list()
    
    for(h in 1:k){
      Lambda[[h]] <- diag(t(B[[P[h]]])%*%Sn[[h]]%*%B[[P[h]]])
    }
    vol <- volRestrCPC(LambdaR,Lambda,n,cvol)
    
    L <- list()
    for(h in 1:k){
      L[[1]] <- Lambda[[h]]
      LambdaR[[h]] <- LambdaRestr(L,1,csh)[[1]]
      LambdaR[[h]] <- LambdaR[[h]]/prod(LambdaR[[h]])^(1/d)
    }
    
    for(j in 1:M){
      r <- which(par$P==j)
      ni <- n[r]
      L <- list()
      LambdaVol <- list()
      for(h in 1:length(r)){
        L[[h]] <- Sn[[r[h]]]
        LambdaVol[[h]] <- vol[r[h]]*LambdaR[[r[h]]]
      }
      B[[j]] <- ejesMCN(L,ni,LambdaVol,B[[j]],tol)
    }
    
    diff <- sum((vol-vol0)^2)
    for(h in 1:k){
      diff <- diff + sum((LambdaR[[h]]-LambdaR0[[h]])^2)
    }
    for(h in 1:M){
      diff <- diff + norm(diag(d)-abs(t(B[[h]])%*%B[[h]]))
    }
    # cat(i,"-",loglikelihood(X,par),"\n")
    #loglik <- c(loglik,loglikelihood(X,par))
    #diff <- abs(loglik[i]-loglik[i-1])
  }
  par$ejes <- B
  for(h in 1:k){
    par$Lambda[[h]] <- vol[h]*LambdaR[[h]]
  }
  i <- i+1
  par
}


# -----------------------------------------------------------------------------

# Clustering C-CPC

# nstart <- number of random initializations
# graph <- if TRUE, displayis the graph of the solution of each random start

GCPC <- function(X,k,M,csh,cvol,niter,tol,nstart,graph){
  
  N <- nrow(X)
  d <- ncol(X)
  
  par_ <- list()
  posterior_ <- rep(0,M)
  
  #loglikF <- rep(0,nstart)
  loglikF <- rep(0,1)
  
  cat("start: ")
  for(rep in 1:nstart){
    cat(rep,"-")
    #cat("\n start",rep,"--------------------------------------------------------- \n\n")
    
    # Initial solution : ------------------------------------------
    
    par <- list()
    par$pi <- rep(1/k,k)
    par$mu <- X[sample(1:N,k,F),]
    par$P <- sample(1:M,k,rep=T)
    for(j in 1:M){
      if(sum(par$P==j)==0){
        hecho <- 0
        while(hecho==0){
          l <- sample(1:k,1)
          if(sum(par$P==par$P[l])>1){
            par$P[l] <- j
            hecho <- 1
          }
        }
      }
    }
    par$Lambda <- list()
    for(i in 1:k){
      par$Lambda[[i]] <- rep(1,d)
    }
    par$ejes <- list()
    for(i in 1:M){
      par$ejes[[i]] <- diag(d)
    }
    
    options(warn=-1)
    muestra <- sample(1:N,trunc(N/2))
    clus <- tclust(X[muestra,],k,alpha=0,nstart=200,restr.fact = min(csh,cvol))
    S <- cov(X)
    centros <- X[sample(1:N,k),]
    r <- ncol(clus$centers)
    centros[1:r,] <-  t(clus$centers) +0.1*mvrnorm(n=r,rep(0,d),S)
    par$mu <- centros
    #G <- clus$cluster
    #cat(clus$obj,"\n")
    options(warn=0)
    
    #sol0 <- GCPCclas(X,G,M,c,niter,tol,nstart=1,graph=F)
    #par <- sol0$par
    
    # ----------------------------------------------------------------------
    loglik <- 10000000
    i <- 1
    diff <- 10
    
    while(i<niter && diff>tol){
      i <- i+1
      
      #cat(" Iteration \t", i-1," \t loglik = :",round(loglik[i-1],2),"\n")
      
      posterior <- pasoE1(X,par)
      res <-  pasoM1(X,posterior)
      par$mu <- res$mu
      par$pi <- res$pi
      Sn <- res$Sn
      n <- res$n
      
      loglik2 <- 10000000
      h <- 1
      diff2 <- 10
      
      while(h<niter && diff2>tol){
        h <- h+1
        P0 <- par$P
        par <- pasoE2(Sn,n,par,csh)
        par <- pasoM2(Sn,n,par,csh,cvol,niter,tol)
        loglik2 <- c(loglik2,loglikelihood(X,par))
        diff2 <- sum((par$P-P0)^2)
        #diff2 <- abs(loglik2[h]-loglik2[h-1])
      }
      loglik <- c(loglik,loglikelihood(X,par))
      diff <- abs(loglik[i]-loglik[i-1])
    }
    loglikF[rep] <- loglikelihood(X,par)
    
    if(graph==T){
      grafElipses(X,posterior,par,min(d,8),componentes=0)
      title(paste("Loglikelihhod=",round(loglikF[rep],2)),line=-6,outer=T)
      title(paste("rep=",rep),line=-8,outer=T)
    }
    
    
    if(rep==1){
      posterior_ <- posterior
      par_ <- par
      loglik_ <- loglikF[rep]
    }else{
      if(loglikF[rep]>loglik_){
        posterior_ <- posterior
        par_ <- par
        loglik_ <- loglikF[rep]
      }
    }
  }
  cat("\n")
  res <- list( "par"= par_,
               "posterior"= posterior_,
               "loglik"=loglik_,
               "df"=k-1 + k*d+ (k*d+M*(d^2-d)/2),
               "BIC" = 2*loglik_- log(N)*(k-1 + k*d+ (k*d+M*(d^2-d)/2)))
  res
}

# ------------------------------------------------------------------------------

# Discriminant analysis G-CPC

# G <- vector of length N with the classification of the observations
# M <- number of classes

GCPCclas <- function(X,labels,M,csh,cvol,niter,tol,nstart,graph){
  
  N <- nrow(X)
  d <- ncol(X)
  k <- length(unique(labels))
  
  Sn <- list()
  n <- rep(0,k)
  mu <- array(NA,c(k,d))
  for(i in unique(labels)){
    ind <- which(labels==i)
    n[i] <- length(ind)
    Sn[[i]] <- cov(X[ind,])
    mu[i,] <- apply(X[ind,],2,mean)
    
  }
  
  par_ <- list()
  posterior_ <- rep(0,M)
  
  loglikF <- rep(0,nstart)
  
  cat("start: ")
  for(rep in 1:nstart){
    cat(rep,"-")
    #cat("start ", rep,":\t loglik =")
    
    # initial solution: ------------------------------------------
    
    par <- list()
    par$pi <- n/N
    par$mu <- mu
    par$Lambda <- list()
    for(i in 1:k){
      par$Lambda[[i]] <- eigen(Sn[[i]])$values
    }
    par$P <- sample(1:M,k,rep=T)
    
    for(j in 1:M){
      if(sum(par$P==j)==0){
        par$P[sample(1:k,1)] <- j
      }
    }
    par$ejes <- list()
    for(i in 1:M){
      par$ejes[[i]] <- eigen(Sn[[sample(which(par$P==i),1)]])$vectors
    }
    
    loglik <- 10000000
    h <- 1
    diff <- 10
    
    while(h<niter && diff>tol){
      P0 <- par$P
      h <- h+1
      par <- pasoE2(Sn,n,par,csh)
      par <- pasoM2(Sn,n,par,csh,cvol,niter,tol)
      loglik <- c(loglik,loglikelihood(X,par))
      #cat("--------------------------------------------------------- \n")
      diff <- sum((par$P-P0)^2)
      if(h<3){diff=1}
    }
    par <- pasoM2(Sn,n,par,csh,cvol,4*niter,tol*10^(-4))
    loglikF[rep] <- loglikelihood(X,par)
    #cat(round(loglikF[rep],2),"\n")
    
    posterior <- pasoE1(X,par)
    
    if(graph==T){
      grafElipses(X,posterior,par,min(d,8),componentes=0)
      title(paste("Loglikelihhod=",round(loglikF[rep],2)),line=-6,outer=T)
      title(paste("rep=",rep),line=-8,outer=T)
    }
    
    
    if(rep==1){
      posterior_ <- posterior
      par_ <- par
      loglik_ <- loglikF[rep]
    }else{
      if(loglikF[rep]>loglik_){
        posterior_ <- posterior
        par_ <- par
        loglik_ <- loglikF[rep]
      }
    }
  }
  cat("\n")
  res <- list( "par"= par_,
               "posterior"= posterior_,
               "loglik"=loglik_,
               "df"=  k*d+ (k*d+M*(d^2-d)/2),
               "BIC" = 2*loglik_- log(N)*( k*d+ (k*d+M*(d^2-d)/2)))
  res
}


################################ 4) G-PROP AND INTERNAL STEPS #####################################################################

# Shape constraint in the proportionality model

LambdaRestrPROP <- function(vol,Lambda, n, c){
  # Lambda autovalores estimados sin la restriccion de proporcionalidad (lista)
  
  k <- length(Lambda)
  d <- length(Lambda[[1]])
  N <- sum(n)
  
  # LambdaP autovalores estimados con la restriccion de proporcionalidad (vector)
  LambdaP <- rep(0,d)
  
  for(j in 1:k){
    LambdaP <- LambdaP + n[j]*Lambda[[j]]/(vol[j]*N)
  }
  
  #cat(sol(Sn,n,LambdaR0,vol,B),"\n")
  #cat(sol(Sn,n,LambdaR,vol,B),"\n")
  #cat(sol(Sn,n,LambdaP,vol,B),"\n")
  
  if(max(LambdaP)/min(LambdaP)<=c){
    return(LambdaP)
  }else{
    # cat("\n ??Restriccion autovalores!! \n")
    e <- sort(c(LambdaP,LambdaP/c))
    f <- c(e-10^(-5),e[2*d]+1)
    
    m <- rep(0,2*d+1)
    for(i in 1:(2*d+1)){
      menores <- as.numeric(LambdaP<f[i])
      #cat("menores:",menores,"\n")
      mayores <- as.numeric(LambdaP>c*f[i])
      #cat("mayores:",mayores,"\n")
      num <- 0
      den <- N*sum(menores+mayores)
      for(j in 1:k){
        num <- num + n[j]*(sum(Lambda[[j]]/vol[j]*menores)+ sum(Lambda[[j]]/vol[j]*mayores)/c)
      }
      m[i] <- num/den
    }
    
    sol <- rep(0,2*d+1)
    for(i in 1:(2*d+1)){
      LambdaP_ <- rep(0,d)
      for(l in 1:d){
        LambdaP_[l] <- min(max(LambdaP[l],m[i]),c*m[i])
      }
      for(j in 1:k){
        for(l in 1:d){
          sol[i] <- sol[i] +n[j]*(log(LambdaP_[l]*vol[j])+ Lambda[[j]][l]/(LambdaP_[l]*vol[j]))
        }
      }
    }
    J <- which.min(sol)
    LambdaPr <- rep(0,d)
    for(l in 1:d){
      LambdaPr[l] <- min(max(LambdaP[l],m[J]),c*m[J])
    }
    return(LambdaPr)
  }
}


# ------------------------------------------------------------------------------------------------

# Compute the best P, it is the step P-V in the theory 

pasoE2P <- function(Sn,n,par){
  
  k <- length(Sn)
  M <- length(par$ejes)
  d <- nrow(Sn[[1]])
  
  # M matrices de proporcionalidad que estamos considerando
  Sprop <- list()
  for(i in 1:M){
    j <- which(par$P==i)[1]
    Sprop[[i]] <- par$ejes[[i]]%*%diag(par$Lambda[[j]])%*%
      t(par$ejes[[i]])/prod(par$Lambda[[j]])^(1/d)
  }
  
  
  # distribucion posterior de las variables latentes
  posterior2 <- matrix(rep(0,M*k),nrow=k)
  
  for(i in 1:k){
    Sigma <- list()
    # Matrices de covarianzas estimadas
    for(j in 1:M){
      Sigma[[j]] <- sum(diag(inv(Sprop[[j]])%*%Sn[[i]]))/d * Sprop[[j]]
    }
    posterior2[i,] <- post(Sn[[i]],Sigma,n[i])
    par$P[i] <-  which.max(posterior2[i,])
  }
  for(j in 1:M){
    if(sum(par$P==j)==0){
      hecho <- 0
      while(hecho==0){
        l <- sample(1:k,1)
        if(sum(par$P==par$P[l])>1){
          par$P[l] <- j
          hecho <- 1
        }
      }
    }
  }
  par
}

# ------------------------------------------------------------------------------

# Compute the best Lambda, ejes (V-C step)

pasoM2P <- function(Sn,n,par,csh,cvol,niter,tol){
  
  k <- length(Sn)
  M <- length(par$ejes)
  d <- nrow(Sn[[1]])
  N <- sum(n)
  
  P <- par$P
  B <- par$ejes
  LambdaR <- par$Lambda
  for(h in 1:k){
    LambdaR[[h]] <- LambdaR[[h]]/prod(LambdaR[[h]])^(1/d)
  }
  vol <- rep(0,k)
  
  diff <- 1
  i <- 1
  
  
  while(i<niter && diff>tol){
    #cat(diff,"\n")
    B0 <- B
    LambdaR0 <- LambdaR
    vol0 <- vol
    
    Lambda <- list()
    LambdaRR <- list()
    for(h in 1:k){
      Lambda[[h]] <- diag(t(B[[P[h]]])%*%Sn[[h]]%*%B[[P[h]]])
      LambdaRR[[h]] <- LambdaR[[P[h]]]
    }
    vol <- volRestrCPC(LambdaRR,Lambda,n,cvol)
    
    for(j in 1:M){
      r <- which(par$P==j)
      ni <- n[r]
      LambdaAux <- list()
      volAux <- rep(0,length(r))
      for(h in 1:length(r)){
        LambdaAux[[h]] <- Lambda[[r[h]]]
        volAux[h] <- vol[r[h]]
      }
      LambdaR[[j]] <- LambdaRestrPROP(volAux,LambdaAux,ni,csh)
      LambdaR[[j]] <- LambdaR[[j]]/prod(LambdaR[[j]])^(1/d)
      
      L <- list()
      LambdaVol <- list()
      for(h in 1:length(r)){
        L[[h]] <- Sn[[r[h]]]
        LambdaVol[[h]] <- vol[r[h]]*LambdaR[[h]]
      }
      B[[j]] <- ejesMCN(L,ni,LambdaVol,B[[j]],tol)
    }
    
    #cat(sol(Sn,n,LambdaR0,vol0,B0),"\n")
    #cat(sol(Sn,n,LambdaR0,vol,B0),"\n")
    #cat(sol(Sn,n,LambdaR,vol,B0),"\n")
    #cat(sol(Sn,n,LambdaR,vol,B),"\n\n")
    
    
    diff <- sum((vol-vol0)^2)
    for(h in 1:M){
      diff <- diff + sum((LambdaR[[h]]-LambdaR0[[h]])^2)
      diff <- diff + norm(diag(d)-abs(t(B[[h]])%*%B[[h]]))
    }
    i <- i+1
    par$ejes <- B
    for(h in 1:k){
      par$Lambda[[h]] <- vol[h]*LambdaR[[P[h]]]
    }
    #cat(i,"-",loglikelihood(X,par),"\n")
    #cat(diff,"\n")
    #cat(B,"\n",LambdaR[[1]],"\n",vol,"\n\n")
  }
  par$ejes <- B
  for(h in 1:k){
    par$Lambda[[h]] <- vol[h]*LambdaR[[P[h]]]
  }
  par
}


# -----------------------------------------------------------------------------

# Clustering G-PROP

GPROP <- function(X,k,M,csh,cvol,niter,tol,nstart,graph){
  
  N <- nrow(X)
  d <- ncol(X)
  
  par_ <- list()
  posterior_ <- rep(0,M)
  
  #loglikF <- rep(0,nstart)
  loglikF <- rep(0,1)
  
  cat("start: ")
  
  for(rep in 1:nstart){
    cat(rep,"- ")
    
    # solucion inicial: ------------------------------------------
    par <- list()
    par$pi <- rep(1/k,k)
    par$mu <- X[sample(1:N,k,F),]
    par$P <- sample(1:M,k,rep=T)
    for(j in 1:M){
      if(sum(par$P==j)==0){
        hecho <- 0
        while(hecho==0){
          l <- sample(1:k,1)
          if(sum(par$P==par$P[l])>1){
            par$P[sample(1:k,1)] <- j
            hecho <- 1
          }
        }
      }
    }
    par$Lambda <- list()
    for(i in 1:k){
      par$Lambda[[i]] <- rep(1,d)
    }
    par$ejes <- list()
    for(i in 1:M){
      par$ejes[[i]] <- diag(d)
    }
    
    options(warn=-1)
    muestra <- sample(1:N,trunc(N/2))
    clus <- tclust(X[muestra,],k,alpha=0,nstart=100,restr.fact = min(csh,cvol))
    S <- cov(X)
    centros <- X[sample(1:N,k),]
    centros[1:ncol(clus$centers),] <-  t(clus$centers) +0.1*mvrnorm(n=ncol(clus$centers),rep(0,d),S)
    par$mu <- centros
    options(warn=0)
    
    
    
    # ----------------------------------------------------------------------
    loglik <- 10000000
    i <- 1
    diff <- 10
    
    while(i<niter && diff>tol){
      i <- i+1
      
      #cat(" Iteration \t", i-1," \t loglik = :",round(loglik[i-1],2),"\n")
      
      posterior <- pasoE1(X,par)
      res <-  pasoM1(X,posterior)
      par$mu <- res$mu
      par$pi <- res$pi
      Sn <- res$Sn
      n <- res$n
      
      loglik2 <- 10000000
      h <- 1
      diff2 <- 10
      
      while(h<niter && diff2>tol){
        P0 <- par$P
        h <- h+1
        par <- pasoE2P(Sn,n,par)
        par <- pasoM2P(Sn,n,par,csh,cvol,niter,tol)
        loglik2 <- c(loglik2,loglikelihood(X,par))
        diff2 <- sum((par$P-P0)^2)
        #diff2 <- abs(loglik2[h]-loglik2[h-1])
      }
      loglik <- c(loglik,loglikelihood(X,par))
      diff <- abs(loglik[i]-loglik[i-1])
    }
    loglikF[rep] <- loglikelihood(X,par)
    #cat(loglikF[rep],"\n")
    if(graph==T){
      grafElipses(X,posterior,par,min(d,8),componentes=0)
      title(paste("Loglikelihhod=",round(loglikF[rep],2)),line=-6,outer=T)
      title(paste("rep=",rep),line=-8,outer=T)
    }
    
    
    if(rep==1){
      posterior_ <- posterior
      par_ <- par
      loglik_ <- loglikF[rep]
    }else{
      if(loglikF[rep]>loglik_){
        posterior_ <- posterior
        par_ <- par
        loglik_ <- loglikF[rep]
      }
    }
  }
  cat("\n")
  res <- list( "par"= par_,
               "posterior"= posterior_,
               "loglik"=loglik_,
               "df"=k-1 + k*d+ (k+M*(d-1)+M*(d^2-d)/2),
               "BIC" = 2*loglik_- log(N)*(k-1 + k*d+ (k+M*(d-1)+M*(d^2-d)/2)))
  res
}

# ------------------------------------------------------------------------------

# discriminant analysis G-PROP

GPROPclas <- function(X,labels,M,csh,cvol,niter,tol,nstart,graph){
  
  N <- nrow(X)
  d <- ncol(X)
  k <- length(unique(labels))
  
  Sn <- list()
  n <- rep(0,k)
  mu <- array(NA,c(k,d))
  for(i in unique(labels)){
    ind <- which(labels==i)
    n[i] <- length(ind)
    Sn[[i]] <- cov(X[ind,])
    mu[i,] <- apply(X[ind,],2,mean)
    
  }
  
  par_ <- list()
  posterior_ <- rep(0,M)
  
  loglikF <- rep(0,nstart)
  
  cat("start: ")
  for(rep in 1:nstart){
    cat(rep, "-")
    #cat("start ",rep,": \t Loglik = ")
    
    # solucion inicial: ------------------------------------------
    
    par <- list()
    par$pi <- n/N
    par$mu <- mu
    par$Lambda <- list()
    for(i in 1:k){
      par$Lambda[[i]] <- eigen(Sn[[i]])$values
    }
    par$P <- sample(1:M,k,rep=T)
    for(j in 1:M){
      if(sum(par$P==j)==0){
        par$P[sample(1:k,1)] <- j
      }
    }
    par$ejes <- list()
    for(i in 1:M){
      par$ejes[[i]] <- eigen(Sn[[sample(which(par$P==i),1)]])$vectors
    }
    
    
    loglik <- 10000000
    h <- 1
    diff <- 10
    
    while(h<niter && diff>tol){
      h <- h+1
      P0 <- par$P
      par <- pasoE2P(Sn,n,par)
      #cat(par$P,"\n")
      par <- pasoM2P(Sn,n,par,csh,cvol,niter,tol)
      loglik <- c(loglik,loglikelihood(X,par))
      diff <- sum((par$P-P0)^2)
      if(h<3){diff<-1}
    }
    par <- pasoM2P(Sn,n,par,csh,cvol,4*niter,tol*10^(-4))
    loglikF[rep] <- loglikelihood(X,par)
    #cat(round(loglikF[rep],2),"\n")
    
    posterior <- pasoE1(X,par)
    
    if(graph==T){
      grafElipses(X,posterior,par,min(d,8),componentes=0)
      title(paste("Loglikelihhod=",round(loglikF[rep],2)),line=-6,outer=T)
      title(paste("rep=",rep),line=-8,outer=T)
    }
    
    
    if(rep==1){
      posterior_ <- posterior
      par_ <- par
      loglik_ <- loglikF[rep]
    }else{
      if(loglikF[rep]>loglik_){
        posterior_ <- posterior
        par_ <- par
        loglik_ <- loglikF[rep]
      }
    }
  }
  cat("\n")
  res <- list( "par"= par_,
               "posterior"= posterior_,
               "loglik"=loglik_,
               "df"=  k*d+ (k+M*(d-1)+M*(d^2-d)/2),
               "BIC" = 2*loglik_ - log(N)*(k*d + k +M*(d-1)+M*(d^2-d)/2))
  res
}



# CROSS VALIDATION FUNCTIONS -------------------------------------------------------

# Given a partition of the observations in train and test, these functions compute the proportion of test
# observations missclassified by the model estimated with the train observations. For the models 2-CPC and
# 2-PROP, we are taking as initial solution the global solution, in order to decrease computational times.


cvGCPC <- function(train,test,labelsTrain,labelsTest,par0,csh,cvol,niter,tol){
  # par0 <- parameters estimated by GCPCclas / GPROPclas
  
  if(is.matrix(test)==F){
    test <- matrix(test,nrow=1)
  }
  
  k <- length(unique(labelsTrain))
  d <- ncol(train)
  
  Sn <- list()
  n <- rep(0,k)
  mu <- array(NA,c(k,d))
  for(i in 1:k){
    ind <- which(labelsTrain==i)
    n[i] <- length(ind)
    Sn[[i]] <- cov(train[ind,])
    mu[i,] <- apply(train[ind,],2,mean)
  }
  
  par <- par0 # We take the fixed value of P estimated in the model, and the initialization of
  # covariance matrix parameters given by the parameters estimated in the model
  par$pi <- n/sum(n)
  par$mu <- mu
  
  par <- pasoM2(Sn,n,par,csh,cvol,niter,tol)
  pred <- apply(pasoE1(as.matrix(test),par),1,which.max)
  mal <- sum(as.numeric(pred)!=as.numeric(labelsTest))/length(labelsTest)
  mal
}

cvGPROP <- function(train,test,labelsTrain,labelsTest,par0,csh,cvol,niter,tol){
  # par0 <- parameters estimated by GCPCclas / GPROPclas
  
  if(is.matrix(test)==F){
    test <- matrix(test,nrow=1)
  }
  
  k <- length(unique(labelsTrain))
  d <- ncol(train)
  
  Sn <- list()
  n <- rep(0,k)
  mu <- array(NA,c(k,d))
  for(i in 1:k){
    ind <- which(labelsTrain==i)
    n[i] <- length(ind)
    Sn[[i]] <- cov(train[ind,])
    mu[i,] <- apply(train[ind,],2,mean)
  }
  
  par <- par0 # We take the fixed value of P estimated in the model, and the initialization of
  # covariance matrix parameters given by the parameters estimated in the model
  par$pi <- n/sum(n)
  par$mu <- mu
  
  par <- pasoM2P(Sn,n,par,csh,cvol,niter,tol)
  pred <- apply(pasoE1(test,par),1,which.max)
  mal <- sum(as.numeric(pred)!=as.numeric(labelsTest))/length(labelsTest)
  mal
}

cvMclust <- function(train,test,labelsTrain,labelsTest,modelo){
  if(is.matrix(test)==F){
    test <- matrix(test,nrow=1)
  }
  sink("NUL")
  Clas <- MclustDA(train, labelsTrain ,modelType = "EDDA", modelName = modelo)
  sink()
  mal <- -1
  if(is.null(Clas)==F){
    pred <- predict.MclustDA(Clas,test)$classification
    mal <-  sum(as.numeric(pred)!=as.numeric(labelsTest))/length(labelsTest)
  }
  mal
}

