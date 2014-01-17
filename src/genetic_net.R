###Genetic Networks Functions
###Author: M.K. Lau (glomus@gmail.com)
###Initiated: Aug2010 (Maybe?)
###Updated: 17Jan2014
###This is a set of functions for conducting 
###genetic distance based network modeling in R.

require('MASS',quiet=TRUE)
require('sna',quiet=TRUE)

genetic.net <- function(d='distnace matrix',n='number of observations',alpha=0.05,partial=TRUE,crit=0.3,fix.na=FALSE,fix.inf=FALSE){

                                        #Gower transformation 
  C_ij <- gower(d)
                                        #Solve for the inverse of the Gower matrix
  p_ij <- ginv(C_ij)
                                        #Convert to correlations
  r_ij <- p_ij*0
  for (i in 1:nrow(p_ij)){
    for (j in 1:ncol(p_ij)){
      r_ij[i,j] <- (-p_ij[i,j])/sqrt(p_ij[i,i]*p_ij[j,j])
    }
  }
                                        #change sign (see Fortuna et al. 2009)
  R <- r_ij*-1
                                        #make diagonal 1 (perfect correlation with self)
  diag(R) <- 1
                                        #calculate tau
  tau <- -n*log((1-R^2))
  if (fix.na == TRUE){tau[is.na(tau)] <- 0}
  if (fix.inf == TRUE){tau[tau == Inf|tau == -Inf] <- 0}
  diag(tau)=0
                                        #Test for significance
  chisq <- qchisq(alpha,1,lower.tail=FALSE)
                                        #Remove non-significant correlations
  R.ind <- R*0
  R.ind[tau >= chisq] <- R[tau >= chisq]
                                        #calculate the pearson correlation matrix
  if (partial == FALSE){
    R.ind <- cov2cor(abs(C_ij))
    R.ind[abs(R.ind) < abs(crit)] <- 0
  }
                                        #re-name and export
  rownames(R.ind) <- colnames(R.ind) <- rownames(d)
  return(R.ind)

}


allele.patch=function(x,missing=-1,method=c('random','mean','abundance')){
                                        #Takes the most abundanct allele from all populations and applies them to missing data
  method=method[1]
  if (method=='mean'){
    
    p.max=0
    
    for (i in 1:ncol(x)){
      y=x[x[,i]!=missing,i]
      p.max[i]=as.numeric(names(table(y))[table(y)==max(table(y))])[1]
    }
    
    for (i in 1:ncol(x)){
      x[x[,i]==missing,i]<-p.max[i]
    }
  }
  
  else if (method=='abundance'){
    for (i in 1:ncol(x)){
      y=x[x[,i]!=missing,i]
      z=sample(names(table(y)),1,prob=(table(y)/sum(table(y))))
      x[x[,i]==missing,i]<-as.numeric(z)
    }
  }
  
  else if (method=='random'){
    for (i in 1:ncol(x)){
      y=x[x[,i]!=missing,i]
      z=sample(names(table(y)),1)
      x[x[,i]==missing,i]<-as.numeric(z)
    }
  }
  
  return(x)
  
}


allele.split=function(x){
                                        #separate the alleles for all the loci into two matrixes
  x1=x[,(1:ncol(x))[!((1:ncol(x)) %% 2 == 0)]]
  x2=x[,(1:ncol(x))[!((1:ncol(x)) %% 2 == 1)]]
  colnames(x2)=colnames(x1)
  for (i in 1:ncol(x1)){
    x1[,i]=paste(colnames(x1)[i],x1[,i],sep='')
  }
  for (i in 1:ncol(x2)){
    x2[,i]=paste(colnames(x2)[i],x2[,i],sep='')
  }
  return(list(x1,x2))
}


vectorize=function(x){
                                        #create a vector for the allele names
  out=character(nrow(x)*ncol(x))
  h=1
  for (i in 1:ncol(x)){
    for (j in 1:nrow(x)){
      out[h]=x[j,i]
      h=h+1
    }
  }
  return(out)	
}


codify=function(x1,x2){
  C=array('A',dim=c(nrow(x1),length(alleles)))
                                        #create the codification matrix
  for (i in 1:length(alleles)){
    for (j in 1:nrow(x1)){
      x3=unlist(c(x1[j,],x2[j,]))
      C[j,i]=length(x3[x3==alleles[i]])
    }
  }
  C=matrix(as.numeric(C),nrow=nrow(C),ncol=ncol(C))
  return(C)
}

centroids=function(C='codification matrix',pop='corresponding population names'){
                                        #Obtain the centroids for populations by averaging the multivariate coding vectors for each population
  pop=factor(pop)
  C.=array(0,c(nlevels(pop),ncol(C)))
  
  for (i in 1:nlevels(pop)){
    if (class(C[pop==levels(pop)[i],])=='numeric'){C.[i,]=C[pop==levels(pop)[i],]}
    else{C.[i,]=apply(C[pop==levels(pop)[i],],2,mean)}
    
  }
  rownames(C.)=levels(pop)
  C.=na.omit(C.)
  return(C.)
}

gdist=function(C.='codification centroid matrix',pk='allelic frequency table',K='allelic richness'){
                                        #Genetic distance
                                        #dij^2 = (1/2)*sum((1/K*pk)*(yik-yjk)^2)
  dij=array(0,c(nrow(C.),nrow(C.)))
  
  for (i in 1:nrow(dij)){
    for (j in 1:ncol(dij)){
      dij[i,j]=sum((1/(K*pk))*(C.[i,]-C.[j,])^2)/2
    }
  }
  
  colnames(dij)=rownames(dij)=rownames(C.)	
  return(as.dist(dij))
  
}

gower <- function(d='distance matrix'){
  d <- as.matrix(d)
  di. <- matrix(rep(apply(d,1,mean),nrow(d)),nrow=nrow(d),ncol=ncol(d),byrow=FALSE)
  d.j <- matrix(rep(apply(d,2,mean),nrow(d)),nrow=nrow(d),ncol=ncol(d),byrow=TRUE)
  out <- 0.5 * (d - di. - d.j + mean(d))
  return(out)
}


d2cov <- function(d='distance matrix'){
  d_i.=apply(d,2,mean)
  d_.j=apply(d,1,mean)
  d..=mean(unlist(d))
  c_ij=d*0
  for (i in 1:nrow(d)){
    for (j in 1:ncol(d)){
      c_ij[i,j]=0.5*(d[i,j]-d_i.[i]-d_.j[j]-d..)
    }
  }
  return(c_ij)
}

d2pcor <- function(d){

  C_ij <- gower(d)
  p_ij <- ginv(C_ij)
  r_ij <- p_ij*0
  
  for (i in 1:nrow(p_ij)){
    for (j in 1:ncol(p_ij)){
      r_ij[i,j] <- (-p_ij[i,j])/sqrt(p_ij[i,i]*p_ij[j,j])
    }
  }

  R <- r_ij* - 1
  diag(R) <- 1

  return(R)
  
}


edge.value=function(x,n1,n2){
  return(x[rownames(x)==n1,colnames(x)==n2])
}
