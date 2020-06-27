# Initialize matrix of coefficients

library(mvtnorm)


n=100                                         #number of subjects
q=10                                         #number of basis functions
r=10                                         #number of covariates
k=2                                          #number of latent factors
rep=1
nrun=11000                                    #number of Gibbs iterations



knots<-seq(0,1,1/q)
bw=1

timevec=seq(0,1,.2)                         #assuming every subject is measured at same time points
ni=rep(length(timevec),n)
Bmat=matrix(1,length(timevec),q)

for(j in 2:q){
  Bmat[,j]=exp(-bw*(timevec-knots[j+1])^2)
}

tildX=matrix(rnorm(n*r),n,r)       #Data matrix
trueBeta=100*matrix(1,r,k)
trueLambda=diag(q)[,1:k]
ymat=Bmat%*%(trueLambda%*%(t(trueBeta)%*%t(tildX)))
yy=c(ymat)+c(rnorm(length(timevec)*n,0,1))

GramX=t(tildX)%*%tildX

rep=1
betaout=matrix(0,r*k,200)
lambdaout=matrix(0,q*k,200)

for(g in 1:rep){
  

  # --- Define hyperparameter values --- # 
  as = .5;bs = .25                          # gamma hyperparameters for diagonal elements of inv(Sigma)
  df = 5                                    # gamma hyperparameters for t_{ij}
  ad1 = 1.5;bd1 = 1                         # gamma hyperparameters for delta_ 1
  ad2 = 1.5;bd2 = 1                         # gamma hyperparameters delta_h, h >= 2
  adf = 1; bdf = 1                          # gamma hyperparameters for ad1 and ad2 or df

# --- Initial values --- #
  apsi = .5; bpsi = .2                      # gamma hyperparameters for psi
  Psiinv = rgamma(1,apsi, bpsi)
  sig = rgamma(q,as,bs)                  # diagonals of sigmainv
  Sigma = diag(1/sig)                      # Sigma (inv gamma diagonal elements)
  Sigmainv = diag(sig)
  Lambda = matrix(0,q,k)                       # loading matrix

# -- Initialize CovOmega -- #
  CovOmega = matrix(rgamma(r*k,.5, 2),r,k)

# -- Initialize Cauchy prior on Beta coefficients -- #
  Beta = matrix(0,r,k)
  
  for(l in 1:k){
    Beta[,l] = rmvnorm(1,rep(0,r), diag(1/(CovOmega[,l])))
  }

# -- Initialize latent factors -- #
  InvEta = diag(k)
  eta=matrix(0,n,k)
  for(j in 1:k){
   eta[,j]=rmvnorm(1,tildX%*%Beta[,j],diag(n))        # latent factors
  }
  THETA=matrix(0,n,q)
  for(i in 1:n){
  THETA[i,]= rmvnorm (1,c(Lambda%*%c(eta[i,])), Sigma)}     # Basis function coefficients
  tm = matrix(rgamma(q*k,df/2,df/2),q,k)              # local shrinkage coefficients
  delta =c(rgamma(1,ad1,1/bd1),rgamma(k-1,ad2,1/bd2)) # gobal shrinkage coefficients multilpliers
  tau = cumprod(delta)                                # global shrinkage coefficients
  Ptht = t(t(tm) * tau)

#####################################
# ----- Start Gibbs sampling ------ #
#####################################
  
  for(iter in 1:nrun){
    # -- Update Lambda -- #
    Lambda = matrix(0,q,k)
    for(j in 1:q){
    Vlam1 = diag(Ptht[j,]) + sig[j]*(t(eta)%*%eta)
    Vlam = solve(Vlam1)
    Vlam = (Vlam + t(Vlam))/2
    Elam = c(sig[j]*Vlam%*%t(eta)%*%(THETA[,j]))                                 
    Lambda[j,] = rmvnorm(1,Elam,Vlam)}
    
    #k = size(Lambda,2)
    
    # -- Update tij's -- #
    #for(i in 1:q){
    #  for(j in 1:k){
    #    tm[i,j]=rgamma(1,df/2+0.5, df/2+tau[j]*Lambda[i,j]^2/2)      # Local shrinkage coefficients
    #  }
    #}
    
    tm=matrix(rgamma(q*k,df/2+0.5,1),q,k)/(df/2+t(tau*t(Lambda^2))/2)
        
    # -- Update delta -- #
      ad = ad1 + q*k/2               
      bd = bd1 + 0.5* (1/delta[1])*sum(tau*rowSums(tm*Lambda^2))
      delta[1] = rgamma(1,ad,bd)
      tau = cumprod(delta)
      for(h in 2:k){
        ad = ad2 + q*(k-h+1)/2     
        temp1 = tau*rowSums(tm*Lambda^2)
        bd = bd2 + 0.5* (1/delta[h])*sum(temp1[h:k])
        delta[h] = rgamma(1,ad,bd)
        tau  = cumprod(delta)
      }
                                                 
        
      # -- Update precision parameters -- #
      Ptht = t(t(tm) * tau)

      
                                        
      # -- Update Sigma -- #         
      THETAtil = t(THETA) - Lambda%*%t(eta)
      sig = rgamma(q,as + n/2)/(bs+0.5*rowSums(THETAtil^2))  
      Sigma = diag(1/sig)                                     
      Sigmainv = diag(sig)
       
      # -- Update Psi -- #
      N=sum(ni[1:n])             #total count of observations, ni for ith subject
      Fvec = rep(0,N)
      for(h in 1:n){
        starth=ifelse(h==1,1,(sum(ni[1:(h-1)])+1))
        Fvec[starth:sum(ni[1:h])]= Bmat%*%THETA[h,]#[(sum(ni[1:(h-1)])+1):sum(ni[1:h]),]
      }
  
        Psiinv = rgamma(1,apsi + (N/2),(bpsi + 0.5*sum((yy - Fvec)^2)))                                                                 
        psiout = Psiinv
        
      # -- Update of Cov. Omega -- #
        #for(l in 1:k){
        #   for(g in 1:r){
        #       CovOmega[g,l] = rgamma(1, .5 + .5*(Beta[g,l])^2);
        #   }
        #}
        
        CovOmega=matrix(rgamma(k*r,1,1),r,k)/(.5+.5*Beta^2)
       
      # -- Update of Beta -- #
        #GramX=t(tildX)%*%tildX
        Beta = matrix(0,r,k)
        for(l in 1:k){
            CovCauchy = diag(CovOmega[,l])
            BetaPostCov = solve(CovCauchy + GramX)          
            BetaPostMean = c(BetaPostCov%*%(t(tildX)%*%eta[, l])) 
            Beta[,l] = rmvnorm(1,BetaPostMean, BetaPostCov) 
        }
                                                                    
        
        # -- Update of eta -- # 
        for(h in 1:n){
          Bh = Bmat#[(sum(ni[1:(h-1)])+1): sum(ni[1:h]),]
          inn=((1/Psiinv)*diag(ni[h]) + Bh%*%Sigma%*%t(Bh))
          invinner = solve(inn)
          invinner = (invinner + t(invinner))/2
          etavar = t(Lambda)%*%t(Bh)%*%(invinner)%*%Bh%*%Lambda + diag(k)
          etav = (etavar + t(etavar))/2
          invetavar = solve(etav)
          invetar = (invetavar + t(invetavar))/2
          starth=ifelse(h==1,1,(sum(ni[1:(h-1)])+1))
          meaneta = t(Beta)%*%(tildX[h,]) + t(Lambda)%*%t(Bh)%*%invinner%*%(yy[starth:sum(ni[1:h])])
          eta[h,] = rmvnorm(1,invetar%*%meaneta, invetar)
        }
   
      # -- Update Theta -- # 
        for(h in 1:n){
            #Bh = Bmat[(sum(ni[1:(h-1)])+1):sum(ni[1:h]),]
            covfun = Sigmainv + Psiinv*t(Bh)%*%Bh
            covfun = (covfun + t(covfun))/2
            Invcovfun = solve(covfun)
            Invcovfun =(Invcovfun + t(Invcovfun))/2
            starth=ifelse(h==1,1,(sum(ni[1:(h-1)])+1))
            meanvec = c(Invcovfun%*%(Sigmainv%*%(Lambda%*%eta[h,]))+ Psiinv*t(Bh)%*%(yy[starth:sum(ni[1:h])]))
            THETA[h,] = c(rmvnorm(1,meanvec, Invcovfun))
        }
       if(iter>10000 & iter%%5==0){
         betaout[,((iter-10000)/5)]=c(Beta)
         lambdaout[,((iter-10000)/5)]=c(Lambda)
       } 
  if(iter%%1000==0){
    print(iter)      
  }
  }
  
}
