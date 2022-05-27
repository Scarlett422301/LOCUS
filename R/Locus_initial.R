Locus_initial <- function(Y,q,V,rho=0.95,R = NULL,maxIter = 100)
{ 
  ICcorr = icaimax(t(Y),nc = q,center=F,maxit=maxIter)
  S_ini = matrix(0,ncol=dim(ICcorr$S)[1],nrow=q)
  theta_ini = list()
  for(i in 1:q)
  {
    Sl = Ltrinv( ICcorr$S[,i],V,F)
    Sl = Sl + diag( rep(mean(ICcorr$S[,i]),V ))
    eigenSl = eigen(Sl)
    orderEigen = order(abs(eigenSl$values),decreasing = T)
    if(is.null(R))
    {
      Rl = 2
      while( TRUE )
      {
        eigenset = orderEigen[1:Rl]
        imgeRL = eigenSl$vectors[,eigenset]%*% diag(eigenSl$values[eigenset])%*% t(eigenSl$vectors[,eigenset])
        # image( imgeRL ) 
        if(cor(Ltrans(imgeRL,F),ICcorr$S[,i]) > rho) break
        Rl = Rl + 1
      }
    }else
    {
      Rl = R[i]; eigenset = orderEigen[1:Rl]
    }
    theta_ini[[i]] = list()
    
    theta_ini[[i]]$lam_l = eigenSl$values[ eigenset ]
    if( theta_ini[[i]]$lam_l[1]<0 ){theta_ini[[i]]$lam_l = -1*theta_ini[[i]]$lam_l}    
    theta_ini[[i]]$X_l = matrix(0,ncol = V, nrow = Rl)
    for(j in 1:Rl)
    {
      theta_ini[[i]]$X_l[j,] = eigenSl$vectors[,eigenset[j]]
    }
    S_ini[i,] = Ltrans( t(theta_ini[[i]]$X_l)%*%diag(theta_ini[[i]]$lam_l)%*%theta_ini[[i]]$X_l,F)
  }
  
  A_ini = Y%*%t(S_ini)%*%solve(S_ini%*%t(S_ini))
  M_ini =  solve(t(A_ini)%*%A_ini)%*%t(A_ini) # Generalized inverse of A
  # Scale Up
  for(l in 1:q)
  {
    scaleL = sqrt(sum(A_ini[,l]^2))
    A_ini[,l] = A_ini[,l] / scaleL
    theta_ini[[l]]$X_l = theta_ini[[l]]$X_l * sqrt(scaleL)
    theta_ini[[l]]$M_l = M_ini[l,]
  }
  return(list(A=A_ini,theta = theta_ini,S = S_ini))
}

