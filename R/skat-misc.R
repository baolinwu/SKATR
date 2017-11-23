### saddlepoint approx: modified from Lumley survey package.
## @export
saddle = function(x,lambda){
  d = max(lambda)
  lambda = lambda/d
  x = x/d
  k0 = function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0 = function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n = length(lambda)
  if (any(lambda < 0)) {
    lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
  } else if (x>sum(lambda)){
    lmin = -0.01
  } else {
    lmin = -length(lambda)/(2*x)
  }
  lmax = min(1/(2*lambda[lambda>0]))*0.99999
  hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
  w = sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta)))
  v = hatzeta*sqrt(kpprime0(hatzeta))
  if(abs(hatzeta)<1e-4){
    return(NA)
  } else{
    return( pnorm(w+log(v/w)/w, lower.tail=FALSE) )
  }
}
## @export
Sadd.pval = function(Q.all,lambda){
  sad = rep(1,length(Q.all))
  if(sum(Q.all>0)>0){
    sad[Q.all>0] = sapply(Q.all[Q.all>0],saddle,lambda=lambda)
  }
  id = which(is.na(sad))
  if(length(id)>0){
    sad[id] = Liu.pval(Q.all[id], lambda)
  }
  return(sad)
}
### modified Liu method from Lee SKAT-O paper
## @export
Liu.pval = function(Q, lambda){
  c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lambda^i) }
  muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
  s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
  } else {
    l = 1/s2;  a = sqrt(l);  d = 0
  }
  muX = l+d;  sigmaX = sqrt(2)*a
  
  Q.Norm = (Q - muQ)/sigmaQ
  Q.Norm1 = Q.Norm*sigmaX + muX
  pchisq(Q.Norm1, df = l,ncp=d, lower.tail=FALSE)
}
## @export
Liu.qval.mod = function(pval, lambda){
  c1 = rep(0,4)
  c1[1] = sum(lambda); c1[2] = sum(lambda^2)
  c1[3] = sum(lambda^3); c1[4] =sum(lambda^4)
  muQ = c1[1]; sigmaQ = sqrt(2 *c1[2])
  s1 = c1[3]/c1[2]^(3/2); s2 = c1[4]/c1[2]^2
  beta1= sqrt(8)*s1; beta2 = 12*s2; type1 = 0
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2)); d = s1 *a^3 - a^2; l = a^2 - 2*d
  } else {
    type1 = 1; l = 1/s2; a = sqrt(l); d = 0
  }
  muX = l+d; sigmaX = sqrt(2) *a
  df = l
  q.org = qchisq(pval,df=df,lower.tail=FALSE)
  (q.org - df)/sqrt(2*df)*sigmaQ + muQ
}
## @export
KAT.pval <- function(Q.all, lambda, acc=1e-9,lim=1e6){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp = davies(Q.all[i],lambda,acc=acc,lim=lim); pval[i] = tmp$Qq
    if((tmp$ifault>0)|(pval[i]<=0)|(pval[i]>=1)) pval[i] = Sadd.pval(Q.all[i],lambda)
  }
  return(pval)
}
## Reproduce the p-value from the SKAT R package
## @export
SKAT.pval <- function(Q.all, lambda){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp = davies(Q.all[i],lambda,acc=1e-6,lim=1e4); pval[i] = tmp$Qq
    if((tmp$ifault>0)|(pval[i]<=0)|(pval[i]>=1)) pval[i] = Liu.pval(Q.all[i],lambda)
  }
  return(pval)
}
#' Fit a null binomial logistic regression model
#'
#' Fit a null binomial model to be used for variant set association test
#' @param  D 0-1 disease outcome
#' @param  X covariates to be adjusted, setting X=NULL with no covariate
#' @keywords KAT.null
#' @export
KAT.null <- function(D,X){
  if(is.null(X)){
    X0 = 1
    gl0 = glm(D~1, family='binomial')
  } else{
    X0 = cbind(1,X)
    gl0 = glm(D~X, family='binomial')
  }
  pi0 = gl0$fitted; Yv = pi0*(1-pi0); llk0=logLik(gl0)
  Yh = sqrt(Yv)
  Ux = svd(Yh*X0,nv=0)$u*Yh
  return(list(U0=D-pi0,pi0=pi0,Yv=Yv,llk0=llk0,Ux=Ux,coef=gl0$coef, Y=D,X=X, mode='D') )
}
#' Fit a null linear regression model
#'
#' Fit a null linear model to be used for variant set association test
#' @param  Y continuous outcome
#' @param  X covariates to be adjusted, setting X=NULL with no covariate
#' @keywords KAT.cnull
#' @export
KAT.cnull <- function(Y,X){
  if(is.null(X)){
    U0 = Y - mean(Y)
    s2 = sum(U0^2)/(length(Y)-1)
    Ux = matrix(1/sqrt(length(Y)), length(Y),1)
  } else{
    X0 = cbind(1,X)
    a0 = svd(X0,nv=0)
    U0 = lm(Y~X)$res
    s2 = sum(U0^2)/(length(Y)-dim(X0)[2])
    Ux = a0$u
  }
  return(list(U0=U0,s2=s2,Ux=Ux, mode='C'))
}

