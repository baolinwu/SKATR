####
#' Sequence kernel association test (SKAT) with linear kernel using variant test statistics
#'
#' Compute accurate SKAT (linear kernel) p-value based on marginal variant score statistics
#'
#' We compute the SKAT based on the variant test statistics (typically of much smaller dimension), which
#' leads to much more efficient computations. Davies' method is used to compute the tail probability
#' of 1-DF chi-square mixtures with more stringent convergence criteria (acc=1e-9,lim=1e6). When it
#' fails, we then switch to the saddlepoint approximation.
#'
#' @param  obj a fitted null model using KAT.null() or KAT.cnull()
#' @param  G genotype matrix, sample in rows, variant in columns
#' @param  W.beta Beta parameters for variant weights
#' @return more accurate SKAT p-value
#' @keywords SKAT
#' @export
#' @references
#' Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics, 89, 82-93.
#'
#' Wu,B., Guan,W., and Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. AHG, 80(2), 123-135.
SKATh <- function(obj,G, W.beta=c(1,25)){
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  tmp = t(obj$Ux)%*%G
  if(obj$mode=='C'){
    Gs = t(G)%*%G - t(tmp)%*%tmp
    Zs = colSums(obj$U0*G)/sqrt(obj$s2)
  } else{
    Gs = t(G*obj$Yv)%*%G - t(tmp)%*%tmp
    Zs = colSums(obj$U0*G)
  }
  R = t(Gs*W)*W
  Z = Zs*W
  lam = svd(R, nu=0,nv=0)$d
  KAT.pval(sum(Z^2), lam)
}

#' Optimal sequence kernel association test (SKAT-O) using marginal variant score statistics
#'
#' Efficient SKAT-O p-value calculation using marginal variant score statistics directly
#'
#' Efficiently compute the SKAT-O significance p-value based on the variant test statistics (typically of much smaller dimension).
#' The individual p-values of weighted SKAT and burden test are comptued more accurately using the SKATh.
#' To obtain more accruate results, the one-dimensional integration is computed based on the convolution of survival function of
#' 1-DF chi-square mixtures and 1-DF chi-square density function. For details, please see the Wu et. al (2016) reference.
#' 
#' @param  obj a fitted null model using KAT.null() or KAT.cnull()
#' @param  G genotype matrix, sample in rows, variant in columns
#' @param  W.beta Beta parameters for variant weights
#' @param  rho proportion weight assigned to burden test statistic
#' @return SKAT-O p-value
#' @keywords SKAT-O
#' @export
#' @references
#' Lee, S., Wu, M. C., and Lin, X. (2012) Optimal tests for rare variant effects in sequencing association studies. Biostatistics, 13, 762-775.
#' 
#' Wu,B., Guan,W., and Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. AHG, 80(2), 123-135.
SKATOh <- function(obj,G, W.beta=c(1,25), rho=c(0,0.1^2,0.2^2,0.3^2,0.4^2,0.5^2,0.5,1)){
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  tmp = t(obj$Ux)%*%G
  if(obj$mode=='C'){
    Gs = t(G)%*%G - t(tmp)%*%tmp
    Zs = colSums(obj$U0*G)/sqrt(obj$s2)
  } else{
    Gs = t(G*obj$Yv)%*%G - t(tmp)%*%tmp
    Zs = colSums(obj$U0*G)
  }
  R = t(Gs*W)*W;  Z = Zs*W
  ##
  K = length(rho); K1 = K
  Qs = sum(Z^2); Qb = sum(Z)^2; Qw = (1-rho)*Qs + rho*Qb
  pval = rep(0,K)
  Rs = rowSums(R); R1 = sum(Rs); R2 = sum(Rs^2); R3 = sum(Rs*colSums(R*Rs))
  RJ2 = outer(Rs,Rs,'+')/N
  ## min-pval
  if(rho[K]>=1){
    K1 = K-1
    pval[K] = pchisq(Qb/R1, 1, lower.tail=FALSE)
  }
  Lamk = vector('list', K1);  rho1 = rho[1:K1]
  tmp = sqrt(1-rho1+N*rho1) - sqrt(1-rho1)
  c1 = sqrt(1-rho1)*tmp;  c2 = tmp^2*R1/N^2
  for(k in 1:K1){
    mk = (1-rho[k])*R + c1[k]*RJ2 + c2[k]
    Lamk[[k]] = abs(eigen(mk, sym=TRUE, only.val=TRUE)$val)
    pval[k] = KAT.pval(Qw[k],Lamk[[k]])
  }
  minP = min(pval)
  qval = rep(0,K1)
  for(k in 1:K1) qval[k] = Liu.qval.mod(minP, Lamk[[k]])
  lam = abs(eigen(R-outer(Rs,Rs)/R1, sym=TRUE, only.val=TRUE)$val)
  tauk = (1-rho1)*R2/R1 + rho1*R1;  vp2 = 4*(R3/R1-R2^2/R1^2)
  MuQ = sum(lam);  VarQ = sum(lam^2)*2
  sd1 = sqrt(VarQ)/sqrt(VarQ+vp2)
  if(K1<K){
    q1 = qchisq(minP,1,lower=FALSE)
    T0 = minP
  } else{
    tmp = ( qval-(1-rho)*MuQ*(1-sd1)/sd1 )/tauk
    q1 = min(tmp)
    T0 = pchisq(q1,1,lower=FALSE)
  }
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    x = (eta1-MuQ)*sd1 + MuQ
    KAT.pval(x,lam)*dchisq(xpar,1)
  }
  p.value = try({ T0 + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=1e-25)$val }, silent=TRUE)
  prec = 1e-4
  while(class(p.value)=='try-error'){
    p.value = try({ T0 + integrate(katint, 0,q1, abs.tol=minP*prec)$val }, silent=TRUE)
    prec = prec*2
  }
  min(p.value, minP*K)
}

