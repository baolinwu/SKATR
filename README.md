# SKATR
The pakcage provides some efficient and accurate algorithms to compute significance p-values for SKAT and SKAT-O.

# Reference
Wu,B., Guan,W., Pankow,J.S. (2016).On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. *AHG*, 80(2), 123-135.

# Sample R codes

```r
library(SKAT)
library(SKATR)
library(CompQuadForm)
X = cbind(rbinom(5e3,1,0.5), rnorm(5e3))
G = matrix(rbinom(5e3*10, 2,0.01), 5e3,10)
Yc = rnorm(5e3) + rowMeans(X) + G[,1]*0.5 - G[,2]*0.3
eta = X[,1]-0.5 + X[,2] + G[,1]*0.5 + G[,2]*0.7 + G[,3]*0.6
Yd = rbinom(5e3,1,1/(1+exp(-eta)))
### continuous outcome
obj1 = SKAT_Null_Model(Yc~X, out_type='C')
SKAT(G,obj1)$p.value
SKAT(G,obj1, method='optimal.adj', r.cor=c((0:5)^2/100,0.5,1))$p.value
objc = KAT.cnull(Yc,X)
SKATh(objc,G)
SKATOh(objc,G, rho=c((0:5)^2/100,0.5,1))
### binary outcome
obj2 = SKAT_Null_Model(Yd~X, out_type='D')
SKAT(G,obj2)$p.value
SKAT(G,obj2, method='optimal.adj', r.cor=c((0:5)^2/100,0.5,1))$p.value
objd = KAT.null(Yd,X)
SKATh(objd,G)
SKATOh(objd,G, rho=c((0:5)^2/100,0.5,1))
```

