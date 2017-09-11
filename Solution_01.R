# Excercise 1

library (microbenchmark)


################ Linear regression part C - Inversion method and your method #################

####### inverse method
####### Assume weights wi are all 1. So I didn't put it on my code

inv_func <- function(X,y){
  X <- Matrix(X, sparse = FALSE)
  beta_h <- solve(t(X) %*% X) %*% t(X) %*% y     #beta = Xt * X * Xt * y
  return (beta_h)
}


####### LU decomposition (A=LU)
# AX=b

LU_decomposition <- function(X,y) {
  X <- Matrix(X, sparse = FALSE)          
  A = t(X) %*% X                           # Let A = t(X)*W*X, b = t(X)*W*y
  b = t(X) %*% y
  LU_deco <- lu(A)                         # Find A = LU   =>  LUX = b
  L <- expand(LU_deco)$L                   # L is lower triangular matrix in LU decomposition
  U <- expand(LU_deco)$U                   # U is lower triangular matrix in LU decomposition
  d = solve(L,b)                           # Let UX = d   =>   Ld = b, find d
  betaLU_h <- solve(U,d)                   # X = inv(U)*d   (our estimate beta) 
  return (betaLU_h)
}


####### cholesky decomposition (A=LLt)
# AX=b

cholesky_decomposition <- function(X,y) {
  X <- Matrix(X, sparse = FALSE)
  A = t(X) %*% X                           # Let A = t(X)*W*X,b = t(X)*W*y
  b = t(X) %*% y
  L <- chol(A)                             # L: Find cholesky decompostion of A=LLt
  Lt <- t(L)                               # Lt: t(L)  =>  L*Lt*X = b  =>  let z = Lt*X
  z = solve(L,b)                           # Lz=b  =>  z = inv(L)*b
  betach_h <- solve(Lt,z)                  # X = inv(Lt)*z  (our estimate beta) 
  return (betach_h)
}


####### QR decomposition (A=QR)
# AX=b

QR_decomposition <- function(X,y) {

  X <- Matrix(X, sparse = FALSE)
  A = t(X) %*% X                           # Let A = t(X)*W*X,b = t(X)*W*y
  b = t(X) %*% y
  AA <- qr(A)                              # Find A = QR   =>  QRX = b
  Q <-qr.Q(AA)                             # Q part
  R <-qr.R(AA)                             # R part
  s = solve(Q, b)                          # Let s=RX => Qs=b, find s
  betaQR_h <- solve(R, s)                  # X = inv(R)*s  (our estimate beta) 
  return(betaQR_h)
}


nn <- c(50,100,500,1000)                   # Test n=50/100/500/1000
result <- list()

for(i in 1:4){
  n=nn[i]
  p=n/5
  sp=0.05
  X = matrix(rnorm(n*p),nrow=n,ncol=p)    
  mask = matrix(rbinom(n*p,1,sp), nrow=n)  # Generate sparse matrix
  X_sparse <- mask * X
  y <- matrix(rnorm(n), nrow = n)
  W=diag(1,nrow=n)

  result[[i]] <- microbenchmark(           # Use benchmark function
    inv_func(X,y),
    LU_decomposition(X,y),
    cholesky_decomposition(X,y),
    QR_decomposition(X,y),
    unit='ms')

}
names(result) <- c("n=50,p=10", "n=100,p=20", "n=500,p=100", "n=1000,p=200")
result


#$`n=50,p=10`
#Unit: milliseconds
#                         expr      min       lq     mean   median       uq      max neval
#               inv_func(X, y) 1.120799 1.151077 2.065534 1.175197 1.318376 25.34525   100
#       LU_decomposition(X, y) 1.175197 1.208041 2.164867 1.229082 1.338903 19.18701   100
# cholesky_decomposition(X, y) 2.499217 2.583893 4.960993 2.645989 3.078605 47.73146   100
#       QR_decomposition(X, y) 3.709311 3.828371 5.296360 3.888927 4.249183 26.59332   100

#$`n=100,p=20`
#Unit: milliseconds
#                         expr      min       lq     mean   median       uq      max neval
#               inv_func(X, y) 1.233701 1.271676 1.662201 1.301954 1.429737 12.66441   100
#       LU_decomposition(X, y) 1.232674 1.270649 1.773244 1.295282 1.376880 20.07277   100
# cholesky_decomposition(X, y) 2.574143 2.641884 4.022724 2.695255 2.910793 71.09273   100
#       QR_decomposition(X, y) 3.834529 3.947429 4.894033 4.053146 4.317437 47.67809   100

#$`n=500,p=100`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq       max neval
#               inv_func(X, y) 10.979109 11.248018 12.185938 11.437897 12.675190  23.19910   100
#       LU_decomposition(X, y)  6.438436  6.681687  7.796235  6.882342  8.143754  27.05313   100
# cholesky_decomposition(X, y)  7.791195  8.031879  8.937910  8.316185  9.567332  18.28688   100
#       QR_decomposition(X, y) 12.235389 12.529445 16.021041 12.884056 14.870087 221.74678   100

#$`n=1000,p=200`
#Unit: milliseconds
#                         expr      min       lq     mean   median       uq       max neval
#               inv_func(X, y) 74.70556 77.29407 89.98629 81.21071 83.29784 270.25418   100
#       LU_decomposition(X, y) 38.67269 42.12182 49.89788 44.51789 46.53574 239.36960   100
# cholesky_decomposition(X, y) 39.65596 42.46463 45.34468 45.63458 47.71299  52.15205   100
#       QR_decomposition(X, y) 61.15026 65.70581 69.23378 69.41717 72.04418  79.87540   100




#########################  Linear regression part D - Sparsity ##############################

cholesky_sparse <- function(X, y){ 
  
  X <- Matrix(X, sparse = TRUE)
  A = t(X) %*% X
  b = t(X) %*% y
  L <- chol(A)
  Lt <- t(L)
  z = solve(L,b)
  betach_h <- solve(Lt,z)
  return (betach_h)
}


result.sp <- list()
sp=c(0.05,0.15,0.25,0.5)      # Test sparsity = .05/.15/.25/.5


for(i in 1:4){
  n=500
  p=n/5
  X = matrix(rnorm(n*p),nrow=n,ncol=p)
  mask = matrix(rbinom(n*p,1,sp[i]), nrow=n) 
  X_sparse <- mask * X
  y <- matrix(rnorm(n), nrow = n)
  
  W=diag(1,nrow=n)
  
  result.sp[[i]] <- microbenchmark(      # Use benchmark function 
    inv_func(X,y),
    cholesky_decomposition(X,y),
    cholesky_sparse(X_sparse, y),
    unit='ms')
  
}
names(result.sp) <- c("sp=0.05", "sp=0.15", "sp=0.25", "sp=0.5")
result.sp


#$`sp=0.05`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq      max neval
#               inv_func(X, y) 11.017084 11.213635 11.865238 11.414804 11.719637 17.51300   100
# cholesky_decomposition(X, y)  7.801459  8.052407  8.667943  8.276156  8.575344 13.03185   100
# cholesky_sparse(X_sparse, y)  1.913159  2.050693  2.449387  2.102525  2.175911 12.29287   100

#$`sp=0.15`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq      max neval
#               inv_func(X, y) 10.951397 11.187462 11.765731 11.359379 11.560548 17.55200   100
# cholesky_decomposition(X, y)  7.798379  8.041629  8.723378  8.307460  8.847846 13.56352   100
# cholesky_sparse(X_sparse, y)  2.925163  3.052433  3.312249  3.122225  3.206902 12.18407   100

#$`sp=0.25`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq        max neval
#               inv_func(X, y) 11.014005 11.211069 13.596851 11.346550 11.622131 202.993930   100
# cholesky_decomposition(X, y)  7.786063  8.047275  8.581071  8.258708  8.533776  12.269259   100
# cholesky_sparse(X_sparse, y)  3.990538  4.087016  4.236138  4.153730  4.246105   7.328301   100

#$`sp=0.5`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq      max neval
#               inv_func(X, y) 10.999636 11.309088 12.387066 11.812010 12.467862 20.57877   100
# cholesky_decomposition(X, y)  7.821986  8.121687  8.961835  8.463983  9.144981 16.61081   100
# cholesky_sparse(X_sparse, y)  7.781957  7.940532  8.189725  8.101672  8.223299 12.98669   100
