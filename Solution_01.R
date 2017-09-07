# Excercise 1

library (Matrix)
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



#########################  Linear regression part D - Sparsity ##############################

cholesky_sparse <- function(X, y){ 
  
  X <- Matrix(X, sparse = TRUE)
  A = crossprod(X)
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

