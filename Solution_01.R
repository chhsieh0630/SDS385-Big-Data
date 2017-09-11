# Excercise 1

#library (Matrix)
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
  A = crossprod(X)                          # Let A = t(X)*W*X, b = t(X)*W*y
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
  A = crossprod(X)                           # Let A = t(X)*W*X,b = t(X)*W*y
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
  A = crossprod(X)                           # Let A = t(X)*W*X,b = t(X)*W*y
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
  y <- matrix(rnorm(n), nrow = n)


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
#                         expr      min       lq     mean   median       uq       max neval
#               inv_func(X, y) 1.161854 1.295282 1.538133 1.422040 1.577022  7.227716   100
#       LU_decomposition(X, y) 1.135169 1.221384 1.577289 1.358918 1.550849 10.931895   100
# cholesky_decomposition(X, y) 1.241911 1.335311 1.616650 1.468739 1.616024  9.817255   100
#       QR_decomposition(X, y) 3.773972 4.265605 5.046695 4.674615 5.163168 13.258683   100

#$`n=100,p=20`
#Unit: milliseconds
#                         expr      min       lq     mean   median       uq       max neval
#               inv_func(X, y) 1.275782 1.389709 1.779012 1.495425 1.845418  7.192819   100
#       LU_decomposition(X, y) 1.193672 1.327614 1.694182 1.425632 1.843879  4.190680   100
# cholesky_decomposition(X, y) 1.294257 1.446673 1.796789 1.531862 1.792047  4.697708   100
#       QR_decomposition(X, y) 3.822212 4.329754 5.252616 4.706432 5.586034 10.249358   100

#$`n=500,p=100`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq      max neval
#               inv_func(X, y) 11.101247 11.795075 14.099569 12.394989 14.607335 38.54953   100
#       LU_decomposition(X, y)  5.887273  6.164907  7.455787  6.523112  8.125792 14.59399   100
# cholesky_decomposition(X, y)  5.733318  6.097680  8.059017  6.351707  7.506377 60.57036   100
#       QR_decomposition(X, y) 11.611354 12.528418 16.175685 13.751854 16.446083 51.66965   100

#$`n=1000,p=200`
#Unit: milliseconds
#                         expr      min       lq     mean   median       uq      max neval
#               inv_func(X, y) 74.62858 79.76352 90.92948 83.56777 92.23908 213.0493   100
#       LU_decomposition(X, y) 36.64561 37.82235 47.80017 40.23483 45.00182 324.4816   100
# cholesky_decomposition(X, y) 35.67671 37.56575 42.94068 40.37237 43.78249 110.2827   100
#       QR_decomposition(X, y) 58.80807 61.17181 73.22935 64.98530 69.36329 398.5036   100






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
sp=c(0.05,0.15,0.25,0.4)      # Test sparsity = .05/.15/.25/.4


for(i in 1:4){
  n=500
  p=n/5
  X = matrix(rnorm(n*p),nrow=n,ncol=p)
  mask = matrix(rbinom(n*p,1,sp[i]), nrow=n) 
  X_sparse <- mask * X
  y <- matrix(rnorm(n), nrow = n)

  
  result.sp[[i]] <- microbenchmark(      # Use benchmark function 
    inv_func(X,y),
    cholesky_decomposition(X,y),
    cholesky_sparse(X_sparse, y),
    unit='ms')
  
}
names(result.sp) <- c("sp=0.05", "sp=0.15", "sp=0.25", "sp=0.4")
result.sp

#$`sp=0.05`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq      max neval
#               inv_func(X, y) 11.060192 11.403514 12.789517 11.763257 12.781419 32.54525   100
# cholesky_decomposition(X, y)  5.834929  5.978108  6.600993  6.206989  6.540047 12.46222   100
# cholesky_sparse(X_sparse, y)  1.953188  2.081997  2.745106  2.196438  2.345775 33.44025   100

#$`sp=0.15`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq       max neval
#               inv_func(X, y) 11.078667 11.316785 11.830731 11.491782 11.685253 15.648077   100
# cholesky_decomposition(X, y)  5.804138  5.931921  6.238991  6.013518  6.210581  9.038238   100
# cholesky_sparse(X_sparse, y)  2.943637  3.044221  3.159165  3.118121  3.222811  5.620417   100

#$`sp=0.25`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq       max neval
#               inv_func(X, y) 11.099194 11.287021 11.815079 11.389145 11.678581 18.552712   100
# cholesky_decomposition(X, y)  5.829797  5.942185  6.280939  6.014544  6.206476 13.366451   100
# cholesky_sparse(X_sparse, y)  3.975142  4.110623  4.355402  4.165534  4.309226  8.616398   100

#$`sp=0.4`
#Unit: milliseconds
#                         expr       min        lq      mean    median        uq       max neval
#               inv_func(X, y) 11.147433 11.300877 11.644414 11.429172 11.633421 16.539996   100
# cholesky_decomposition(X, y)  5.814402  5.938079  6.255249  5.988884  6.114615 13.356187   100
# cholesky_sparse(X_sparse, y)  5.991964  6.078692  6.254161  6.128984  6.258307  8.620504   100
