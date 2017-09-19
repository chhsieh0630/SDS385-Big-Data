
################################## Generalized linear models ###################################

rm(list=ls())

data = read.csv("C:/R/bigdata/wdbc.csv", header=FALSE)

y = data[,2]
y = ifelse(y=='B',0,1)
y = as.matrix(y)
X = data[,c(3:12)]
X = scale(X)
X = as.matrix(cbind(1,X))


w <- function(X,beta){                        # w = 1 / [ 1 + e^(-Xb) ]
  weight = 1 / (1+ exp(-X %*% beta))
  return(weight)
}

ll <- function(X,y,beta) {                    # The MLE function
  logl <- -sum( y*log(w(X,beta)+1E-5) + (1-y)*log(1-w(X,beta)+1E-5)  )
  return(logl)                                # add 1E-5 is to prevent log0 => error
}

grad <- function(X,y,beta) {
  g <- -t(X) %*% ( y - w(X,beta)       )      # calculate the gradient = Xt*(y-w)
  return(g)
}


beta <- matrix(0,nrow=11,ncol=1)
mle <- ll(X,y,beta)
beta.log <-matrix(0,nrow=11,ncol=1)

alpha = 0.01            

for (i in 1:1000) {

  beta.new = beta - alpha * grad(X,y,beta)    # Beta_k+1 = Beta_k + alpha_k * grad
  ll(X,y,beta.new)                            # Find the new beta and calculate MLE again and again and again.....
  mle <- c(mle, ll(X,y,beta.new)) 
  beta.log <- cbind(beta, beta.new) 
  beta = beta.new

}
  

mle
plot(mle,log='xy', ylim=c(70,400), type='l', col='blue')


#[1] 394.28696 128.30137 106.93759  98.13894  93.51579  90.57664  88.35646  86.54833  85.03425
#[10]  83.75215  82.66070  81.72878  80.93170  80.24926  79.66462  79.16357  78.73392  78.36519
#[19]  78.04830  77.77541  77.53978  77.33566  77.15814  77.00309  76.86704  76.74710  76.64086
#[28]  76.54631  76.46179  76.38590  76.31750  76.25560  76.19940  76.14821  76.10144  76.05858
#[37]  76.01922  75.98296  75.94950  75.91853  75.88983  75.86315  75.83831  75.81513  75.79346
#[46]  75.77315  75.75409  75.73616  75.71927  75.70332  75.68823  75.67394  75.66037  75.64746



################################## Newton Method ####################################

hess <- function(X,beta) {                  # Find the Hessian Matrix to calculate by Newton method
                                            # Hessian Matrix = \/2 * f(x_k)
  
  H = crossprod(diag(as.vector(sqrt(w(X,beta) * (1 - w(X,beta))))) %*% X)
#  H = crossprod(diag(as.vector(    w(X,beta) * (1 - w(X,beta))    )) %*% X)
  return(H)
}

beta <- matrix(0,nrow=11,ncol=1)
nt <- ll(X,y,beta)
beta.log <-matrix(0,nrow=11,ncol=1)


for (i in 1:1000) {
  
  beta.new = beta - solve(hess(X,beta)) %*% grad(X,y,beta)   # Beta_k+1 = Beta_k + inv(Hess) * grad
  ll(X,y,beta.new)                          # Find the new beta and calculate MLE again and again and again.....
  nt <- c(nt, ll(X,y,beta.new))
  beta.log <- cbind(beta, beta.new) 
  beta =beta.new
  
}

nt

#[1] 394.38937 161.55837 108.22240  86.59502  77.83202  73.99941  73.09228  73.05701  73.05694
#[10]  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694
#[19]  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694
#[28]  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694
#[37]  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694
#[46]  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694
#[55]  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694
#[64]  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694
#[73]  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694  73.05694


par(new=TRUE)                               # Plot two curve in one chart

plot(nt,log='xy', ylim=c(70,400), type='l', col='red')
