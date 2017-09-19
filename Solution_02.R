##################################### Exercise 2 #########################################
################### Stochastic gradient descent - Robbins Monro ##########################

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

alpha = 0.001            

for (i in 1:10000) {

  beta.new = beta - alpha * grad(X,y,beta)    # Beta_k+1 = Beta_k + alpha_k * grad
  ll(X,y,beta.new)                            # Find the new beta and calculate MLE again and again and again.....
  mle <- c(mle, ll(X,y,beta.new)) 
  beta.log <- cbind(beta, beta.new) 
  beta = beta.new

}
  

mle
plot(mle,log='xy', xlim=c(1,10000),ylim=c(70,400), type='l', col='blue', xlab="iterations", 
     ylab="negative loglike")


################################## Newton Method ####################################

hess <- function(X,beta) {                  # Find the Hessian Matrix to calculate by Newton method
                                            # Hessian Matrix = \/2 * f(x_k)
  
  H = crossprod(diag(as.vector(sqrt(w(X,beta) * (1 - w(X,beta))))) %*% X)
  return(H)
}

beta <- matrix(0,nrow=11,ncol=1)
nt <- ll(X,y,beta)
beta.log <-matrix(0,nrow=11,ncol=1)


for (i in 1:10000) {
  
  beta.new = beta - solve(hess(X,beta)) %*% grad(X,y,beta)   # Beta_k+1 = Beta_k + inv(Hess) * grad
  ll(X,y,beta.new)                          # Find the new beta and calculate MLE again and again and again.....
  nt <- c(nt, ll(X,y,beta.new))
  beta.log <- cbind(beta, beta.new) 
  beta =beta.new
  
}

nt
lines(nt, type="l", col="red")

######################################## Robbins-Monro ################################

rm_step <- function(C,t,alpha) {            # step r = C(t+t0)^(-a)
  r <- C*(t+2)^(-alpha)
  return(r)
}  

C=0.01
alpha=0.5
beta_SGD <- matrix(0,nrow=11,ncol=1)
rm <- ll(X,y,beta_SGD)
rm_record <- rm_step(C,0,alpha)

for (i in 1:10000) {
  
  beta_SGD = beta_SGD - rm_step(C,i,alpha) * grad(X,y,beta_SGD)  
  rm_record <- c(rm_record,rm_step(C,i,alpha))
  rm <- c(rm, ll(X,y,beta_SGD))
  
}

lines(rm, type="l", col="forestgreen")
legend(500, 350, c("GLM","Newton","Robbins-Monro"), col=c("blue","red","forestgreen"),lty=1)
