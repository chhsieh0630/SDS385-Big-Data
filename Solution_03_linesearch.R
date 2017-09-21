################################ Exercise 3 ##################################
################################ LineSearch ##################################

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


hess <- function(X,beta) {                  # Find the Hessian Matrix to calculate by Newton method
                                            # Hessian Matrix = \/2 * f(x_k)
  H = crossprod(diag(as.vector(sqrt(w(X,beta) * (1 - w(X,beta))))) %*% X)
  return(H)
}


rm_step <- function(C,t,alpha) {            # step r = C(t+t0)^(-a)
  r <- C*(t+2)^(-alpha)
  return(r)
}  



linesearch <- function(X,y,beta,rho,c){   # f(beta-alpha*p_k) < f(beta)+c*alpha*grad^2
  alpha = 1
  r1 <- ll(X,y,beta)
  r2 <- c * alpha * (  -t(grad(X,y,beta)) %*% grad(X,y,beta)  ) # direction = -grad

  beta.new = beta - alpha * grad(X,y,beta)                      # b.new = b + a*direction
  left <- ll(X,y,beta.new)
  
  while(left > r1 + r2) {
    alpha = rho * alpha
    r2 <- c * alpha * (  -t(grad(X,y,beta)) %*% grad(X,y,beta)  )
    beta.new = beta - alpha * grad(X,y,beta)
    left <- ll(X,y,beta.new)
  }
  return(alpha)
}


beta <- matrix(0,nrow=11,ncol=1)
backline <- ll(X,y,beta)
rho = 0.5
c = 0.5
alpha <- linesearch(X,y,beta,rho,c)
alpha_record <- alpha                     # trace alpha for debugging

for (i in 1:5000) {
  alpha <- linesearch(X,y,beta,rho,c)
  beta = beta - alpha * grad(X,y,beta)    # Beta_k+1 = Beta_k + alpha_k * grad
  alpha_record <- c(alpha_record,alpha)   # trace alpha for debugging
  backline <- c(backline, ll(X,y,beta)) 
}


plot(backline,log='xy', xlim=c(1,5000),ylim=c(70,400), type='l', col='blue')

