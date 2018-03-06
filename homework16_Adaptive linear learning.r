##### (a)
n = 1000
# training set
x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

# training test set
x1.test <- runif(n, min=-1, max=1)
x2.test <- runif(n, min=-1, max=1)

# validation set
y <- 3*x1 - 2*x2
y.test <- 3*x1.test - 2*x2.test

 

##### (b)
n.run = 10000  # number of iteration

X <- cbind(1, x1, x2)  # add bias term
W <- matrix(0, nrow=n.run, ncol=3)  # record weights in each iteration
MSE <- numeric(n.run)  # record MSE in each iteration

W[1, ] = runif(3, min=-0.5, max=0.5)  # generate initial values from U(-0.5, 0.5)
MSE[1] = t((y - X%*%W[1, ]))%*%(y - X%*%W[1, ])/n  # initial MSE

a = 0.5  # learning rate

for (i in 2:n.run) {
    W[i, ] = W[(i-1), ] + a*t(X)%*%(y - X%*%W[(i-1), ])/n
    MSE[i] = t((y - X%*%W[i, ]))%*%(y - X%*%W[i, ])/n
}

plot(MSE, type='l', xlab='Iteration')
opti.w <- W[which.min(MSE), ]  # optimal weights that minimize MSE



##### (c)
# We now use these weights to predict the test set.
X1 <- cbind(1, x1.test, x2.test)
y.predict <- X1%*%opti.w
MSE.test <- t((y.test - y.predict))%*%(y.test - y.predict)/n  # MSE for the test set

MSE.test

