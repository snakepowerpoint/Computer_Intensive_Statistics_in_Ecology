###### data setting ###################################################
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

rawdata <- read.table("modeldata.txt")

x.test <- as.matrix(cbind(1, rawdata[, 1:2]))
y.valid <- as.matrix(rawdata[, 3:5])

n <- dim(rawdata)[1]  # number of data point
n.run <- 1000  # number of learning step
rate <- 0.001  # learning rate
# tanh(x) is the activation function
dtanh <- function(x) {1 - tanh(x)^2}  # derivative of tanh(x)

# There are 5 hidden neurons, and we suppose that each hidden neuron
# contains 2 units. Remember to add bias term.

#############################################################
# Start from here.
# hidden unit value
S <- matrix(0, nrow=n, ncol=2)  # hidden neuron 1
H <- matrix(0, nrow=n, ncol=2)  # hidden neuron 2
P <- matrix(0, nrow=n, ncol=2)  # hidden neuron 3
Q <- matrix(0, nrow=n, ncol=2)  # hidden neuron 4
R <- matrix(0, nrow=n, ncol=2)  # hidden neuron 5

# coefficient 1
a1 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 1, unit 1
b1 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 2, unit 1
c1 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 3, unit 1
d1 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 4, unit 1
e1 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 5, unit 1

# coefficient 2
a2 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 1, unit 2
b2 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 2, unit 2
c2 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 3, unit 2
d2 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 4, unit 2
e2 <- matrix(0, nrow=n.run, ncol=2+1)  # hidden neuron 5, unit 2

# final output weights
f1 <- matrix(0, nrow=n.run, ncol=2+1)
f2 <- matrix(0, nrow=n.run, ncol=2+1)
f3 <- matrix(0, nrow=n.run, ncol=2+1)

# final output
y.result <- matrix(0, nrow=n.run, ncol=3)

# mse
mse = matrix(0, nrow=n.run, ncol=3)

# set initial values
a1[1, ] = runif(3, -0.5, 0.5)
b1[1, ] = runif(3, -0.5, 0.5)
c1[1, ] = runif(3, -0.5, 0.5)
d1[1, ] = runif(3, -0.5, 0.5)
e1[1, ] = runif(3, -0.5, 0.5)

a2[1, ] = runif(3, -0.5, 0.5)
b2[1, ] = runif(3, -0.5, 0.5)
c2[1, ] = runif(3, -0.5, 0.5)
d2[1, ] = runif(3, -0.5, 0.5)
e2[1, ] = runif(3, -0.5, 0.5)

f1[1, ] = runif(3, -0.5, 0.5)
f2[1, ] = runif(3, -0.5, 0.5)
f3[1, ] = runif(3, -0.5, 0.5)

# Now we compute MSE for step 1
S = tanh(x.test%*%cbind(a1[1, ], a2[1, ]))
H = tanh(cbind(1, S)%*%cbind(b1[1, ], b2[1, ]))
P = tanh(cbind(1, H)%*%cbind(c1[1, ], c2[1, ]))
Q = tanh(cbind(1, P)%*%cbind(d1[1, ], d2[1, ]))
R = tanh(cbind(1, Q)%*%cbind(e1[1, ], e2[1, ]))
y.result = tanh(cbind(1, R)%*%cbind(f1[1, ], f2[1, ], f3[1, ]))

mse[1, ] = diag(t(y.valid - y.result)%*%(y.valid - y.result))/n


# In backward pass process, we first compute the next time values mannually.
delta1 <- matrix(0, nrow=n, ncol=3)
delta2 <- matrix(0, nrow=n, ncol=2)
delta3 <- matrix(0, nrow=n, ncol=2)
delta4 <- matrix(0, nrow=n, ncol=2)
delta5 <- matrix(0, nrow=n, ncol=2)
delta6 <- matrix(0, nrow=n, ncol=2)

delta1 = dtanh(cbind(1, R)%*%cbind(f1[1, ], f2[1, ], f3[1, ]))*(y.valid - y.result)
delta2 = dtanh(cbind(1, Q)%*%cbind(e1[1, ], e2[1, ]))*
    (delta1%*%t(cbind(f1[1, ], f2[1, ], f3[1, ])[-1, ]))
delta3 = dtanh(cbind(1, P)%*%cbind(d1[1, ], d2[1, ]))*
    (delta2%*%t(cbind(e1[1, ], e2[1, ])[-1, ]))
delta4 = dtanh(cbind(1, H)%*%cbind(c1[1, ], c2[1, ]))*
    (delta3%*%t(cbind(d1[1, ], d2[1, ])[-1, ]))
delta5 = dtanh(cbind(1, S)%*%cbind(b1[1, ], b2[1, ]))*
    (delta4%*%t(cbind(c1[1, ], c2[1, ])[-1, ]))
delta6 = dtanh(x.test%*%cbind(a1[1, ], a2[1, ]))*
    (delta5%*%t(cbind(b1[1, ], b2[1, ])[-1, ]))

f1[2, ] = rate*delta1[, 1]%*%cbind(1, R)/n + f1[1, ]
f2[2, ] = rate*delta1[, 2]%*%cbind(1, R)/n + f2[1, ]
f3[2, ] = rate*delta1[, 3]%*%cbind(1, R)/n + f3[1, ]

e1[2, ] = rate*delta2[, 1]%*%cbind(1, Q)/n + e1[1, ]
d1[2, ] = rate*delta3[, 1]%*%cbind(1, P)/n + d1[1, ]
c1[2, ] = rate*delta4[, 1]%*%cbind(1, H)/n + c1[1, ]
b1[2, ] = rate*delta5[, 1]%*%cbind(1, S)/n + b1[1, ]
a1[2, ] = rate*delta6[, 1]%*%x.test/n + a1[1, ]

e2[2, ] = rate*delta2[, 2]%*%cbind(1, Q)/n + e2[1, ]
d2[2, ] = rate*delta3[, 2]%*%cbind(1, P)/n + d2[1, ]
c2[2, ] = rate*delta4[, 2]%*%cbind(1, H)/n + c2[1, ]
b2[2, ] = rate*delta5[, 2]%*%cbind(1, S)/n + b2[1, ]
a2[2, ] = rate*delta6[, 2]%*%x.test/n + a2[1, ]

# Thus, we have new coefficients now.
# Run the loop for the following steps.
for (i in 2:(n.run-1)) {
    S = tanh(x.test%*%cbind(a1[i, ], a2[i, ]))
    H = tanh(cbind(1, S)%*%cbind(b1[i, ], b2[i, ]))
    P = tanh(cbind(1, H)%*%cbind(c1[i, ], c2[i, ]))
    Q = tanh(cbind(1, P)%*%cbind(d1[i, ], d2[i, ]))
    R = tanh(cbind(1, Q)%*%cbind(e1[i, ], e2[i, ]))
    y.result = tanh(cbind(1, R)%*%cbind(f1[i, ], f2[i, ], f3[i, ]))
    
    mse[i, ] = diag(t(y.valid - y.result)%*%(y.valid - y.result))/n
    
    delta1 = dtanh(cbind(1, R)%*%cbind(f1[i, ], f2[i, ], f3[i, ]))*(y.valid - y.result)
    delta2 = dtanh(cbind(1, Q)%*%cbind(e1[i, ], e2[i, ]))*
        (delta1%*%t(cbind(f1[i, ], f2[i, ], f3[i, ])[-1, ]))
    delta3 = dtanh(cbind(1, P)%*%cbind(d1[i, ], d2[i, ]))*
        (delta2%*%t(cbind(e1[i, ], e2[i, ])[-1, ]))
    delta4 = dtanh(cbind(1, H)%*%cbind(c1[i, ], c2[i, ]))*
        (delta3%*%t(cbind(d1[i, ], d2[i, ])[-1, ]))
    delta5 = dtanh(cbind(1, S)%*%cbind(b1[i, ], b2[i, ]))*
        (delta4%*%t(cbind(c1[i, ], c2[i, ])[-1, ]))
    delta6 = dtanh(x.test%*%cbind(a1[i, ], a2[i, ]))*
        (delta5%*%t(cbind(b1[i, ], b2[i, ])[-1, ]))
    
    f1[(i+1), ] = rate*delta1[, 1]%*%cbind(1, R)/n + f1[i, ]
    f2[(i+1), ] = rate*delta1[, 2]%*%cbind(1, R)/n + f2[i, ]
    f3[(i+1), ] = rate*delta1[, 3]%*%cbind(1, R)/n + f3[i, ]
    
    e1[(i+1), ] = rate*delta2[, 1]%*%cbind(1, Q)/n + e1[i, ]
    d1[(i+1), ] = rate*delta3[, 1]%*%cbind(1, P)/n + d1[i, ]
    c1[(i+1), ] = rate*delta4[, 1]%*%cbind(1, H)/n + c1[i, ]
    b1[(i+1), ] = rate*delta5[, 1]%*%cbind(1, S)/n + b1[i, ]
    a1[(i+1), ] = rate*delta6[, 1]%*%x.test/n + a1[i, ]
    
    e2[(i+1), ] = rate*delta2[, 2]%*%cbind(1, Q)/n + e2[i, ]
    d2[(i+1), ] = rate*delta3[, 2]%*%cbind(1, P)/n + d2[i, ]
    c2[(i+1), ] = rate*delta4[, 2]%*%cbind(1, H)/n + c2[i, ]
    b2[(i+1), ] = rate*delta5[, 2]%*%cbind(1, S)/n + b2[i, ]
    a2[(i+1), ] = rate*delta6[, 2]%*%x.test/n + a2[i, ]
}

# Remarkably, it is very very fast XDDDD.
y.result  # seems not resonable
plot(mse[, 1])  # MSE declines slowly

# try different rate, and run all of the things above
rate = 0.5

#-----------------------------------------------
# After run all of the things above....
y.result  # seems resonable now
plot(mse[, 1])
plot(mse[, 2])
plot(mse[, 3])
plot(apply(mse, 1, sum)[-n.run], type='l', xlab='Step', ylab='MSE')

# find the minimum MSE
min.mse <- which.min(apply(mse, 1, sum)[-n.run])

# compute the final codings with minimal MSE
S = tanh(x.test%*%cbind(a1[min.mse, ], a2[min.mse, ]))
H = tanh(cbind(1, S)%*%cbind(b1[min.mse, ], b2[min.mse, ]))
P = tanh(cbind(1, H)%*%cbind(c1[min.mse, ], c2[min.mse, ]))
Q = tanh(cbind(1, P)%*%cbind(d1[min.mse, ], d2[min.mse, ]))
R = tanh(cbind(1, Q)%*%cbind(e1[min.mse, ], e2[min.mse, ]))
y.result = tanh(cbind(1, R)%*%cbind(f1[min.mse, ], f2[min.mse, ], f3[min.mse, ]))

# coefficients with minimal MSE
a1[min.mse, ]
b1[min.mse, ]
c1[min.mse, ]
d1[min.mse, ]
e1[min.mse, ]

a2[min.mse, ]
b2[min.mse, ]
c2[min.mse, ]
d2[min.mse, ]
e2[min.mse, ]

f1[min.mse, ]
f2[min.mse, ]
f3[min.mse, ]

# plot the original data, and compare it with the output data
y1 = ifelse(rawdata[, 3]==1, 1, 0)
y2 = ifelse(rawdata[, 4]==1, 2, 0)
y3 = ifelse(rawdata[, 5]==1, 3, 0)
y <- y1 + y2 + y3
plot(x.test[, -1], col=y, pch=19, xlab='', ylab='', main="Original data")

y1 = ifelse(floor(y.result)[, 1]==0, 1, 0)
y2 = ifelse(floor(y.result)[, 2]==0, 2, 0)
y3 = ifelse(floor(y.result)[, 3]==0, 3, 0)
y <- y1 + y2 + y3
plot(x.test[, -1], col=y, pch=19, xlab='', ylab='', main="Output data")
# We see that some points are classified into wrong classes, especially the upper-right corner.
# Try longer steps.
n.run = 5000
rate = 0.5

