###### data setting ###################################################
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

########## Exercise
#####  Metropolis-Hastings method

# generate data from N(m=0.06, 1)
d <- rnorm(20, mean=0.06, sd=1)

# set up posterior distribution (target distribution)
# the constant terms are ignored
post <- function(m, data){
    n = length(data)
    mu = mean(data)
    exp(-n*((m - mu)^2)/2)/(1 + m^2)
}

# set up proposal distribution (assume Normal)
# ignore the constant terms
prop <- function(x, mu, sig){1/sig*exp(-0.5*((x-mu)/sig)^2)}

# Generate 10000 samples in the chain.
# Set up the constants.
n = 10000
sig = 2
x = numeric(n)
x[1] = rnorm(1)  # generate the starting point

for (k in 2:n) {
    # generate a candidate from the proposal distribution
    # which is the normal in this case. This will be a normal
    # with mean given by the previous value in the chain and 
    # standard deviation of 'sig'
    y = x[k-1] + sig*rnorm(1)
    # generate a uniform for comparison
    u = runif(1)
    alpha = min(c(1, post(m=y, data=d)*prop(x[k-1], y, sig)/
                      (post(m=x[k-1], data=d)*prop(y, x[k-1], sig))))
    if (u <= alpha) {
        x[k] = y
    } 
    else {
        x[k] = x[k-1]
    }
}

plot(x, type='l', xlab="Step")  # the sequence converges

# try sig = 0.5
sig = 0.5
x1 = numeric(n)
x1[1] = rnorm(1)  # new starting value

for (k in 2:n) {
    y = x1[k-1] + sig*rnorm(1)
    u = runif(1)
    alpha = min(c(1, post(m=y, data=d)*prop(x1[k-1], y, sig)/
                      (post(m=x1[k-1], data=d)*prop(y, x1[k-1], sig))))
    if (u <= alpha) {
        x1[k] = y
    } 
    else {
        x1[k] = x1[k-1]
    }
}

plot(x1, type='l', xlab="Step")  # the result is similar to the previous one

# To compute the mean and the variance of posterior distribution, we need to determine 
# when to stop the sequence first. That is, from which step will the sequence converge?



##### Gelman-Rubin method

# We will generate 10000 iterations for the chain. (4 chains)
numchain = 4
# Set up the vectors to store the samples.
# This is 4 chains, 10000 samples.
X = matrix(0, numchain, n)
# This is 4 sequences (rows) of summaries.
nu = matrix(0, numchain, n)
# Track the rhat for each iteration
rhat = numeric(n)
# Get the starting values for the chain
# Use over-dispersed starting points
X[1,1] = -8
X[2,1] = 8
X[3,1] = -5
X[4,1] = 5

# The following implements the chains. Note that each column of our matrices X and nu
# is one iteration of the chains, and each row contains one of the chains.
# The X matrix keeps the chains, and the matrix nu is the sequence of scalar summaries
# of each chain

source("csgelrub.R")

### for sample mean

# Run the chain
for (j in 2:n) {
    for (i in 1:numchain) {
        # Generate variate from proposal distribution.
        y = rnorm(1)*sig + X[i, j-1]
        # Generate variate from uniform
        u = runif(1)
        # Calculate alpha
        alpha = min(c(1, post(m=y, data=d)*prop(X[i, j-1], y, sig)/
                          (post(m=X[i, j-1], data=d)*prop(y, X[i, j-1], sig))))
        if (u <= alpha) {
            # Then set the chain to the y
            X[i, j] = y
        } 
        else {
            X[i, j] = X[i, j-1]
        }
    }
    # Get the scalar summary: mean of each row
    nu[, j] = apply(X[, 1:j], 1, mean)
    rhat[j] = csgelrub(nu[, 1:j])
}
# The function csgelrub will return the estimated for a given set of 
# sequence of scalar summaries

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
for (i in 1:4) {
    plot(nu[i, ], type="l", xlim=c(1,n))
}
# each sequence converges to the same value

par(mfrow=c(1,1)) 
par(mar=c(5,5,5,5))
plot(rhat, type="l", xlab="Iteration of Chain", ylab="R-hat")
# Thus, we compute the mean of posterior distribution using
# the steps after 2000th.

### for sample variance

Y = matrix(0, numchain, n)
nu.y = matrix(0, numchain, n)
rhat.y = numeric(n)

# starting values
Y[1,1] = -8
Y[2,1] = 8
Y[3,1] = -5
Y[4,1] = 5

for (j in 2:n) {
    for (i in 1:numchain) {
        # Generate variate from proposal distribution.
        y = rnorm(1)*sig + Y[i, j-1]
        # Generate variate from uniform
        u = runif(1)
        # Calculate alpha
        alpha = min(c(1, post(m=y, data=d)*prop(Y[i, j-1], y, sig)/
                          (post(m=Y[i, j-1], data=d)*prop(y, Y[i, j-1], sig))))
        if (u <= alpha) {
            # Then set the chain to the y
            Y[i, j] = y
        } 
        else {
            Y[i, j] = Y[i, j-1]
        }
    }
    # Get the scalar summary: variance of each row
    nu.y[, j] = apply(Y[, 1:j], 1, var)
    rhat.y[j] = csgelrub(nu.y[, 1:j])
}

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
for (i in 1:4) {
    plot(nu.y[i, ], type="l", xlim=c(1,n))
}
# each sequence converges to the same value

par(mfrow=c(1,1)) 
par(mar=c(5,5,5,5))
plot(rhat.y, type="l", xlab="Iteration of Chain", ylab="R-hat")
# Thus, we compute the variance of posterior distribution using
# steps after 4000th.

mean(x1[2000:10000])
var(x1[4000:10000])

mean(d)  # original data sample mean
