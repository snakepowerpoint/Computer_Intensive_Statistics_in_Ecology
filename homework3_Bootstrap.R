setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

########## Exercise1
envi <- read.csv("enviANDdensity.csv", header=T)

##### normal theory
### fish
mean(envi$FishDensity..ind..1000m3.) # mean
sqrt(var(envi$FishDensity..ind..1000m3.)/34) # SE(mean)

### copepod
mean(envi$CopepodDensity.ind..m3.) # mean 
sqrt(var(envi$CopepodDensity.ind..m3.)/34) # SE(mean)

##### non-parametric bootstrap
### fish
# we first create a function for sampling with replacement
sam <- function(x){
    l <- length(x) # record length of x
    m <- matrix(0, nrow=l, ncol=l) # create 0 matrix
    m <- apply(m, 2, FUN=function(x){runif(l, 0, 1)}) 
    # generate random variables from U(0,1)
    m <- apply(m, 2, FUN=function(x){order(x)})
    # find their order, then we can get randomly selected number
    x[m[1, ]] # sample x using random number(each row can be a random number)
}

f <- envi$FishDensity..ind..1000m3. # record fish density
fish.b <- matrix(0, nrow=34, ncol=999) # 0 matrix
fish.b <- apply(fish.b, 2, FUN=function(x){sam(f)}) # bootstrap it with my function
fish.b <- cbind(fish.b, f) # the last column is the original sample

mean(apply(fish.b, 2, mean)) # mean of "mean of each bootstrap data" 
mean(f) # compare with the original sample mean, it seems reliable

sqrt(var(apply(fish.b, 2, mean))) # SE(mean) under bootstrap
sqrt(var(f)/34) # compare with the original SE(mean), it seems reliable

### copepod
# repeat all the things done above
c <- envi$CopepodDensity.ind..m3.
cop.b <- matrix(0, nrow=34, ncol=999)
cop.b <- apply(cop.b, 2, FUN=function(x){sam(c)})
cop.b <- cbind(cop.b, c)

mean(apply(cop.b, 2, mean))
mean(c)

sqrt(var(apply(cop.b, 2, mean)))
sqrt(var(c)/34)

##### histogram of bootstrapped mean
hist(apply(fish.b, 2, mean), main="Bootstrapped mean of fish",
     xlab="mean") # for fish
hist(apply(cop.b, 2, mean), main="Bootstrapped mean of copepod",
     xlab="mean") # for copepod



########## Exercise2
# we use the generated data above
##### fish
median(f) # median of the original fish density
sqrt(var(apply(fish.b, 2, median))) # SE(median) of bootstrapped median of fish

##### copepod
median(c) # median of the original copepod density
sqrt(var(apply(cop.b, 2, median))) # SE(median) of bootstrapped median of copepod

##### histogram of bootstrapped median
hist(apply(fish.b, 2, median), main="Bootstrapped median of fish",
     xlab="median") # for fish
hist(apply(cop.b, 2, median), main="Bootstrapped median of copepod",
     xlab="median") # for copepod



########## Exercise3
##### plot and calculate regression coefficients
plot(c, f, pch=20, xlab="Copepod", ylab="Fish", main="Fish vs. Copepod")
X <- cbind(rep(1, 34), c) # generate variable matrix X
b <- solve(t(X)%*%X)%*%t(X)%*%f # regression coefficients
abline(a=b[1], b=b[2], col="red") # regression line

##### bootstrap beta0 and beta1
pair <- cbind(f, c) # bind fish and copepod data
b.b <- matrix(0, nrow=2, ncol=999) # create a 0 matrix to record coefficients

# write a function to calculate regression coefficients
coeff <- function(x){
    data <- pair[sam(1:34), ]
    Y1 <- data[, 1]
    X1 <- cbind(rep(1, 34), data[, 2])
    solve(t(X1)%*%X1)%*%t(X1)%*%Y1
}

b.b <- apply(b.b, 2, FUN=function(x){x=coeff()}) # record bootstrapped data
b.b <- cbind(b.b, b) # add the original data into bootstrapped data

sd(b.b[1, ]) # SE(beta0)
sd(b.b[2, ]) # SE(beta1)

##### histogram for regression coefficients
hist(b.b[1, ], main="Bootstrapped beta0", xlab="Beta0") # intercept
hist(b.b[2, ], main="Bootstrapped beta1", xlab="Beta1") # slope

