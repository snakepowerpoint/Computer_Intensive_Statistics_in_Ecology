setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

########## Exercise1
envi <- read.csv("enviANDdensity.csv", header=T) # load data
f <- envi$FishDensity..ind..1000m3. # fish density
c <- envi$CopepodDensity.ind..m3. # copepod density

##### fish
n <- length(f) # sample size
jack.f <- matrix(0, nrow=n-1, ncol=n) # 0 matrix for jackknife sample

# generate jackknife sample
for (i in 1:n){
    jack.f[, i] = f[-i]
}

mean(apply(jack.f, 2, mean)) # mean for the fish using jackknife
mean(f) # compare with the original sample

sqrt(var(apply(jack.f, 2, mean))*(n-1)*(n-1)/n) # SE(mean) under jackknife
sqrt(var(f)/n) # compare with SE(mean) under normal theory

##### copepod
jack.c <- matrix(0, nrow=n-1, ncol=n) # 0 matrix for jackknife sample

# generate jackknife sample
for (i in 1:n){
    jack.c[, i] = c[-i]
}

mean(apply(jack.c, 2, mean)) # mean for the copepod using jackknife
mean(c) # compare with the original sample

sqrt(var(apply(jack.c, 2, mean))*(n-1)*(n-1)/n) # SE(mean) under jackknife
sqrt(var(c)/n) # compare with SE(mean) under normal theory

##### histogram
hist(apply(jack.f, 2, mean), main="Jackknife means of fish", xlab="mean")
hist(apply(jack.c, 2, mean), main="Jackknife means of copepod", xlab="mean")



########## Exercise2
X <- cbind(rep(1, n), c) # regressor
b <- solve(t(X)%*%X)%*%t(X)%*%f # regression coefficients

jack.b <- matrix(0, nrow=2, ncol=n) # 0 matrix for jackknife coefficients

# calculate regression coefficients in each jackknife subsample
for (i in 1:n){
    jack.b[, i] = solve(t(X[-i, ])%*%X[-i, ])%*%t(X[-i, ])%*%f[-i]
}

sqrt(var(jack.b[1, ])*(n-1)*(n-1)/n) # SE(beta0) under jackknife
sqrt(var(jack.b[2, ])*(n-1)*(n-1)/n) # SE(beta1) under jackknife

##### histogram
hist(jack.b[1, ], main="Jackknife beta0", xlab="Beta0")
hist(jack.b[2, ], main="Jackknife beta1", xlab="Beta1")



########## Exercise3
##### exercise1
### fish
boot.f <- matrix(0, nrow=n, ncol=999) # 0 matrix for bootstrap
# generate bootstrapped sample
boot.f <- apply(boot.f, 2, FUN=function(x){x=f[ceiling(runif(34, 0, 34))]})
boot.f <- cbind(boot.f, f) # combine with the orignial sample

mean(f) # mean under normal theroy
mean(apply(jack.f, 2, mean)) # mean under jackknife
mean(apply(boot.f, 2, mean)) # mean under bootstrap

sqrt(var(f)/n) # SE(mean) under normal theroy
sqrt(var(apply(jack.f, 2, mean))*(n-1)*(n-1)/n) # SE(mean) under jackknife
sqrt(var(apply(boot.f, 2, mean))) # SE(mean) under bootstrap

### copepod
boot.c <- matrix(0, nrow=n, ncol=999) # 0 matrix for bootstrap
# generate bootstrapped sample
boot.c <- apply(boot.c, 2, FUN=function(x){x = c[ceiling(runif(34, 0, 34))]})
boot.c <- cbind(boot.c, c) # combine with the orignial sample

mean(c) # mean under normal theroy
mean(apply(jack.c, 2, mean)) # mean under jackknife
mean(apply(boot.c, 2, mean)) # mean under bootstrap

sqrt(var(c)/n) # SE(mean) under normal theroy
sqrt(var(apply(jack.c, 2, mean))*(n-1)*(n-1)/n) # SE(mean) under jackknife
sqrt(var(apply(boot.c, 2, mean))) # SE(mean) under bootstrap

##### exercise2
pair <- cbind(f, c) # bind fish and copepod data

coeff <- function(x){
    data <- pair[ceiling(runif(34, 0, 34)), ]
    X1 <- cbind(rep(1, 34), data[, 2])
    Y1 <- data[, 1]
    solve(t(X1)%*%X1)%*%t(X1)%*%Y1
}

boot.b <- matrix(0, nrow=2, ncol=999) # 0 matrix for bootstrap
boot.b <- apply(boot.b, 2, FUN=function(x){x = coeff()})
boot.b <- cbind(boot.b, b)

b # regression coefficients under normal theory
apply(jack.b, 1, mean) # regression coefficients under jackknife
apply(boot.b, 1, mean) # regression coefficients under bootstrap

u <- f - X%*%b # error term under linear regression model

sqrt(c(t(u)%*%u/(n-2))*diag(solve(t(X)%*%X)))[1] # SE(beta0) under normal theory
sqrt(var(jack.b[1, ])*(n-1)*(n-1)/n) # SE(beta0) under jackknife
sd(boot.b[1, ]) # SE(beta0) under bootstrap

sqrt(c(t(u)%*%u/(n-2))*diag(solve(t(X)%*%X)))[2] # SE(beta1) under normal theory
sqrt(var(jack.b[2, ])*(n-1)*(n-1)/n) # SE(beta1) under jackknife
sd(boot.b[2, ]) # SE(beta1) under bootstrap





