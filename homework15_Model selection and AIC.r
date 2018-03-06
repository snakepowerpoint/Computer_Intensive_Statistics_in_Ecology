###### data setting ###################################################
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

ceta <- read.table("cetaceancorrect.dat")

x <- ifelse(ceta$V2 > 3, 1, 0)  # turn V2 into binary data
ceta1 <- data.frame(cbind(x, ceta[, 3:4]))
colnames(ceta1) <- c("x", "y", "z")

########## Exercise1
##### Full model

# -log(L(par|x)) function, which should be minimized
neg.loglike.1 <- function(par){
    b0 <- par[1]
    b1 <- par[2]
    b2 <- par[3]
    b3 <- par[4]
    sum(log(1 + exp(b0*1 + b1*x + b2*y + b3*y^2)) - z*(b0*1 + b1*x + b2*y + b3*y^2))
}

x <- ceta1$x
y <- ceta1$y
z <- ceta1$z

# set up different inital values
value <- seq(-0.1, 0.1, by=0.005)
ini <- matrix(rep(value, each=4), nrow=4)

# record the -loglikelihood value under each initial value
likelihood <- numeric(41)  # 41 is the length of "value"

for (i in 1:41) {
    likelihood[i] = with(ceta1, optim(par=ini[, i], fn=neg.loglike.1))$value
}

table(round(likelihood, 2))  # 365.74 is the most frequent, No.variable = 4
2*4 + 2*as.numeric(names(table(round(likelihood, 2))[1]))  # AIC = 739.48

# findout which initial value works better than others
which.min(likelihood)  # 21
b1 <- with(ceta1, optim(par=ini[, 21], fn=neg.loglike.1))$par  # regression coefficient

d <- matrix(c(rep(1, 1000), x, y, y^2), nrow=1000)  # data with intercept term
pi1 <- 1/(1 + exp(-d%*%b1))  # each pi

plot(x=y, y=pi1, pch=19, cex=0.8, main="Full model")


##### Model without beta1

neg.loglike.2 <- function(par){
    b0 <- par[1]
    b2 <- par[2]
    b3 <- par[3]
    sum(log(1 + exp(b0*1 + b2*y + b3*y^2)) - z*(b0*1 + b2*y + b3*y^2))
}

ini <- matrix(rep(value, each=3), nrow=3)

for (i in 1:41) {
    likelihood[i] = with(ceta1, optim(par=ini[, i], fn=neg.loglike.2))$value
}

table(round(likelihood, 2))  # 368.98 is the most frequent, No.variable = 3
2*3 + 2*as.numeric(names(table(round(likelihood, 2))[1]))  # AIC = 743.96

which.min(likelihood)  # 12
b2 <- with(ceta1, optim(par=ini[, 12], fn=neg.loglike.2))$par  

d <- matrix(c(rep(1, 1000), y, y^2), nrow=1000) 
pi2 <- 1/(1 + exp(-d%*%b2))

plot(x=y, y=pi2, pch=19, cex=0.8, main="Model without beta1")


##### Model without beta2

neg.loglike.3 <- function(par){
    b0 <- par[1]
    b1 <- par[2]
    b3 <- par[3]
    sum(log(1 + exp(b0*1 + b1*x + b3*y^2)) - z*(b0*1 + b1*x + b3*y^2))
}

for (i in 1:41) {
    likelihood[i] = with(ceta1, optim(par=ini[, i], fn=neg.loglike.3))$value
}

table(round(likelihood, 2))  # 374.13 is the most frequent, No.variable = 3
2*3 + 2*as.numeric(names(table(round(likelihood, 2))[1]))  # AIC = 754.26

which.min(likelihood)  # 19
b3 <- with(ceta1, optim(par=ini[, 19], fn=neg.loglike.3))$par  

d <- matrix(c(rep(1, 1000), x, y^2), nrow=1000) 
pi3 <- 1/(1 + exp(-d%*%b3))

plot(x=y, y=pi3, pch=19, cex=0.8, main="Model without beta2")


##### Model without beta3

neg.loglike.4 <- function(par){
    b0 <- par[1]
    b1 <- par[2]
    b2 <- par[3]
    sum(log(1 + exp(b0*1 + b1*x + b2*y)) - z*(b0*1 + b1*x + b2*y))
}

for (i in 1:41) {
    likelihood[i] = with(ceta1, optim(par=ini[, i], fn=neg.loglike.4))$value
}

table(round(likelihood, 2))  # 369.91 is the most frequent, No.variable = 3
2*3 + 2*as.numeric(names(table(round(likelihood, 2))[1]))  # AIC = 745.82

which.min(likelihood)  # 9
b4 <- with(ceta1, optim(par=ini[, 9], fn=neg.loglike.4))$par  

d <- matrix(c(rep(1, 1000), x, y), nrow=1000) 
pi4 <- 1/(1 + exp(-d%*%b4))

plot(x=y, y=pi4, pch=19, cex=0.8, main="Model without beta3")


##### Full model adding interaction term x*y

neg.loglike.5 <- function(par){
    b0 <- par[1]
    b1 <- par[2]
    b2 <- par[3]
    b3 <- par[4]
    b4 <- par[5]
    sum(log(1 + exp(b0*1 + b1*x + b2*y + b3*y^2 + b4*(x*y))) - z*(b0*1 + b1*x + b2*y + b3*y^2 + b4*(x*y)))
}

ini <- matrix(rep(value, each=5), nrow=5)

for (i in 1:41) {
    likelihood[i] = with(ceta1, optim(par=ini[, i], fn=neg.loglike.5))$value
}

table(round(likelihood, 2))  # 365.62 is the smallest, No.variable = 5
2*5 + 2*as.numeric(names(table(round(likelihood, 2))[1]))  # AIC = 741.24

which.min(likelihood)  # 9
b5 <- with(ceta1, optim(par=ini[, 9], fn=neg.loglike.5))$par  

d <- matrix(c(rep(1, 1000), x, y, y^2, x*y), nrow=1000) 
pi5 <- 1/(1 + exp(-d%*%b5))

plot(x=y, y=pi5, pch=19, cex=0.8, main="Full model adding interaction term x*y")



########## Exercise2
##### Model 1

# separate data
train1 <- ceta1[-c(1:200), ]
train2 <- ceta1[-c(201:400), ]
train3 <- ceta1[-c(401:600), ]
train4 <- ceta1[-c(601:800), ]
train5 <- ceta1[-c(801:1000), ]
train <- list(train1, train2, train3, train4, train5)

test1 <- ceta1[c(1:200), ]
test2 <- ceta1[c(201:400), ]
test3 <- ceta1[c(401:600), ]
test4 <- ceta1[c(601:800), ]
test5 <- ceta1[c(801:1000), ]
test <- list(test1, test2, test3, test4, test5)

# result matrix
result1 <- matrix(0, nrow=5, ncol=2)
beta1 <- matrix(0, nrow=5, ncol=4)

ini <- matrix(rep(value, each=4), nrow=4)

for (i in 1:5) {
    x = train[[i]]$x
    y = train[[i]]$y
    z = train[[i]]$z
    
    for (j in 1:41) {
        likelihood[j] = optim(par=ini[, j], fn=neg.loglike.1)$value
    }
    
    a = which.min(likelihood)
    b = optim(par=ini[, a], fn=neg.loglike.1)
    
    result1[i, 1] = b$value
    beta1[i, ] = b$par
    
    x = test[[i]]$x
    y = test[[i]]$y 
    z = test[[i]]$z   
    d = matrix(c(rep(1, 200), x, y, y^2), nrow=200)
    
    result1[i, 2] = sum(log(1 + exp(d%*%c(b$par))) - z*(d%*%c(b$par)))
}


##### Model2

result2 <- matrix(0, nrow=5, ncol=2)
beta2 <- matrix(0, nrow=5, ncol=3)

ini <- matrix(rep(value, each=3), nrow=3)

for (i in 1:5) {
    x = train[[i]]$x
    y = train[[i]]$y
    z = train[[i]]$z
    
    for (j in 1:41) {
        likelihood[j] = optim(par=ini[, j], fn=neg.loglike.2)$value
    }
    
    a = which.min(likelihood)
    b = optim(par=ini[, a], fn=neg.loglike.2)
    
    result2[i, 1] = b$value
    beta2[i, ] = b$par
    
    x = test[[i]]$x
    y = test[[i]]$y 
    z = test[[i]]$z   
    d = matrix(c(rep(1, 200), y, y^2), nrow=200)
    
    result2[i, 2] = sum(log(1 + exp(d%*%c(b$par))) - z*(d%*%c(b$par)))
}


##### Model3

result3 <- matrix(0, nrow=5, ncol=2)
beta3 <- matrix(0, nrow=5, ncol=3)

for (i in 1:5) {
    x = train[[i]]$x
    y = train[[i]]$y
    z = train[[i]]$z
    
    for (j in 1:41) {
        likelihood[j] = optim(par=ini[, j], fn=neg.loglike.3)$value
    }
    
    a = which.min(likelihood)
    b = optim(par=ini[, a], fn=neg.loglike.3)
    
    result3[i, 1] = b$value
    beta3[i, ] = b$par
    
    x = test[[i]]$x
    y = test[[i]]$y 
    z = test[[i]]$z   
    d = matrix(c(rep(1, 200), x, y^2), nrow=200)
    
    result3[i, 2] = sum(log(1 + exp(d%*%c(b$par))) - z*(d%*%c(b$par)))
}


##### Model4

result4 <- matrix(0, nrow=5, ncol=2)
beta4 <- matrix(0, nrow=5, ncol=3)

for (i in 1:5) {
    x = train[[i]]$x
    y = train[[i]]$y
    z = train[[i]]$z
    
    for (j in 1:41) {
        likelihood[j] = optim(par=ini[, j], fn=neg.loglike.4)$value
    }
    
    a = which.min(likelihood)
    b = optim(par=ini[, a], fn=neg.loglike.4)
    
    result4[i, 1] = b$value
    beta4[i, ] = b$par
    
    x = test[[i]]$x
    y = test[[i]]$y 
    z = test[[i]]$z   
    d = matrix(c(rep(1, 200), x, y), nrow=200)
    
    result4[i, 2] = sum(log(1 + exp(d%*%c(b$par))) - z*(d%*%c(b$par)))
}


##### Model5

result5 <- matrix(0, nrow=5, ncol=2)
beta5 <- matrix(0, nrow=5, ncol=5)

ini <- matrix(rep(value, each=5), nrow=5)

for (i in 1:5) {
    x = train[[i]]$x
    y = train[[i]]$y
    z = train[[i]]$z
    
    for (j in 1:41) {
        likelihood[j] = optim(par=ini[, j], fn=neg.loglike.5)$value
    }
    
    a = which.min(likelihood)
    b = optim(par=ini[, a], fn=neg.loglike.5)
    
    result5[i, 1] = b$value
    beta5[i, ] = b$par
    
    x = test[[i]]$x
    y = test[[i]]$y 
    z = test[[i]]$z   
    d = matrix(c(rep(1, 200), x, y, y^2, x*y), nrow=200)
    
    result5[i, 2] = sum(log(1 + exp(d%*%c(b$par))) - z*(d%*%c(b$par)))
}


### Summay
apply(result1, 2, mean)
apply(result2, 2, mean)
apply(result3, 2, mean)
apply(result4, 2, mean)
apply(result5, 2, mean)

beta1
beta2
beta3
beta4
beta5


