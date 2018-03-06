###### data setting ############################################
##### For exercise1
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

envi <- read.csv("enviANDdensity.csv", header=T) # load data
f <- envi$FishDensity..ind..1000m3. # fish density
c <- envi$CopepodDensity.ind..m3. # copepod density
n <- length(f) # sample size

##### For exercise2
cope <- read.csv("copepod_datasheet.csv", skip=1, header=T)
cope_cut <- cope[c(which(cope[, 2] == "Oncaea venusta"), which(cope[, 2] == "Canthocalanus pauper")), ]
rownames(cope_cut) <- cope_cut[, 2]
cope_cut <- cope_cut[, -c(1, 2)]
oncaea <- as.numeric(matrix(cope_cut[1, ]))
canthocalanus <- as.numeric(matrix(cope_cut[2, ]))



########## Exercise1

X <- cbind(rep(1, n), c) # regressor
b <- solve(t(X)%*%X)%*%t(X)%*%f # regression coefficients

pair <- cbind(f, c) # bind fish and copepod data
# write a function for randomization (on copepod), and computing regression coeffitions
rand.coeff <- function(x){
    X1 <- cbind(rep(1, n), c[order(runif(n, 0, 1))])
    b1 <- solve(t(X1)%*%X1)%*%t(X1)%*%f
    b1[2]
}

rand.b <- matrix(0, 1, 5000) # 0 matrix for randomized coeffitions
rand.b <- apply(rand.b, 2, FUN=function(x){x=rand.coeff()}) # generate 5000 randomized beta1

hist(rand.b, main="Histogram of randomized beta1", xlab="Beta1")
sum((b[2] - abs(rand.b)) <= 0)/5000 # p-value under 2-tailed test



########## Exercise2

diff <- mean(oncaea) - mean(canthocalanus) # original difference
# write a function for randomization (on all data), and computing the difference
rand.diff <- function(x){
    data <- c(oncaea, canthocalanus)
    n <- length(data)
    data <- data[order(runif(n, 0, 1))] # randomize the original data first
    newdata <- matrix(data[order(runif(n, 0, 1))], nrow=n/2) # randomly fill up the matrix
    mean(newdata[, 1]) - mean(newdata[, 2])
}

rand.d <- matrix(0, 1, 5000) # 0 matrix for randomized differences
rand.d <- apply(rand.d, 2, FUN=function(x){x=rand.diff()})

hist(rand.d, main="Histogram of randomized difference", xlab="Difference")
sum((diff - abs(rand.d)) <= 0)/5000 # p-value under 2-tailed test


