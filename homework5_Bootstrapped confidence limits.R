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

# write a function for calculating regression coeffitions
coeff <- function(x){
    data <- pair[ceiling(runif(34, 0, 34)), ]
    X1 <- cbind(rep(1, 34), data[, 2])
    Y1 <- data[, 1]
    solve(t(X1)%*%X1)%*%t(X1)%*%Y1
}

boot.b <- matrix(0, nrow=2, ncol=999) # 0 matrix for bootstrap
boot.b <- apply(boot.b, 2, FUN=function(x){x = coeff()}) # fill the matrix with regression coeffitions
boot.b <- cbind(boot.b, b) # combind the original coeffitions

boot.b1.sort <- sort(boot.b[2, ]) # sort the bootstrapped beta1
# draw a CDF picture
plot(x=boot.b1.sort, y=1:1000, pch=19, cex=0.5,
     main="Cummulative distribution of beta1", xlab="Beta1", ylab="i")



##### percentile CI method

c(boot.b1.sort[1000*0.025], boot.b1.sort[1000*0.975]) # 95% confidence interval under percentile



##### Bias Correction (BC) method

# find the order of the original coeffition among all bootstrapped coeffitions
position <- which(b[2]==boot.b1.sort)
z0 <- qnorm(position*(1/1000)) # z score

# write a function to find the closest probability
close.p <- function(x, n=1000){
    # x is value, n is length of observation. For example, if n=1000, the empirical cummulative
    # probability should be 1/1000
    p <- c(1:n)/n # probability sequence
    diff <- abs(p-x) # absolute difference between p and x
    closest1 <- which(sort(diff)[1]==diff) # the closest value
    closest2 <- which(sort(diff)[2]==diff) # the second closest value
    c(closest1, closest2) # print the result
}

# Once we know the relative order of pnorm(2*z0 +- 1.96) among empirical CDF,
# we can find out the corresponding bootstrapped coefficient.
# Here, we take average quantile since empirical CDF is not continuous.
lower <- (boot.b1.sort[close.p(x=pnorm(2*z0 - 1.96))[1]] + boot.b1.sort[close.p(x=pnorm(2*z0 - 1.96))[2]])/2
upper <- (boot.b1.sort[close.p(x=pnorm(2*z0 + 1.96))[1]] + boot.b1.sort[close.p(x=pnorm(2*z0 + 1.96))[2]])/2
c(lower, upper) # 95% confidence interval under BC



##### Bias Correction Accelerate (BCa) method

# we first generate jackknife regression coeffitions, and then calculate a
jack.b <- matrix(0, nrow=2, ncol=n) # 0 matrix for jackknife coefficients

# calculate regression coefficients in each jackknife subsample
for (i in 1:n){
    jack.b[, i] = solve(t(X[-i, ])%*%X[-i, ])%*%t(X[-i, ])%*%f[-i]
}

num <- sum((mean(jack.b[2, ]) - jack.b[2, ])^3) # numerator of a
den <- 6*((sum((mean(jack.b[2, ]) - jack.b[2, ])^2))^(3/2)) # denominator of a
a = num/den

# calculate the adjusted quantile
quan <- c(z0 + (z0 - 1.96)/(1 - a*(z0 - 1.96)), z0 + (z0 + 1.96)/(1 - a*(z0 + 1.96)))

# Similar with previous method, we first find out the relative order of pnorm(quan) among empirical CDF,
# then we can know which bootstrapped coefficient should be choose to define the confidence interval.
# Here, we also take average since the empirical CDF is not continuous.
lower.a <- (boot.b1.sort[close.p(x=pnorm(quan[1]))[1]] + boot.b1.sort[close.p(x=pnorm(quan[1]))[2]])/2
upper.a <- (boot.b1.sort[close.p(x=pnorm(quan[2]))[1]] + boot.b1.sort[close.p(x=pnorm(quan[2]))[2]])/2
c(lower.a, upper.a) # 95% confidence interval under BCa

# We summarize the results under three different methods
c(boot.b1.sort[1000*0.025], boot.b1.sort[1000*0.975]) # 95% confidence interval under percentile
c(lower, upper) # 95% confidence interval under BC
c(lower.a, upper.a) # 95% confidence interval under BCa





########## Exercise2

##### generate bootstrap data for Oncaea Venusta
boot.o <- matrix(0, nrow=n, ncol=999) # 0 matrix for bootstrap
# generate bootstrapped sample
boot.o <- apply(boot.o, 2, FUN=function(x){x=oncaea[ceiling(runif(34, 0, 34))]})
boot.o <- cbind(boot.o, oncaea) # combine with the orignial sample

##### generate bootstrap data for Canthocalanus Pauper
boot.c <- matrix(0, nrow=n, ncol=999) # 0 matrix for bootstrap
# generate bootstrapped sample
boot.c <- apply(boot.c, 2, FUN=function(x){x=canthocalanus[ceiling(runif(34, 0, 34))]})
boot.c <- cbind(boot.c, canthocalanus) # combine with the orignial sample

# calculate differences of the mean with bootstrapped data
diff.oc <- apply(boot.o, 2, mean) - apply(boot.c, 2, mean)
# plot histogram for the differences
hist(diff.oc, main="Distribution of the difference of the means", xlab="Difference")
# it seems that the difference is not equal to 0

##### CI method

c(sort(diff.oc)[1000*0.025], sort(diff.oc)[1000*0.975])
# 0 is not in the range above, so we can reject the null hypothesis that the difference=0

# plot the distribution of the differences
plot(density(diff.oc), main="Density of the differences", xlab="Difference")
abline(v=sort(diff.oc)[1000*0.025], col="blue")
abline(v=sort(diff.oc)[1000*0.975], col="blue")
abline(v=0, col="red")



##### BCa method

# we first generate jackknife differences, and then calculate a
jack.diff <- matrix(0, nrow=1, ncol=34) # 0 matrix for jackknife coefficients

# calculate difference in each jackknife subsample
for (i in 1:34){
    jack.diff[, i] = mean(oncaea[-i]) - mean(canthocalanus[-i])
}

num.diff <- sum((mean(jack.diff[1, ]) - jack.diff[1, ])^3) # numerator of a
den.diff <- 6*((sum((mean(jack.diff[1, ]) - jack.diff[1, ])^2))^(3/2)) # denominator of a
a.diff = num.diff/den.diff

diff.oc.sort <- sort(diff.oc) # sort the bootstrapped differences

# draw a CDF picture
plot(x=diff.oc.sort, y=1:1000, pch=19, cex=0.5,
     main="Cummulative distribution of differences", xlab="Difference", ylab="i")

# find the order of the original coeffition among all bootstrapped coeffitions,
# and calculate z0
z0.diff <- qnorm(which(diff.oc[1000]==diff.oc.sort)*(1/1000)) # z score
# calculate the adjusted quantile
quan.diff <- c(z0.diff + (z0.diff - 1.96)/(1 - a.diff*(z0.diff - 1.96)),
               z0.diff + (z0.diff + 1.96)/(1 - a.diff*(z0.diff + 1.96)))

lower.a.diff <- (diff.oc.sort[close.p(x=pnorm(quan.diff[1]))[1]] + 
                     diff.oc.sort[close.p(x=pnorm(quan.diff[1]))[2]])/2
upper.a.diff <- (diff.oc.sort[close.p(x=pnorm(quan.diff[2]))[1]] + 
                     diff.oc.sort[close.p(x=pnorm(quan.diff[2]))[2]])/2
c(lower.a.diff, upper.a.diff) # 95% confidence interval under BCa
# 0 is not in the interval, so we reject the null hypothesis

# summarize the CI under 2 methods
c(sort(diff.oc)[1000*0.025], sort(diff.oc)[1000*0.975])
c(lower.a.diff, upper.a.diff) # 95% confidence interval under BCa
