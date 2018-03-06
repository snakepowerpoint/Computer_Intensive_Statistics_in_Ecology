###### data setting ###################################################
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

ceta <- read.table("cetaceancorrect.dat")

x <- ifelse(ceta$V2 > 3, 1, 0)  # turn V2 into binary data
ceta1 <- data.frame(cbind(x, ceta[, 3:4]))
colnames(ceta1) <- c("x1", "y", "z")

########## Exercise1

# -log(L(par|x)) function, which should be minimized
neg.loglike <- function(par){
    b1 <- par[1]
    b2 <- par[2]
    b3 <- par[3]
    b4 <- par[4]
    sum(log(1 + exp(b1*1 + b2*x1 + b3*y + b4*y^2)) - z*(b1*1 + b2*x1 + b3*y + b4*y^2))
}

x1 <- ceta1$x1
y <- ceta1$y
z <- ceta1$z

b <- with(ceta1, optim(par=c(0, 0, 0, 0), fn=neg.loglike))  # minimization result
b$value  # minimum value of -log(L(par|x))
# or optim(par=c(0, 0, 0 ,0), fn=neg.loglike)$par

# compare with built-in functions
ceta2 <- data.frame(cbind(ceta1, ceta1$y^2))  # create y^2 data
colnames(ceta2) <- c("x1", "y", "z", "y2")

fit <- glm(z ~ x1 + y + y2, data=ceta2, family=binomial(link = "logit"))
fit  
b$par  # consistent with built-in functions

