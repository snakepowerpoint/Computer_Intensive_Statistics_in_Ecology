########## data setting
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

zoo <- read.csv("VidalTvsDuration.csv", header=T)



########## Exercise1

# plot the function first
x <- seq(-100, 100, 0.1)
y <- sin(x)/x - 0.6
plot(x, y, type='l')   # symmetric
abline(h=0)   # probably there are 2 solutions

# more detailed
x <- seq(-10, 10, 0.01)
y <- sin(x)/x - 0.6
plot(x, y, type='l')
abline(h=0)   # the solutions are probably around -2 & 2

# write a function to find the roots
# original function
f1 <- function(x){
    sin(x)/x - 0.6
}

# first derivative
f2 <- function(x){
    (x*cos(x) - sin(x))/(x^2)
}

find.root <- function(x0 = 0.5, tolerance = 0.000001){
    x1 = x0 - f1(x0)/f2(x0)
    improv = abs(x1 - x0)/abs(x0)
    xn <- numeric(1)
    
    while (tolerance <= improv){
        xn = x1 - f1(x1)/f2(x1)
        improv = abs(xn - x1)/abs(x1)
        x1 = xn
    }
    print(xn)
}

find.root()
find.root(x0 = 1)
find.root(x0 = 2)   # the function find.root works well

find.root(x0 = -1)
find.root(x0 = -2)   
# since f1 is symmetric, the solutions should be the same except for sign



########## Exercise1

temp <- zoo[, 1]

# Residual function
f3 <- function(x, ...){
    a = x[1]
    alpha = x[2]
    sum((stage - a*((temp - alpha)^(-2.05)))^2)
}

# minimize residual under C2 ~ C5
C2 <- optim(c(1, 1), f3, stage=zoo[, 2])$par
C3 <- optim(c(1, 1), f3, stage=zoo[, 3])$par
C4 <- optim(c(1, 1), f3, stage=zoo[, 4])$par
C5 <- optim(c(1, 1), f3, stage=zoo[, 5])$par

# Belehradek's equation
f4 <- function(a, alpha){
    a*((temp - alpha)^(-2.05))
}

plot(temp, zoo[, 5], pch=19, cex=0.8, col='red', xlab="Temperature", ylab="Stage duration",
     ylim=c(min(zoo[, 2:5]), max(zoo[, 2:5])))
points(temp, zoo[, 4], pch=19, col='red', cex=0.8)
points(temp, zoo[, 3], pch=19, col='red', cex=0.8)
points(temp, zoo[, 2], pch=19, col='red', cex=0.8)

points(temp, f4(a=C5[1], alpha=C5[2]), type='l', col='blue')
points(temp, f4(a=C4[1], alpha=C4[2]), type='l', col='blue')
points(temp, f4(a=C3[1], alpha=C3[2]), type='l', col='blue')
points(temp, f4(a=C2[1], alpha=C2[2]), type='l', col='blue')
