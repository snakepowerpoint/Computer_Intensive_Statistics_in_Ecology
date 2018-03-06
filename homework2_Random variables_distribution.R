##### exercise1
### (a)
N <- 20 + sqrt(10)*rnorm(n=10000)
hist(N, freq=F, main="N(20,10)")

## compare it with rnorm(10000, mean=20, sd=sqrt(10))
hist(rnorm(10000, mean=20, sd=sqrt(10)), freq=F, main="N(20,10)")

### (b)
B <- matrix(0, 40, 10000)
B <- apply(B, 2, FUN=function(x){X=runif(n=40)})
B1 <- apply(B, 2, FUN=function(x){ifelse(x>0.5, 1, 0)})
B2 <- apply(B1, 2, sum)
hist(B2, freq=F, main="Binomial(40,0.5)")

## compare it with rbinom(n=10000, size=40, prob=0.5)
hist(rbinom(n=10000, size=40, prob=0.5), freq=F, main="Binomial(40,0.5)")

##### exercise2
# a function to sample students
candidate <- function(exclude=c()){
    all <- 1:19 ## the number of students
    No.student <- subset(all, !(all %in% exclude)) ## exclude someone
    ## Create a table to record the remaining students and generate random number
    ## from uniform distribution
    candid <- data.frame(No.=No.student, Order=runif(length(No.student)))
    candid$Order <- order(candid$Order) ## Order it, we can get
    ## an order for every week's speaker
    print(candid[order(candid$Order), ], row.names=F)
}

## try the function
candidate(exclude=c())
candidate(exclude=c(2, 5))
candidate(exclude=c(7, 18))
