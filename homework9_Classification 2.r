###### data setting ###################################################
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

domi_cop <- read.csv("dominant_copepod.csv", header=T)
rownames(domi_cop) <- domi_cop[, 1]
domi_cop <- domi_cop[, -1]

library(vegan)
library(boot) 
library(energy) 
library(MASS) 
library(nortest) 

#######################################################################
########## Exercise1

##### using built-in functions
# I choose MRPP method
y = domi_cop[, -44]   # drop group column
grp = domi_cop[, 44]   # group vector 
table(grp)   # check group size

y.std = scale(y)   # column standardization
y.eucl = dist(y.std, method='euclidean')   # distance matrix
y.mrpp = mrpp(y.eucl, grp)   # MRPP analysis
y.mrpp   # result shows that each group is significantly different

##### write my own functions
grp.sizes = table(grp)   # group size 13, 13, 8

# write a function for computing delta
mymrpp <- function(data){
    rand = order(runif(34, 0, 1))
    data = data[rand, ]
    y.dist = as.matrix(dist(data, diag=T, upper=T))
    delta1 = (sum(y.dist[1:13, 1:13])/2)/(13*12/2)
    delta2 = (sum(y.dist[14:26, 14:26])/2)/(13*12/2)
    delta3 = (sum(y.dist[27:34, 27:34])/2)/(8*7/2)
    delta1*13/34 + delta2*13/34 + delta3*8/34
}

delta <- matrix(0, ncol=999)
delta <- apply(delta, 2, FUN=function(x){x=mymrpp(y)})

# compute the original delta
y.sort = y[order(grp), ]
y.dist = as.matrix(dist(y.sort, diag=T, upper=T))
delta0 <- ((sum(y.dist[1:13, 1:13])/2)/(13*12/2))*13/34 +
    ((sum(y.dist[14:26, 14:26])/2)/(13*12/2))*13/34 +
    ((sum(y.dist[27:34, 27:34])/2)/(8*7/2))*8/34

delta <- c(delta, delta0)   # combine the permutated delta and the original
which(sort(delta) == delta0)/1000   # p-value



########## Exercise2

boxplot(X3~grp, data = y)   # check equal variance assumption
boxplot(X5~grp, data = y)
boxplot(X14~grp, data = y)
boxplot(X15~grp, data = y)
# it seems that the variations are different among 3 groups with respect to
# these 4 variables

grp.sizes = table(grp)   # group sample sizes 
y.sort = y[order(grp), ]   # sorted data matrix
eqdist.etest(y.sort, sizes=grp.sizes, distance=FALSE, R=999)
# 3 groups come from different distributions

pairs(y[1:8])   # check multivariate normality
# it seems that these variables don't follow multivariate normality

y.mat = as.matrix(y)
summary(lm(y.mat~as.factor(grp)))
# group 3 does not significantly different

# Derive discriminant functions
y.lda <- lda(y, grouping=grp) 
y.lda 

# Assessing the importance of the canonical functions
y.lda.pred = predict(y.lda)   # classification based on DA
y.lda.pred 

scores = y.lda.pred$x
summary(lm(scores~as.factor(grp)))
# R^2 gives squared canonical correlation coef
# we have high R^2, which means the discriminant result is good

y.table = table(grp, y.lda.pred$class)   # classification table
y.table
sum(diag(y.table))/sum(y.table)   # correct classification rate, perfect!

spy <- order(abs(y.lda$scaling[, 1]), decreasing=T)[1:4]
colnames(y)[spy] 
# X165, X80, X55, X54 are the first 4 species which are most distinct among
# 3 clusters, this results can also be verified by boxplot
boxplot(X165~grp, data = y)
boxplot(X80~grp, data = y)
boxplot(X55~grp, data = y)
boxplot(X54~grp, data = y)



########## Exercise3

library(rpart)

z = rpart(grp~., data=y, method='class', parms=list(split='gini')) 
z
summary(z)

plot(z, margin=0.1)
text(z)

1 - sum(residuals(z)^2)/sum((grp - mean(grp, na.rm=T))^2, na.rm=T)

# Tree prune
z = rpart(grp~., data=y, method='class', parms=list(split='gini'), cp=.001) 

cp1 = numeric(50)
for(i in 1:50){
    z = rpart(grp~., data=y, method='class', parms=list(split='gini'), cp=.001) 
    a = data.frame(printcp(z))
    b = which.min(a[, 4])
    cp1[i] = a[b, 1]
}
cp1

plotcp(z)

y1 = numeric(50) + 1
t = cbind(cp1, y1)
xtabs(y1~cp1, t)

zp = prune.rpart(z, 0.3333)
plot(zp, margin=0.1)
text(zp)
