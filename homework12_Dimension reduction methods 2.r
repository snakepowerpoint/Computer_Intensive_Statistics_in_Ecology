###### data setting ###################################################
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

copepod_composition <- read.table("copepod_composition.txt", header=T)
env <- read.csv("enviANDdensity.csv", header=T)

# find dominant species
domi <- apply(copepod_composition >= 2, 1, sum)  # ever greater than 2%
domi_cop <- copepod_composition[which(domi > 0), ]
domi_cop <- data.frame(t(domi_cop))  # let row be site


########## Exercise1
library(vegan)

# normalize the domi_cop data
rsum <- apply(domi_cop, 1, sum)
y.domi <- data.frame(matrix(0, nrow(domi_cop), ncol(domi_cop)))

for (i in 1:nrow(domi_cop)){
    y.domi[i, ] = domi_cop[i, ]/rsum[i]
}
rownames(y.domi) = rownames(domi_cop)

# normalize the env data
x.env <- data.frame(scale(env[, -1]))



##### (1)
y.dca = decorana(domi_cop)
y.dca  # nonlinear, hence we use CCA



##### (2)
### correlation method
round(as.dist(cor(x.env)), 2) 

# We first keep depth because it is not highly correlated with other variables 
# Then, keep the variable with the larger constrained eigenvalue
cca(y.domi ~ Temperature, data=x.env)  # 0.4769
cca(y.domi ~ Salinity, data=x.env)  # 0.2151
cca(y.domi ~ DissolvedOxygen, data=x.env)  # 0.4549
cca(y.domi ~ maxT, data=x.env)  # 0.4807
cca(y.domi ~ MaxDO, data=x.env)  # 0.3596
# Hence, we remain maxT

# Salinity is highly correlated with maxS. Since we choose maxT,
# which is highly correlated with salinity, we drop salinity and remain maxS.

cca(y.domi ~ Fluorescence, data=x.env)  # 0.2031
cca(y.domi ~ maxF, data=x.env)  # 0.2077
# Hence, we remain maxF

cca(y.domi ~ FishDensity..ind..1000m3., data=x.env)  # 0.2229
cca(y.domi ~ CopepodDensity.ind..m3., data=x.env)  # 0.2690
# Hence, we remain copepod density

# To summarize, we have variables depth, maxT, maxS, maxF, and copepod density



### variance inflation factors method
# choose the variable with smaller VIF
vif.cca(cca(y.domi ~ ., data=x.env))   
# We keep depth, fish density, copepod density, maxDO, maxS



### stepwise variable selection, (AIC) method
# Stepwise variable selection: forward
y.full = cca(y.domi ~ ., data=x.env)   # rda with full independent variables
y.red = cca(y.domi ~ 1, data=x.env)   # rda with 1 independent variable

# step-wise choose model from simple to full
y.cca = step(y.red, scope=list(lower = ~ 1, upper=formula(y.full)))
y.cca$anova
# choose maxT, copepod density, fluorescence, maxS

# Stepwise variable selection: backward
# step-wise choose model from full to simple
y.cca = step(y.full, scope=list(lower=formula(y.red), upper=formula(y.full)))
y.cca$anova
# choose maxS, depth, fluorescence, fish density, maxT

# comparing 3 methods, we keep depth, maxT, maxS, copepod density
x.new <- x.env[, c(1,6,7,11)]



##### (3)
# Constrained ordination
y.cca = cca(y.domi ~ ., data=x.new)
summary(y.cca)

# Compare differences in WA and LC sample scores
plot(y.cca, choices=c(1,2), type='points', display='lc', scaling=1)
points(y.cca, choices=c(1,2), display='wa', pch=19, scaling=1)

# Triplot
plot(y.cca, choices=c(1,2), display=c('wa','sp','bp'), scaling=1) 
plot(y.cca, choices=c(1,2), display=c('wa','sp','bp'), scaling=2)
plot(y.cca, choices=c(1,2), display=c('wa','sp','bp'), scaling=3)  # this one is better

# combine with cluster analysis
library(cluster)

dist.man = dist(y.domi, method='manhattan')
y.man.ward = hclust(dist.man, method='ward')

plot(y.man.ward, main='Wards-linkage Dendrogram', xlab='Station',
     labels=rownames(domi_cop))
rect.hclust(y.man.ward, k=3)

gp <- as.numeric(cutree(y.man.ward, k=3))  # group label

plot(y.cca, choices=c(1,2), display=c('wa','sp','bp'), scaling=3)
ordihull(y.cca, groups=gp, col='green')



##### (4)
# We can see that the results from cluster analysis is similar with
# CCA. In CCA plot, sites are separated according to species composition
# and environmental variables. Some sites are closer more than others
# because they have similar species composition or environmental conditions.



########## Exercise2
# We have biological variable "copepod density" and 
# physical variables "depth", "maxT", "maxS"

phy <- x.new[, 1:3]
bio <- x.new[, 4]

iner.total <- sum(y.cca$CCA$eig)

y.cca.bio = cca(y.domi, phy, bio)  # partial out biological effects
plot(y.cca.bio, choices=c(1,2), display=c('wa','sp','bp'), scaling=3) 
ordihull(y.cca.bio, groups=gp, col='green')
# After partial out the biological effect, it seems that the pattern becomes obscure
iner.phy <- sum(y.cca.bio$CCA$eig)  # conditional physical effects

y.cca.phy = cca(y.domi, bio, phy)  # partial out physical effects
plot(y.cca.phy, choices=c(1,2), display=c('wa','sp','bp'), scaling=3) 
ordihull(y.cca.phy, groups=gp, col='green')
# The overlapping areas are smaller than the case under partilling out biological effects
iner.bio <- sum(y.cca.phy$CCA$eig)  # conditional biological effects

iner.total
iner.phy
iner.bio
iner.total - iner.phy - iner.bio  # interactive effects

iner.total - iner.phy  # biological effects
cca(y.domi, bio)$CCA$eig  # biological effects, coincide with above results


