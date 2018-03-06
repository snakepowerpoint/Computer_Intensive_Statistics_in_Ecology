###### data setting ###################################################
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

copepod_composition <- read.table("copepod_composition.txt", header=T)
env <- read.csv("enviANDdensity.csv", header=T)

library(vegan)
library(vegan3d) 
library(MASS) 
library(fields) 
library(scatterplot3d) 
library(mgcv) 
library(akima) 
library(ellipse)

#######################################################################
########## PCA

# check which species are dominant species
domi <- apply(copepod_composition >= 2, 1, sum)
domi_cop <- copepod_composition[which(domi > 0), ]  # find out them

# since the number of dominant species are still too large, we reduce it
# by examining whether the species is ever greater than 10% 
# in the corresponding site

domi <- apply(domi_cop >= 10, 1, sum)
domi_cop <- domi_cop[which(domi>0), ]
domi_cop <- data.frame(t(domi_cop))   # now we have 34 sites and 17 species

# denote the seasons of each sites as 1, 2, 3
gp <- ifelse(substring(rownames(domi_cop), 1, 1) == "p", 1, 0) +
    ifelse(substring(rownames(domi_cop), 1, 1) == "s", 2, 0) +
    ifelse(substring(rownames(domi_cop), 1, 1) == "w", 3, 0)

# apply PCA on the dataframe with sites putting in the row

domi.pca <- prcomp(domi_cop, scale=T)   # remove the group variable
domi.pca$x   # PC scores

screeplot(domi.pca, type="l")   # determine how many PC to retain

# Check assumptions
hist(domi_cop[, 1])   # Normality
par(mfrow=c(1, 3))
hist(domi.pca$x[, 1])   # Use PC scores to check multivariate normality
hist(domi.pca$x[, 2])   
hist(domi.pca$x[, 3])   
dev.off()

qqnorm(domi_cop[, 9])
qqline(domi_cop[, 9])

pairs(domi_cop[, 1:8])   # Linearity
# Although it seems that many of assumptions do not meet, we still can do
# PCA on our data in a sense of descriptive statistic. Note that if we
# want to do any further inference on the PCA results, we should be cautious
# about all assumptions.

domi.rda = rda(domi_cop, scale=T)   # rda

plot(domi.pca$x[, 1], domi.pca$x[, 2], xlab='PC1', ylab='PC2', col=gp)
text(domi.pca$x[, 1], domi.pca$x[, 2], pos=1, col=gp, labels=rownames(domi_cop), cex=0.8)
legend("topright", legend=c("Spring", "Summer", "Winter"), text.col=c(1,2,3), cex=0.8)

biplot(domi.pca, choices=c(1,2), col=c('blue','red'), cex=c(0.8,1))   # biplot on PCA

ordiplot(domi.rda, choices=c(1,2), type='text', scaling=2)   # ordiplot on RDA
# it seems that species prefer some sites than others

# fitting vectors
p = ordiplot(domi.rda, choices=c(1,2), type='text')
ef = envfit(domi.rda, env[, -c(1,11,12)], permu=1000, choices=c(1,2))
ef   # most of environmental variables are significant
plot(ef, p.max=0.1)
# it shows that the reason why species prefer some sites than others is
# related to the environmental conditions



########## CA & DCA

domi.ca <- cca(domi_cop)
domi.ca

ordiplot(domi.ca, choices=c(1,2), type='text', scaling=2)   # seems to be unimodal

domi.dca = decorana(domi_cop)   # detrended correspondence analysis 
domi.dca  # print results, the length 1 is greater than 3, probably to be unimodal

domi.dca = decorana(domi_cop, iweigh=1)  # downweighting rare species 
domi.dca  # print new results 

ordiplot(domi.dca)  # ordination plot 

text(domi.dca, display='species', labels=names(domi.dca$adotj), cex=0.8, pos=1, col='red')
text(domi.dca, display='sites', labels=names(domi.dca$aidot), cex=0.8, pos=1)
# We can see that the species can be partitioned into 3 groups according to seasons.
# Some species prefer Spring, others prefer Winter, and still others prefer Summer.



########## MDS

domi.dist <- dist(domi_cop, method='manhattan')
domi.mds <- cmdscale(domi.dist)

plot(domi.mds, xlab="MDS1", ylab="MDS2", col=gp)
text(domi.mds, labels=rownames(domi_cop), pos=1, cex=0.8, col=gp)
legend("topright", legend=c("Spring", "Summer", "Winter"), text.col=c(1,2,3), cex=0.8)
# the sites can be partitioned into 3 groups, which is consistent with previous exercises



########## NMDS

domi.nmds <- metaMDS(domi_cop, distance='bray', k=3, trymax=50, autotransform=F)
stressplot(domi.nmds)   # Shepard plot

ordiplot(domi.nmds)
ordiplot(domi.nmds, type='text')
# the result is quite similar to previous exercise
