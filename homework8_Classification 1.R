###### data setting ###################################################
setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

cop_density <- read.table("cop_density.txt", header=T)
copepod_composition <- read.table("copepod_composition.txt", header=T)

library(vegan)
library(cluster)
library(fpc)
library(TWIX)
#######################################################################
# check which species are dominant species
domi <- apply(copepod_composition >= 2, 1, sum)
domi_cop <- copepod_composition[which(domi > 0), ] # find out them

# transpose the data because we are interested in the patterns among stations
domi_cop <- data.frame(t(domi_cop))



##### Non-hierarchical clustering
# compute distance matrix
dist.eucl = dist(domi_cop, method='euclidean')

# Determine number of clusters
withinss = (nrow(domi_cop) - 1)*sum(apply(domi_cop, 2, var))
si = numeric(15)

for (i in 2:15){
    withinss[i] = sum(kmeans(domi_cop, centers=i)$withinss)
    si[i] = summary(silhouette(pam(domi_cop, i)))$avg.width
}

plot(1:15, withinss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 
par(new=TRUE)
plot(2:15, si[2:15], type="l", lty=2, col=2, lwd=2, xaxt="n", yaxt="n",
     xlab="", ylab="")
axis(4)
mtext("si", side=4)
legend("topright", col=c("black","red"), lty=c(1,2), pch=c(21, NA_integer_), 
       legend=c("WSS","Si"))

# compute k-means NHC around 7 mediods (clusters)
y.pam = pam(dist.eucl, k=7) 
y.clara = clara(domi_cop, k=7) 

# cluster plots of results
plot(y.pam, which.plot=1, labels=2, lines=2) 
plot(y.pam, which.plot=2)

# try 5 mediods (clusters)
y.pam.1 = pam(dist.eucl, k=5) 

plot(y.pam.1, which.plot=1, labels=2, lines=2) 
plot(y.pam.1, which.plot=2)

# try 4 mediods (clusters)
y.pam.2 = pam(dist.eucl, k=4) 

plot(y.pam.2, which.plot=1, labels=2, lines=2) 
plot(y.pam.2, which.plot=2)

# try 3 mediods (clusters)
y.pam.3 = pam(dist.eucl, k=3) 

plot(y.pam.3, which.plot=1, labels=2, lines=2) 
plot(y.pam.3, which.plot=2) # better than k=7, k=5, and k=4


##### Hierarchical clustering

y.eucl.ave = hclust(dist.eucl, method='average')
y.eucl.sin = hclust(dist.eucl, method='single')
y.eucl.com = hclust(dist.eucl, method='complete')
y.eucl.ward = hclust(dist.eucl, method='ward.D')

# alternative algorithms
y.eucl.ave = agnes(dist.eucl,method='average')
y.eucl.sin = agnes(dist.eucl,method='single')
y.eucl.com = agnes(dist.eucl,method='complete')
y.eucl.ward = agnes(dist.eucl,method='ward')

### dendrogram
# under average-linkage
plot(y.eucl.ave, main='Average-linkage Dendrogram', xlab='Station',
     labels=rownames(domi_cop))
rect.hclust(y.eucl.ave, k=3) # cut by #clusters

# under single-linkage
plot(y.eucl.sin, main='Single-linkage Dendrogram', xlab='Station',
     labels=rownames(domi_cop))
rect.hclust(y.eucl.sin, k=3) # cut by #clusters

# under complete-linkage
plot(y.eucl.com, main='Complete-linkage Dendrogram', xlab='Station',
     labels=rownames(domi_cop))
rect.hclust(y.eucl.com, k=3) # cut by #clusters

# under wards-linkage
plot(y.eucl.ward, main='Wards-linkage Dendrogram', xlab='Station',
     labels=rownames(domi_cop))
rect.hclust(y.eucl.ward, k=3) # cut by #clusters
# We can see that wards-linkage maximizes the among-group distance.
# Clearly, spring, summer, and winter are almost partitioned into 3
# different groups

# aggolmerative coefficient
summary(y.eucl.ave)$ac 
summary(y.eucl.ward)$ac 
# the value is very high, which suggests that wards-linkage performs well

# cophenetic correlation
cor(dist.eucl, cophenetic(y.eucl.ave))
cor(dist.eucl, cophenetic(y.eucl.ward))
plot(dist.eucl, cophenetic(y.eucl.ave), xlab="Original distances",
     ylab="Cophenetic distances")
plot(dist.eucl, cophenetic(y.eucl.ward), xlab="Original distances",
     ylab="Cophenetic distances")
