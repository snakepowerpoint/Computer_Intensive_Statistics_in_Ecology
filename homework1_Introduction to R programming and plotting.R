setwd("D:/MyGitHub/Computer_Intensive_Statistics_in_Ecology/data")

##### read data file
cop_density <- read.table("cop_density.txt", header=T)
copepod_composition <- read.table("copepod_composition.txt", header=T)

##### exercise1
com_perc <- copepod_composition/100 ## transform into percentage
den <- as.data.frame(t(c(cop_density[, 1])*t(com_perc))) ## calculate density
den <- round(den, digits=4) # round the table
write.table(den, file="cop_composition_density.txt")

##### exercise2
### richness
richness <- apply(copepod_composition > 0, 2, sum) ## calculate richness
richness

### Shannon diversity index
## H = -sum(p_i*ln(p_i)), where p_i is the proportion of the ith species'
## abundance compared with total community.
## We first replace entry 0 with 0.0001
s <- as.data.frame(lapply(den, function(x){replace(x, x==0, 0.0001)}))
## calculate p_i for each species
s <- as.data.frame(apply(s, 2, function(x){x/sum(x)}))
## calculate Shannon diversity index
si <- apply(s, 2, function(x){-sum(x*log(x))})
si

## interesting relationship between richness and Shannon diversity
plot(1:34, 10*si, type='l', col='red', ylim=c(0, max(range(richness))))
points(1:34, richness, col='blue', type='l')

##### exercise3
## check species if they once became dominant species
domi <- apply(copepod_composition > 2, 1, sum)
domi_cop <- den[which(domi > 0), ] ## find out them

## recognize how many cruises belong to spring, summer or winter
table(substr(colnames(domi_cop), start=1, stop=1))

## create a table to record the result
avg_den <- data.frame(row.names=1:43)
avg_den$spring <- apply(domi_cop[, 1:10], 1, mean)
avg_den$summer <- apply(domi_cop[, 11:25], 1, mean)
avg_den$winter <- apply(domi_cop[, 26:34], 1, mean)
avg_den

