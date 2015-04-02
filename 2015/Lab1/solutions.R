library(CLL)
#data(package="CLL") # To see the data in the package.

# To load the data:
data(CLLbatch, sCLLex,disease, package="CLL")


################################################################
################################################################
################################################################
################################################################
# Question 1: Do an MA plot of the two groups against each other,
# using the preprocessed data. Add a loess line to the plot.


# First let's define the samples for each of the phenotype:
# The phenotype is contained in the data frame "disease":
head(disease)
# It contains "progres." and "stable" and one missing value ("NA")
progressive.samples <- na.omit(disease$SampleID[disease$Disease=="progres."])
stable.samples      <- na.omit(disease$SampleID[disease$Disease=="stable"])
	# The "na.omit" command removes missing values. 

# Notice that the sample names of the sCLLex dataset are not the same:
colnames(sCLLex)
	# They have ".CEL" appended.
progressive.samples <- paste0(progressive.samples, ".CEL")
stable.samples      <- paste0(stable.samples, ".CEL")
	# "paste0" is a version of "paste" where the separator is set to sep=""

# Now let's look at the processed data:
matrix.processed <- exprs(sCLLex)
progressive <- matrix.processed[, colnames(matrix.processed) %in% progressive.samples]
stable      <- matrix.processed[, colnames(matrix.processed) %in% stable.samples]


# For each group, we want the mean of the samples:
progressive.mean <- rowMeans(progressive, na.rm=TRUE)
stable.mean      <- rowMeans(stable, na.rm=TRUE)

# Let's do an MA plot:
A <- progressive.mean + stable.mean # Average
M <- progressive.mean -  stable.mean # Difference
plot(A,M, pch=20, cex=0.5)
abline(h=0, col="grey", lwd = 4) # Draw a horitontal line at 0

# Smoother scatter plot:
#smoothScatter(A,M)

# To compute the loess curve:
lowess.curve <- lowess(x = A, y = M, f = 0.1) # f controls the bandwidth
lines(lowess.curve, col = "deeppink3", lwd = 4) 









################################################################
################################################################
################################################################
################################################################
# Question 2: Do the same, but to the CLLbatch:
matrix.raw.all <- exprs(CLLbatch)
matrix.raw     <- pm(CLLbatch) # Only the perfect match probes. 

progressive.raw <- matrix.raw[, colnames(matrix.raw) %in% progressive.samples]
stable.raw      <- matrix.raw[, colnames(matrix.raw) %in% stable.samples]

progressive.raw.mean <- rowMeans(progressive.raw, na.rm=TRUE)
stable.raw.mean      <- rowMeans(stable.raw, na.rm=TRUE)


A <- log2(progressive.raw.mean) + log2(stable.raw.mean)
M <- log2(progressive.raw.mean) - log2(stable.raw.mean)


plot(A,M, pch=20, cex=0.5)
abline(h=0, col="grey", lwd = 4) # Draw a horitontal line at 0
lowess.curve <- lowess(x = A, y = M, f = 0.1) # f controls the bandwidth
lines(lowess.curve, col = "deeppink3", lwd = 4) 




################################################################
################################################################
################################################################
################################################################
# Question 3: Make boxplots and density plots of the marginal distributions
# of the probe(set) values to see how these marginal distributions
# change as a result of the preprocessing


# Boxplot of the raw data:
boxplot(log2(matrix.raw))

# Boxplot of the processed data:
boxplot(matrix.processed)

# Marginal densities of the raw data:
hist(CLLbatch)
# or
plot(density(log2(matrix.raw[,1])))
for (i in 1:ncol(matrix.raw)){
	lines(density(log2(matrix.raw[,i])), col=i)
}

# Marginal densities of the processed data:


plot(density(matrix.processed[,1]))
for (i in 1:ncol(matrix.processed)){
	lines(density(matrix.processed[,i]), col=i)
}

# Will put the two plots to each other:
par(mfrow=c(2,1))

