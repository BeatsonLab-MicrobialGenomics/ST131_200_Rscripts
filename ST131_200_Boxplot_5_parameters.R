###########################################################
## AUTHOR: Melinda Ashcroft                               #
## AFFILIATIONS: University of Queensland - SCMB          #
## STRATEGY: Obtain summary statistics, visualise and     #
## calculate distribution of data. Plot each variable     #
## using boxplots                                         #
## INPUT: Takes in a csv file of 5 quality metrics        #
###########################################################

# Import libraries
library(car)
library(RSvgDevice)

# Read in the data
st131 = read.csv(file.choose(), header=T, sep=',', row.names=1)
attach(st131)

# Show the names of each column and check the data
names(st131)
head(st131)

# Summary of descriptive statistics
# Unmapped bases
summary(st131$no_unmapped_Ns)    # Obtain a summary of the data
std.err <- function(x) sd(x)/sqrt(length(x)) # Function to obtain std error
std.err(st131$no_unmapped_Ns)    # Obtain standard error
sd(st131$no_unmapped_Ns)         # Obtain standard deviation
var(st131$no_unmapped_Ns)        # Obtain variance

# Uncalled bases
summary(st131$no_low_cov.mixed)  # Obtain a summary of the data
std.err(st131$no_low_cov.mixed)  # Obtain standard error
sd(st131$no_low_cov.mixed)       # Obtain standard deviation
var(st131$no_low_cov.mixed)      # Obtain variance

# Genome size
summary(st131$genome_size)       # Obtain a summary of the data
std.err(st131$genome_size)       # Obtain standard error
sd(st131$genome_size)            # Obtain standard deviation
var(st131$genome_size)           # Obtain variance

# Number of scaffolds
summary(st131$no_scaffolds)      # Obtain a summary of the data
std.err(st131$no_scaffolds)      # Obtain standard error
sd(st131$no_scaffolds)           # Obtain standard deviation
var(st131$no_scaffolds)          # Obtain variance

# Coverage
summary(st131$coverage)          # Obtain a summary of the data
std.err(st131$coverage)          # Obtain standard error
sd(st131$coverage)               # Obtain standard deviation
var(st131$coverage)              # Obtain variance

# Boxplots to identify outliers based on all 5 parameters

# Save image as an svg for journal publication (width = 2 columns)
# Change dimensions as required
# Dimensions: width = 6.88 inches, height = 3.98 inches
devSVG(file="boxplot_outliers_5parameters.svg", width=6.88, height=3.98) # Open plotting device
  par(mfrow=c(1,5), cex=0.6) # Separate plot space into 3 columns, 1 row
  Boxplot(st131$coverage, data=st131, main = "Coverage", labels=row.names(st131), cex=1)
  Boxplot(st131$no_unmapped_Ns, data=st131, main = "No. of unmapped bases", labels=row.names(st131), cex=1)
  Boxplot(st131$no_low_cov.mixed, data=st131, main = "No. of uncalled bases", labels=row.names(st131), cex=1)
  Boxplot(st131$no_scaffolds, data=st131, main = "No. of scaffolds", labels=row.names(st131), cex=1)
  Boxplot(st131$genome_size, data=st131, main = "Genome size", labels=row.names(st131), cex=1)
dev.off() # Turn off plot device
par(mfrow=c(1,1)) # Convert plot space back to 1 column, 1 row

# Boxplots to identify outliers based on the 2 parameters
# that actually identify the 14 outliers (No. of scaffolds, No. of uncalled bases)

# Save image as an svg for journal publication (width = 2 columns)
# Change dimensions as required
# Dimensions: width = 6.88 inches, height = 3.98 inches
devSVG(file="boxplot_outliers_2parameters.svg", width=6.88, height=3.98) # Open plotting device
  par(mfrow=c(1,4), cex=0.8) # Separate plot space into 3 columns, 1 row
  Boxplot(st131$no_low_cov.mixed, data= st131, main= "No. of uncalled bases", labels=row.names(st131), cex=1)
  Boxplot(st131$no_scaffolds, data=st131, main= "No. of scaffolds", labels=row.names(st131), cex=1)
dev.off() # Turn off plot device
par(mfrow=c(1,1)) # Convert plot space back to 1 column, 1 row
