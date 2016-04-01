############################################################
## AUTHOR: Melinda Ashcroft                                #
## AFFILIATIONS: University of Queensland - SCMB           #
## STRATEGY: Plot a linear model of the root-to-tip        #
## distance of a Maximum Likelihood phylogenetic tree      #
## against date of strain isolation and colour according   #
## to phylogenetic clade                                   #
## INPUT: Takes a tab delimited file (export as data from  #
## Path-O-Gen) of strain name, isolation date, root-to-tip #
## distance and residual. Clade column was added manually  #
############################################################

# Import libraries
library(ggplot2)
library(RSvgDevice)

# Read in the data
pathogen = read.csv(file.choose(), header=T, sep='\t')

# Attach variable (column) names
attach(pathogen)

# Show the names of each column and check the data
names(pathogen)
head(pathogen)

# Fit the linear regression model and obtain a summary
fit1 <- lm(distance ~ date, data = pathogen)
fit1
fitsum <- summary(fit1)
fitsum

# Obtain r2, P value and correlation co-efficient values
r2 = fitsum$adj.r.squared       # Obtain R2 value
r2
fitsum$coefficients
pval = fitsum$coefficients[2,4] # Obtain p value
pval
cor(distance,date)              # Obtain correlation co-efficient

###############

# Basic plot - no colour
plot(distance ~ date, data = pathogen)
abline(fit1)

###############

# Set up ggplot with colour and save to sc.plot variable
sc.plot <- ggplot(pathogen, aes(x=date, y=distance)) +  geom_point(aes(color = factor(clade)), shape=1) + xlim(1935, 2012) + ylim(-0.001, 0.04) +
  xlab("Strain Isolation Date (Years)") +
  ylab("Genetic Distance (SNPs/Site)") +
  scale_color_manual(values=c("#FF8601", "#008700"), name = "Clade", labels=c("Clade B", "Clade C")) +
  geom_smooth(method="lm", se=FALSE, fullrange=TRUE, color = "black") +
  theme_bw(base_size = 8, base_family = "Helvetica") + theme(panel.border = element_rect(color = "dark grey", size = .3), 
                                                             panel.grid = element_blank(),
                                                             axis.line = element_line(size=.7, color = "black"),
                                                             axis.text = element_text(size=14),
                                                             legend.position = c(0.9, 0.2),
                                                             legend.text = element_text(size=14),
                                                             legend.title = element_text(size=14),
                                                             text = element_text(size=14))
# Plot sc.plot variable
sc.plot

# Plot Image
# Save image as an svg for journal publication (width = 2 columns)
# Change dimensions as required
# Dimensions: width = 6.88 inches, height = 3.98 inches
devSVG(file="molecular_clock.svg", width=6.88, height=3.98) # Open plotting device
  sc.plot
dev.off() # Turn off plotting device
