## distri-difference-generic.R
##
## Author: Nouri BEN ZAKOUR
## Affiliations: University of Queensland
##
## Script used in Ben Zakour et al., 2016 Mbio (DOI:10.1128/mBio.00347-16)
## to analyse and plot antibiotic resistance genes of a collection of strains
## from different year of isolation, clades and countries.

## current noteboook is relative to:
    # antibiotic prevalence analysis
    # metadata trend analysis (geography ~ clade ~ time)

## data used in this script is in tab-separated format, structured as follows:
    # Strain	Clade	Year	Country	Nb_res_gene

# exploratory analytical methods contained here:

    # Linear models with different interactions:
        # Nb_res_gene ~ Year
        # Nb_res_gene ~ Year * Clade
        # Nb_res_gene ~ Year * (Clade + Country)
    # ANOVA:
        # Nb_res_gene ~ Clade
        # Nb_res_gene ~ Clade*Country
    # ANOVA plotting for each:
        # residuals vs fitted
        # normal q-q plot (standardized residuals vs theoretical quantiles)
        # scale-location (Vstandardised residuals vs fitted values)
        # residuals vs leverage
    # Tukey multiple comparisons of means:
        # on ANOVA Nb_res_gene ~ Clade
        # on ANOVA Nb_res_gene ~ Clade*Country

# exploratory visualisation methods contained here:

    # Plot not nudged - Nb res gene over time colored by clade
    # Plot nudged - Nb res gene over time colored by clade
    # Plot residuals
    # Boxplot on all data - Nb_res_gene ~ Clade
    # Boxplot with India, Korea and Portugal removed - Nb_res_gene ~ Clade * Country
    # Testing ggplot output

# load file
data <- read.csv(file.choose(), header=T, sep=',', row.names=1)
summary(data)

#loading libraries
library(ggplot2)
library(reshape)

#linear models with different interactions 1/3
# Nb_res_gene ~ Year
mod <- lm(Nb_res_gene ~ Year, data)
summary.lm(mod)

#linear models with different interactions 2/3
# Nb_res_gene ~ Year * Clade
mod2 <-lm(Nb_res_gene ~ Year * Clade, data)
summary.lm(mod2)

#linear models with different interactions 3/3
# Nb_res_gene ~ Year * (Clade + Country)
mod3 <-lm(Nb_res_gene ~ Year * (Clade + Country), data)
summary.lm(mod3)

# Plot not nudged - Nb res gene over time colored by clade
dataplot <- data.frame(Nb_res_gene=data$Nb_res_gene, Year=data$Year, Clade = data$Clade)
#plot(dataplot$Year, dataplot$Nb_res_gene, col = dataplot$Clade)
plot(dataplot$Year, dataplot$Nb_res_gene, col = dataplot$Clade)

# Plot nudged - Nb res gene over time colored by clade
ey <- rnorm(nrow(data), mean = 0, sd = 0.5)
en <- rnorm(nrow(data), mean = 0, sd = 0.5)
dataplot <- data.frame(Nb_res_gene=data$Nb_res_gene + ey, Year=data$Year + en, Clade = data$Clade)
plot(dataplot$Year, dataplot$Nb_res_gene, col = dataplot$Clade)

# Plot residuals
plot(dataplot$Year, residuals(mod3))
hist(residuals(mod3),)

#-----BOXPLOTS-------#
boxplot(Nb_res_gene ~ Clade, data, col = c("red", "orange", "chartreuse3", "dark green"), xlab="Clade", ylab="Nb_res_gene")

#remove India, Korea and Portugal
data2 <- data[data$Country != 'India' & data$Country != 'Korea' & data$Country != 'Portugal',]
write.table(data2,"data2.txt")
data3 <- read.table("data2.txt", header = TRUE)
#boxplot(Nb_res_gene ~ Clade * Country, data3, col = c(2:5), xlab="Clade.Country", ylab="Nb_res_gene")
boxplot(Nb_res_gene ~ Clade * Country, data3, col = c("red", "orange", "chartreuse3", "dark green"), ylab="Nb_res_gene", las=2)

#cboxplot(Nb_res_gene ~ Country * Clade, data3, col = c(2:7), xlab="Country.Clade", ylab="Nb_res_gene")
#boxplot(Nb_res_gene ~ Country * Clade, data3, col = c(2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5), xlab="Country.Clade", ylab="Nb_res_gene")
boxplot(Nb_res_gene ~ Country * Clade, data3, col = c("red","red","red","red","red","red","orange","orange","orange","orange","orange","orange","chartreuse3","chartreuse3","chartreuse3","chartreuse3","chartreuse3","chartreuse3","dark green","dark green","dark green","dark green","dark green","dark green"), ylab="Nb_res_gene", las=2)


# Distributions of strains per country over time
#boxplot(Year ~ Country, data, col = c(2:6), xlab="Country", ylab="Year", horizontal=TRUE)
boxplot(Year ~ Country, data, col = c("darkorange2","chartreuse3","maroon","red","gold","darkorchid","darkmagenta","dodgerblue","dark green"), xlab="Year", ylab="Country", horizontal=TRUE, las=1)

# testing ggplot output
data3.m= melt(data3)
ggplot(melt(data3), aes(x=Clade, y=Country)) + geom_boxplot(aes(fill=variable))

#----ANOVA----# http://www.statmethods.net/stats/anova.html
fit <- aov(Nb_res_gene ~ Clade, data)
summary(fit)

plot(fit)

fit2a <- aov(Nb_res_gene ~ Clade*Country, data3)
summary(fit2a)

plot(fit2a)
# fit2b <- aov(Nb_res_gene ~ Country*Clade, data3) #order matters! http://afni.nimh.nih.gov/sscc/gangc/SS.html
# summary(fit2b)
# drop1(fit2a,~.,test="F")
# drop1(fit2b,~.,test="F")

TukeyHSD(fit)

TukeyHSD(fit2a)

# highlighting significant ones
Tukey_CC <- as.data.frame(TukeyHSD(fit2a)$`Clade:Country`)
Tukey_CC[Tukey_CC$`p adj`<0.05,]

