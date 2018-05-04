#' ---
#' title: "Bioinformatics class 5"
#' author: "Laura"
#' ---

#plots
x <- rnorm(1000,0)
summary(x)
#graph data
boxplot(x)
#make histogram
hist(x)

#sect.1 from lab sheet
baby <- read.table("bggn213_05_rstats/bggn213_05_rstats/weight_chart.txt", header=TRUE)
plot(baby, type="b", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab= "Age (months)", ylab= "Weight (kg)")

#sect.1b
feats <- read.table("bggn213_05_rstats/bggn213_05_rstats/feature_counts.txt", sep="\t", header=TRUE)
par(mar=c(5,11,4,2))
barplot(feats$Count, names.arg = feats$Feature, horiz = TRUE, las=2)

#sect.2
#read.delim= header=TRUE, sep='\t'
mfcount <- read.delim("bggn213_05_rstats/bggn213_05_rstats/male_female_counts.txt")
mycols <- cm.colors(nrow(mfcount))
barplot(mfcount$Count, col=mycols)

#sect.2b
expr <- read.delim("bggn213_05_rstats/bggn213_05_rstats/up_down_expression.txt")
palette(c("blue", "purple", "red"))
plot(expr$Condition1, expr$Condition2, col=expr$State)

     