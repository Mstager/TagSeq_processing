#This script tests for differential expression using the edgeR package and a previously generated read counts table

library(edgeR)

reads<-read.delim("Junco_read_counts.txt", header=TRUE, row.names=1, skip=1) #gene expression with gene names of "date" genes fixed
nrow(reads) #17249 loci

reads<-reads[,-c(1:5)] #remove non-numeric fields

##Filtering and normalization
y <- DGEList(counts=reads2) #turn it into a DGElist for edgeR
keep <- rowSums(cpm(y)>1) >= 6 #filter lowly expressed genes, ie genes that are not found in at least 6 individuals (one treatment group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) #default method weighted trimmed mean of M-values (TMM)
plotMDS(y) #outliers?
nrow(y)  #number of genes in the dataset

###Define experimental design
#treatment assignments were in a csv file
treat<-read.csv("./Treatment_assignments.csv", colClasses=c("factor","factor","factor","factor")) #group assignment data
treat <- treat[match(row.names(y$samples), treat$Combo),] #only the samples w/ data in the right order

#Here I have 8 treatment groups: 4 sampling time points (weeks 1, 2, 3, and 6) and Cold and Control treatments for each
#I'm not interestd in differences among my control groups, so I combine all control groups and test for differences between treatment groups
Treatment3 = NA
for (i in 1:nrow(treat)) {
	Treatment3[i]<-paste(as.character(treat$Treatment[i]),as.character(treat$Period[i]), sep=".")}
Treatment3[grep("Control", Treatment3)] = "Control"
Treatment3 <- relevel(as.factor(Treatment3), ref="Control") #I want the Control treatment to serve as a the reference
design3 <- model.matrix(~Treatment3, data=y$samples)

y2 <- estimateGLMCommonDisp(y, design3, verbose=TRUE) #Estimate common dispersion 
y2 <- estimateGLMTrendedDisp(y2, design3) #Estimate trended disp. for each transcript
y2 <- estimateGLMTagwiseDisp(y2, design3) #Estimate empirical Bayes dispersion for each transcript

junco.norm <- cpm(y2, log=TRUE, prior.count=2, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset

##################################	
##################################
###GLMS in edgeR
#Test for differences between temperature treatments
fit3 <- glmFit(y2, design3) #Fit GLM to read counts for each gene
sum_cold <- topTags(glmLRT(fit3, coef=2:5), n=nrow(junco.norm))[topTags(glmLRT(fit3, coef=2:5), n=nrow(junco.norm))$table$FDR<0.05,]

#calculate mean log fold change for each gene across the 4 cold treatments
mean_fc = data.frame(row.names(sum_cold), NA)
names(mean_fc) = c("Gene", "FC")
for (i in 1:nrow(sum_cold)){
	mean_fc$FC[i] = mean(as.numeric(sum_cold$table[i,1:4]))
}
mean_fc = mean_fc[mean_fc$FC > 2 | mean_fc$FC < -2,]


#Here I'm testing for increasing or decreasing expression in the Cold across my 4 sampling time points
#I do this by asking if there is a linear trend in Fold Change across the treatments?
x = data.frame(rep(NA,nrow(sum_cold)), rep(NA,nrow(sum_cold)), rep(NA,nrow(sum_cold)), rep(NA,nrow(sum_cold))) #create a dummy data frame to store these results
names(x) = c("Gene","beta","p","FDR")
for (i in 1:nrow(sum_cold)){
	if (summary(lm(as.numeric(sum_cold$table[i,1:4])~c(1,2,3,6)))$coef[2,4] < 0.05){ #fit a linear model of expression ~ acclimation time
		x$Gene[i] = as.character(row.names(sum_cold)[i]) #store gene name
		x$beta[i] = summary(lm(as.numeric(sum_cold$table[i,1:4])~c(1,2,3,6)))$coef[2,1] #store beta coefficient
		x$p[i] = summary(lm(as.numeric(sum_cold$table[i,1:4])~c(1,2,3,6)))$coef[2,4] #store p-value
		x$FDR[i] = sum_cold$table[i,8] #store FDR
	}
}
x = x[!is.na(x$Gene),]

barplot(as.numeric(sum_cold$table[row.names(sum_cold)=="SLN",][1:4]), names.arg=c(1,2,3,6), col="gray")
