
#This is basically the R code from FGT_T3, FGT_T4 and FGT_T5 combined into one script.
#Running the code in RStudio or MobaXterm is easy-load and run or source the file on one of the servers
#Modifying the code requires that you understand the R script and the FGT classes!
#Note that this code often uses the input files directly from /shared_files.  If you do not intend to edit a copy is not needed
#We do not need to load data between the tutorials - we just use the same data across the whole script

#Problems (features that might be improved)
#
#The code has no interactivity or feedback to the user
#The code has some duplication that could be addressed using a better design!
#Comments are messy and the code layout is poor in places
#Note that there are enhancements given in class that are not in the skeleton such as 
#SimpleAffy, PMA calls, targets file use etc
#Most of the missing features would results in better outputs and more sustainable code
#The code will plot figures one after another, OK in RStudio but for SSH sessions they should be saved

#But the skeleton code does run and reproduce the class practical!!


# Affymetrix Microarray Analysis Basic (Skeleton) Workflow 
##Load the required libraries & load the files for the workflow
library(limma)
library(affy)
library(annotate)
library(mouse4302.db)# load chip-specific annotation

#install.packages("scatterplot3d",repo="http://cran.ma.imperia#l.ac.uk")
#Then load the library
library(scatterplot3d)

## Load the main data- commented code is just for information
system("tar -xvf /shared_files/FGT_T3_GSE10806_RAW.tar")
system("cp  /shared_files/FGT_T3_targets.txt .")
# Load the target file into an AnnotatedDataFrame object
adf<-read.AnnotatedDataFrame("/shared_files/FGT_T3_targets.txt",header=TRUE,row.names=1,as.is=TRUE)
# Load the expression values of all the CEL files in the targets file
#mydata <- ReadAffy(filenames=pData(adf)$FileName,phenoData=adf)

# Or just to quickly load all CEL files in the R working directory
mydata <- ReadAffy()
# View a summary of the example data
mydata

## Build Quality Control Plots

# Quality control plots
hist(mydata)

# And a boxplot with different colour per sample group
colours <- c(rep("yellow",3),rep("red",2),rep("blue",2), "red", rep("green",3))

boxplot(mydata, col=colours, las=2)

## Normalise the data using RMA

eset <- rma(mydata)
eset
# To obtain a matrix of the expression values, use exprs() 
values <- exprs(eset)

## Plot Normalised Data


# Boxplot to observe the results of normalisation
# Notice differences with the boxplot from the raw data
boxplot(values, col=colours,las=2)

# MA plot of the samples 1 and 4
mva.pairs(values[, c(1,4)])
# The same plot for the non-normalised raw data
# Note that the mva.pairs call below only plots a few of the  #samples – you may wish to plot them all but this is slow
mva.pairs(pm(mydata)[, c(1,4)])

## Plot Heatmap

# To facilitate interpretation, let’s replace the columns # # header,currently
# displaying the filename, to show the name of each sample 
# (if you have a targets file)
colnames(values) <- rownames(pData(adf))
# Performs hierarchical clustering with average linkage based on
# Pearson’s Correlation Coefficient
hc<-hclust(as.dist(1-cor(values, method="pearson")), method="average")
plot(hc)
## Perform PCA

pca <- prcomp(t(values), scale=T)
# Plot the PCA results

s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=rainbow(1))
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(values),pos = 3,offset = 0.5)

## Perform fold filtering

#obtaining a matrix of expression values
exprsvals <- exprs(eset)
#RMA outputs log2 data while MAS5 outputs linear data
#To convert from log…
exprsvals10 <-2^exprsvals
#check conversion
exprsvals[1:10,]
#converted
exprsvals10[1:10,]

#More fold filtering
#check order of sample names
mysamples <- sampleNames(eset)
#display the list
mysamples
#it is useful to obtain a vector of ProbeIDs here
probesets <- probeNames(mydata)
#display the first 10 ProbeSets
probesets[1:10]

#Build final fold table
#Calculate the means
#Note mean of the log is not the same as the log of the mean!!
ES.mean <- apply(exprsvals10[,c("GSM272753.CEL", "GSM272836.CEL","GSM272837.CEL")],1,mean)
iPS_OK.mean <- apply(exprsvals10[,c("GSM272839.CEL", "GSM272846.CEL","GSM272890.CEL")],1,mean)
iPS_4F.mean <- apply(exprsvals10[,c("GSM279200.CEL", "GSM279201.CEL","GSM279202.CEL")],1,mean)
NSC.mean <- apply(exprsvals10[,c("GSM272847.CEL","GSM272848.CEL")],1,mean)
#calculate some fold changes
ES_iPS_OK <-ES.mean /iPS_OK.mean
ES_iPS_4F <-ES.mean /iPS_4F.mean
ES_NSC <-ES.mean /NSC.mean
#build a summary table to hold all the data
all.data= cbind(ES.mean,iPS_OK.mean,iPS_4F.mean, NSC.mean, ES_iPS_OK,
ES_iPS_4F, ES_NSC)
#check the column names
colnames(all.data)
#write the table of means as an output
write.table(all.data,file="group_means.txt", quote=F,
sep="\t",col.names=NA)

## Beginning statistical analysis

#Check original sample order
sampleNames(eset)
#Rename the samples
sampleNames(eset) <-c("ESC.1","ESC.2","ESC.3","iPS2.2","iPS2.3","NSC.1","NSC.2","iPS2.1","iPS4.1","iPS4.2","iPS.3")
#Check the samples have renamed
sampleNames(eset)


##Building annotation for differential gene identification
#establish annotation for MOE430v2
#which annotation do we need
#modified from #http://gettinggeneticsdone.blogspot.co.uk/2012/01/annotating-limma-#results-with-gene.html

eset@annotation


#packages in the annotation package
ls("package:mouse4302.db")

#build an annotation table
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA #fix padding with NA characters 
#assign as feature data of the current Eset
fData(eset) <- tmp

## Statistical analysis using Limma

#Build the design matrix
design <- model.matrix(~-1+factor(c(1,1,1,2,2,3,3,2,4,4,4)))
colnames(design) <- c("ESC","iPS2","NSC","iPS4")
#Check it makes sense
sampleNames(eset)
#output the design matrix
design

#This instructs Limma which comparisons to make
contrastmatrix <- makeContrasts(ESC-iPS2,ESC-NSC,ESC-iPS4,
levels=design)
contrastmatrix

#issue these commands to fit the model
#and make the contrasts
fit <- lmFit(eset, design)

fit2 <- contrasts.fit(fit, contrastmatrix)

#this last part essentially moderates the t-statistic using 
#the borrowed variance approach described in class
fit2 <- eBayes(fit2)

#get the results
topTable(fit2,coef=1,adjust="fdr")
myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(eset))
write.table(myresults,"myresults.txt")

#make a venn diagram
clas <- classifyTestsF(fit2)
vennDiagram(clas)


## Carry out Functional Enrichment analysis

Mm.H <- readRDS("/shared_files/MSigDB/Mm.h.all.v7.1.entrez.rds") 

#Check that you have the required objects
ls()

#Show the full contents of the annotation package
ls("package:mouse4302.db")

#Show the annotation keys in this database
keytypes(mouse4302.db) 

sampleNames(eset)
## Process annotation for functional enrichment

#Here we select from the annotation a number of keys with the primary key being PROBEID
res <- select(mouse4302.db, keys = rownames(eset), columns = c("ENTREZID", "ENSEMBL","SYMBOL"), keytype="PROBEID")
#View the top of the table
head(res)
#find the index of each row of the expression set in the #annotation object res
idx <- match(rownames(eset), res$PROBEID)
#Use the index to set the phenotypic data in the ExpressionSet
fData(eset) <- res[idx, ]
head(fData(eset), 10)
#Find all rows that don’t have an EntrezID and remove then
eset_t<-eset[is.na(fData(eset)$ENTREZID)==0,]


## Functional Enrichment Analysis

#convert to indexes
H.indices <- ids2indices(Mm.H,fData(eset_t)$ENTREZID)
#Pick the most suitable enrichment analysis tool to find #enrichment signatures in the data and run this tool So:-

#I just run mroast here as an example- justify the selection of this method!

#if you want to run mroast
results <-mroast(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#if you want to run camera
#results <-camera(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#if you want to run romer
#results <-romer(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#View the results
results
#Use help for other parameters. Note we might decide to use #exactly the same model as our differential gene analysis for #the enrichment analysis- in this case we can extract it from #the fit
#sv <- squeezeVar(fit$sigma^2,df=fit$df.residual)

write.table(results,"enrichment.txt",sep="\t")
#You can then examine the results in “enrichment.txt”.  It is a text file.  It can be downloaded to view in a spreadsheet such as Excel.

## Session Information

sessionInfo()


