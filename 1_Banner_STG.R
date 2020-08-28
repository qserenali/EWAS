# Banner AD data Age matched (with age between 60 and 89 samples)

# Step 1. Load and Filter data by champ.load
library("RColorBrewer")  
library("ChAMP")  

dataDirectory <- "~/Banner_epigenetics/Banner_AD/data" # set data directory

myLoad <- champ.load(dataDirectory, method="minfi", arraytype="EPIC", force=TRUE)  
champ_filtered <- row.names(myLoad$mset)  
write.table(champ_filtered, file="~/Banner_AD/Banner_ADNDonly_STGonly_60to89_filtered_probe_list.csv", sep=",", row.names=FALSE) # 将过滤后剩余的probes写出

# Step 2. Data normalization & estimate cell counts
library("minfi")  
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")  
library("IlluminaHumanMethylationEPICmanifest")  
library("RColorBrewer")  
library("missMethyl")  
library("matrixStats")  
library("wateRmelon")  

dataDirectory <- "~/Banner_epigenetics/Banner_AD/data"   
list.files(dataDirectory, recursive=TRUE)  

# get the EPIC annotation data
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)  

targets <- read.metharray.sheet(dataDirectory, pattern="Banner_ADNDonly_STGonly_60to89_samplesheet.csv")  

#read in methyarray data
rgSet <- read.metharray.exp(targets=targets, force=TRUE)  

# filter samples based on detection p values
detP <- detectionP(rgSet)
keep <- colMeans(detP) < 0.01
table(keep)
rgSet <- rgSet[,keep]
rgSet

BiocManager::install("FlowSorted.DLPFC.450k") 
library("FlowSorted.DLPFC.450k")  
pData(rgSet)$Slide <- as.numeric(pData(rgSet)$Slide)  

# estimate cell counts
cellCounts <- estimateCellCounts(rgSet, compositeCellType = "DLPFC", cellTypes=c("NeuN_neg", "NeuN_pos"))  

targets<-cbind(targets, cellCounts)  
 
mSetRaw <- preprocessRaw(rgSet)  
 
detach("package:ChAMP", unload=TRUE)
detach("package:DMRcate", unload=TRUE)
detach("package:DSS", unload=TRUE)
detach("package:bsseq", unload=TRUE)
mSetDa.dasen<-dasen(mSetRaw) #watermelon dasen normalization
 
mSetDa.ratio<-ratioConvert(mSetDa.dasen) # Converting methylation data from methylation and unmethylation channels, to ratios (Beta and M-values)

mSetDa.gr<-mapToGenome(mSetDa.ratio) # Mapping Ilumina methylation array data to the genome

# filter based on SEX
plotSex(mSetDa.gr)
predictedSex <- getSex(mSetDa.gr, cutoff = -2)$predictedSex
predictedSex[predictedSex=="F"] = 2
predictedSex[predictedSex=="M"] = 1
predictedSex <- as.numeric(predictedSex)
length(predictedSex)
length(targets$gender)
row.names(targets.m.pheno2)<-paste(targets.m.pheno2$Slide, targets.m.pheno2$Array, sep = "_")
targets.m.pheno3<-targets.m.pheno2[row.names(targets),]
table(predictedSex == targets$gender)

#only keep probes that passed Champ filtering step
Champ_filter<-read.csv("~/Banner_epigenetics/Banner_AD/Banner_ADNDonly_STGonly_60to89_filtered_probe_list.csv") #读取刚才的filter结果
dim(Champ_filter) # 736996      1

keep<-rownames(mSetDa.gr) %in% Champ_filter[,1]  
table(keep)  

mSetDaFlt <-mSetDa.gr[keep,]  

# calculate M-values for statistical analysis
DamVals <- getM(mSetDaFlt)  
DabVals <- getBeta(mSetDaFlt)  
 
#Sample methylation level correlations
res <- cor(DabVals,method = "pearson")  
round(res, 3)

mean(res[lower.tri(res)])  

# Step 3. Use SVA to account for hidden confounding factors
targets$AD <- factor(targets$DX, levels=c("AD", "ND"))  
targets$Sex <- factor(targets$gender)
targets$Age<-targets$expired_age
Braak<-targets$Braak
targets$NeuN_neg<-targets$NeuN_neg
targets$NeuN_pos<-targets$NeuN_pos
group.t <- factor(paste(targets$DX, targets$gender, sep = "."))
targets$Braak.stage<-ifelse(targets$Braak.score=="I",1,ifelse(targets$Braak.score=="II",2,
                                            ifelse(targets$Braak.score=="III",3,ifelse(targets$Braak.score=="IV",4,
                                            ifelse(targets$Braak.score=="V",5,6)))))
 
library("sva") # load packages

# include gender as a covariate
mod <- model.matrix(~ Sex+Age+NeuN_pos+AD, data=targets) 
mod0<- model.matrix(~ Sex+Age+NeuN_pos, data=targets) 

# Linear regression with Braak stage
# mod <- model.matrix(~ Sex+Age+NeuN_pos+Braak.stage, data=targets) 
# mod0<- model.matrix(~ Sex+Age+NeuN_pos, data=targets) 

dim(DamVals) #736806     127
tmp <- DamVals[!is.infinite(rowSums(DamVals)),]
dim(tmp) #736803     127
DamVals <- tmp
sva.obj <- sva(DamVals, mod, mod0, method="irw")  
sv<- as.data.frame(sva.obj$sv)  
colnames(sv) <- paste("SV_Neup", 1:sva.obj$n.sv, sep = "") # change column names
# colnames(sv) <- paste("SVB", 1:sva.obj$n.sv, sep = "") #Braak.stage 
dim(targets)  
dim(sv) 

targets <- cbind(targets, sv)  
dim(targets)  

 #Step 4. limma to detect DMPs
library("limma") 

design.d <- model.matrix(~ 0 + AD + Sex + Age + NeuN_pos + SV_Neup1+ SV_Neup2 + SV_Neup3 + SV_Neup4 + SV_Neup5, data = targets) 
dim(design.d)  #127  10
design.d[1:2,]

# design.f <- model.matrix(~ 0 + Braak.stage + Sex + Age + NeuN_pos + SVB1 + SVB2+ SVB3 + SVB4 + SVB5, data = targets)  
# dim(design.f)#127  10 Sex was coded as factor
# design.f[1:2,]

anyDuplicated(targets$DonorID) #0

fit <- lmFit(DamVals, design.d, method="robust", maxit = 5000) #time consuming

cont.matrix <- makeContrasts(
  ADvsND = ADAD - ADND,  
  levels = colnames(design.d))
dafit.d <- contrasts.fit(fit, cont.matrix)
dafit2.d <- eBayes(dafit.d)

summary(decideTests(dafit2.d, adjust.method = "fdr"))
#        ADvsND
# Down     5713
# NotSig 726013
# Up       5077 #results from lmFit(DamVals, design.d, method="robust", maxit = 5000) same as output from NaN 

# Braak.stage
# fit <- lmFit(DamVals, design.f)
# fit.braak.robust <- lmFit(DamVals, design.f, method="robust", maxit = 5000)
# dafit2.f <- eBayes(fit)
# dafit2.f.braak.robust <- eBayes(fit.braak.robust)
 
summary(decideTests(dafit2.f.braak.robust, adjust.method = "fdr"))
#        Braak.stage   Sex1   Sex2    Age NeuN_pos   SVB1   SVB2   SVB3   SVB4   SVB5
# Down          1045 221886 222323   3378     9389  57988 218775  84650 114943  82866
# NotSig      734828 125991 125678 727297   720513 631187 313344 574548 525056 561778
# Up             930 388926 388802   6128     6901  47628 204684  77605  96804  92159

annEPICSub <- annEPIC[match(rownames(DamVals),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))] # get annotation
DMP.Braak.stage.robust<-topTable(dafit2.f.braak.robust, coef="Braak.stage", adjust="fdr", num=nrow(DamVals), genelist=annEPICSub)
 
for (contrastName in colnames(cont.matrix)) {
  resToptable.temp <- topTable(dafit2.d, coef=contrastName, adjust="fdr", num=nrow(DamVals), genelist=annEPICSub)
  print(contrastName)
  print(dim(resToptable.temp[resToptable.temp$adj.P.Val<0.05,])[1])
  print(dim(resToptable.temp[resToptable.temp$adj.P.Val<0.05&resToptable.temp$logFC>=0,])[1]) #up-regulation
  print(dim(resToptable.temp[resToptable.temp$adj.P.Val<0.05&resToptable.temp$logFC<= -0,])[1])
  assign(contrastName, resToptable.temp)
} 
 
 
 

 
 
 
 