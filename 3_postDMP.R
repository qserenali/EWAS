# STG
load("bacon_correct_5SVA_robust.RData")
source("scripts/function.R")
library("gap")
data1<-Banner.ADvsND.STG.robust.corrected[,c("chr", "pos","P.Value")]
MHP.plot(data1, "EWAS of AD vs ND in STG")

# IFG 
data1<-IFG.CC.robust.corrected[,c("chr", "pos","P.Value")]
MHP.plot(data1, "EWAS of AD vs ND in IFG")

data1<-DMP.Braak.stage.robust.corrected[,c("chr", "pos","P.Value")] 
MHP.plot(data1, "EWAS of Braak Stage in STG")

data1<-DMP.Braak.stage.IFG.robust.corrected[,c("chr", "pos","P.Value")]
MHP.plot(data1, "EWAS of Braak Stage in IFG")

# BiocManager::install("methylGSA")
library("methylGSA")
library("missMethyl")
 
source("scripts/function.R")
enrichment(Banner.ADvsND.STG.robust.corrected, annEPIC$Name, "Banner.ADvsND.STG.robust.corrected",0.01)
enrichment(DMP.Braak.stage.robust.corrected, annEPIC$Name, "DMP.Braak.stage.robust.corrected",0.01)
enrichment(IFG.CC.robust.corrected, annEPIC$Name, "IFG.CC.robust.corrected",0.01)
enrichment(DMP.Braak.stage.IFG.robust.corrected, annEPIC$Name, "DMP.Braak.stage.IFG.robust.corrected",0.01)
 