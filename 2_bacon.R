BiocManager::install("bacon")
library("bacon")
library("BiocParallel")
register(MulticoreParam(1, log=TRUE))
 
list0<-getBaconCorrected(Banner.ADvsND.STG.robust,"Banner.ADvsND.STG.robust_bacon")
# sigma.0 
# 1.30417 
# mu.0 
# -0.02742883 
# [1] 29 51
bc<-list0[[1]]
Banner.ADvsND.STG.robust.corrected<-list0[[2]]
dim(Banner.ADvsND.STG.robust.corrected)
png("Banner.ADvsND.STG.robust_bacon.png")
plot(bc, type="qq")
dev.off()

list0<-getBaconCorrected(DMP.Braak.stage.robust,"DMP.Braak.stage.robust_bacon")
# sigma.0 
# 1.221936 
# mu.0 
# -0.01031648
# [1]  5 51
bc<-list0[[1]]
DMP.Braak.stage.robust.corrected<-list0[[2]]
DMP.Braak.stage.robust.corrected[1:2,]
png("DMP.Braak.stage.robust_bacon.png")
plot(bc, type="qq")
dev.off()

list0<-getBaconCorrected(IFG.CC.robust, "IFG.CC.robust_bacon")
# sigma.0 
# 1.142586 
# mu.0 
# 0.001257152 
# sigma.0 NOTE the number is not deterministic as it varies from run to run
# 1.142849 
# mu.0 
# 0.001195565 
# [1]  0 44
# sigma.0 
# 1.14214 
# mu.0 
# 0.001075769 
# [1]  0 44
bc<-list0[[1]]
IFG.CC.robust.corrected<-list0[[2]]
png("IFG.CC.robust_bacon.png")
plot(bc, type="qq")
dev.off()
    
list0<-getBaconCorrected(DMP.Braak.stage.IFG.robust, "DMP.Braak.stage.IFG.robust_bacon")
# sigma.0 
# 1.03943 
# mu.0 
# 0.008878424 
# [1] 14 44
bc<-list0[[1]]
DMP.Braak.stage.IFG.robust.corrected<-list0[[2]]
png("DMP.Braak.stage.IFG.robust_bacon.png")
plot(bc, type="qq")
dev.off()
 