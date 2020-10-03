library("xlsx")

# This is the results from Smith et al., 2020
metaSummary<-read.xlsx("data/media-2.xlsx", sheetName = "Sup Table 1", header = TRUE, startRow = 3)  #prefrontal cortex
dim(metaSummary) #365  23
metaSummary<-metaSummary[is.na(metaSummary$NA..1)==F,]
dim(metaSummary) #236  23

metaTG<-read.xlsx("data/media-2.xlsx", sheetName = "Sup Table 3", header = TRUE, startRow = 3)  #TG
dim(metaTG) #225  26
metaTG<-metaTG[is.na(metaTG$NA..1)==F,]
dim(metaTG) #95 26

metaEC<-read.xlsx("data/media-2.xlsx", sheetName = "Sup Table 5", header = TRUE, startRow = 3)  #EC
dim(metaEC) #86  26
metaEC<-metaEC[is.na(metaEC$NA..1)==F,]
dim(metaEC) #10 26

meta<-read.xlsx("data/media-2.xlsx", sheetName = "Sup Table 7", header = TRUE, startRow = 3)  #cross-cortex
dim(meta) #998  31
meta<-meta[is.na(meta$NA..1)==F,]
dim(meta) #220  31

load("bacon_correct_5SVA_robust.RData") #this is the one in final manuscript
DMP.Braak.stage.robust.corrected[DMP.Braak.stage.robust.corrected$Name=="cg26263477",] #1.209804e-11 matched to the manuscript
DMP.Braak.stage.IFG.robust.corrected[DMP.Braak.stage.IFG.robust.corrected$Name=="cg26263477",] #0.833385

Braak.stage.STG.robust.corrected.m.IFG<-merge(DMP.Braak.stage.robust.corrected, DMP.Braak.stage.IFG.robust.corrected, by = "Name", all = T, suffix = c(".STG",".IFG")) 
dim(Braak.stage.STG.robust.corrected.m.IFG) #741941     94
top.Braak.stage<-Braak.stage.STG.robust.corrected.m.IFG[(is.na(Braak.stage.STG.robust.corrected.m.IFG$P.Value.STG)==F&
                                                           Braak.stage.STG.robust.corrected.m.IFG$P.Value.STG<6.79e-08)|
                                                          (is.na(Braak.stage.STG.robust.corrected.m.IFG$P.Value.IFG)==F&
                                                             Braak.stage.STG.robust.corrected.m.IFG$P.Value.IFG<6.79e-08) ,]
dim(top.Braak.stage) #19 94
cor.test(top.Braak.stage$logFC.STG, top.Braak.stage$logFC.IFG) #0.5029725 p-value = 0.02816

Braak.stage.STG.n.meta <-merge(Braak.stage.STG.robust.corrected.m.IFG, meta, by.x = "Name", by.y = "NA..1")
dim(Braak.stage.STG.n.meta) #203 124
Braak.stage.STG.n.meta$analysis<-"Cross Cortex"
plot(Braak.stage.STG.n.meta$logFC.STG, Braak.stage.STG.n.meta$ES)
cor.test(Braak.stage.STG.n.meta$logFC.STG, Braak.stage.STG.n.meta$ES) #p-value < 2.2e-16 cor=0.784987

Braak.stage.STG.n.metaSummary <-merge(Braak.stage.STG.robust.corrected.m.IFG, metaSummary, by.x = "Name", by.y = "NA..1")
dim(Braak.stage.STG.n.metaSummary) #216 116
Braak.stage.STG.n.metaSummary$analysis<-"Prefrontal Cortex"
plot(Braak.stage.STG.n.metaSummary$logFC.STG, Braak.stage.STG.n.metaSummary$ES)
cor.test(Braak.stage.STG.n.metaSummary$logFC.STG, Braak.stage.STG.n.metaSummary$ES) #0.7773728 p-value < 2.2e-16

Braak.stage.STG.n.metaTG <-merge(Braak.stage.STG.robust.corrected.m.IFG, metaTG, by.x = "Name", by.y = "NA..1")
dim(Braak.stage.STG.n.metaTG) #89 108
Braak.stage.STG.n.metaTG$analysis<-"TG"
plot(Braak.stage.STG.n.metaTG$logFC.STG, Braak.stage.STG.n.metaTG$ES)
cor.test(Braak.stage.STG.n.metaTG$logFC.STG, Braak.stage.STG.n.metaTG$ES) #0.7690421 p-value < 2.2e-16

Braak.stage.STG.n.metaEC <-merge(Braak.stage.STG.robust.corrected.m.IFG, metaEC, by.x = "Name", by.y = "NA..1")
dim(Braak.stage.STG.n.metaEC) #9 112
Braak.stage.STG.n.metaEC$analysis<-"EC"
plot(Braak.stage.STG.n.metaEC$logFC.STG, Braak.stage.STG.n.metaEC$ES)
cor.test(Braak.stage.STG.n.metaEC$logFC.STG, Braak.stage.STG.n.metaEC$ES) #0.3234379 p-value = 0.3959

Braak.stage.STG.n.metaAll<-rbind(Braak.stage.STG.n.metaEC[,c("Name","logFC.STG","logFC.IFG","ES","analysis")],
                                 Braak.stage.STG.n.metaTG[,c("Name","logFC.STG","logFC.IFG","ES","analysis")],
                                 Braak.stage.STG.n.metaSummary[,c("Name","logFC.STG","logFC.IFG","ES","analysis")],
                                 Braak.stage.STG.n.meta[,c("Name","logFC.STG","logFC.IFG","ES","analysis")])
dim(Braak.stage.STG.n.metaAll) #517   5

library("ggplot2")
library("ggpubr")
p1<-ggplot(top.Braak.stage, aes(x=logFC.STG, y=logFC.IFG)) +
  geom_point() + labs(x = expression(paste(ES[STG], " r = 0.503 p-value = 0.03")), 
                      y = expression(ES[IFG]))
p1
p2<-ggplot(Braak.stage.STG.n.metaAll, aes(x=logFC.STG, y=ES, shape=analysis, color=analysis)) +
  geom_point() + labs(x = expression(paste(ES[STG], " ", r[TG], " = 0.769  p-value < 2.2e-16")),
                      y = "ES")
p2
p3<-ggplot(Braak.stage.STG.n.metaAll, aes(x=logFC.IFG, y=ES, shape=analysis, color=analysis)) +
  geom_point() + labs(x = expression(paste(ES[IFG], " ", r[PFC], " = 0.768  p-value < 2.2e-16")), 
                      y = "ES")
p3
pdf("Figure_2_R2.pdf", width = 14)
ggarrange(p1, p2, p3,
          labels = c("A", "B","C"),
          ncol = 1, nrow = 3)
dev.off()
