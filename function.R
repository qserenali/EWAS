getBaconCorrected <-function(res, fig) {
  bc <- bacon(res$t)
  print(inflation(bc))
  print(bias(bc))
  pval.df<-as.data.frame(pval(bc)) 
  colnames(pval.df)<-c("bacon.corrected.pvalue")
  corrected<-cbind(res, pval.df)
  corrected$old.P.Value<-corrected$P.Value
  corrected$old.adj.P.Val<-corrected$adj.P.Val
  corrected$P.Value<-corrected$bacon.corrected.pvalue
  corrected$adj.P.Val<-p.adjust(corrected$P.Value, method = "fdr", n = length(corrected$P.Value))
  print(dim(corrected[corrected$P.Value<0.05/length(corrected$P.Value),]))
  return (list(bc, corrected))
}

getLambda <- function(p.val, stat_type) {
  if (stat_type == "PVAL")
    z = qnorm(p.val/2)
  ## calculates lambda
  lambda <- round(median(z^2, na.rm = T) / qchisq(0.5, 1), 3)
  print(lambda)
}

getCoordinate <- function(annEPIC, symbol, flank) {
  library("org.Hs.eg.db")
  library("TxDb.Hsapiens.UCSC.hg19.knownGene")
  anno<-select(org.Hs.eg.db, symbol, c("ENTREZID","GENENAME"), "SYMBOL")
  entrezid<-anno$ENTREZID
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  cols <- c("TXID", "TXCHROM","TXSTART","TXEND","TXNAME")
  transcript.list<-select(txdb, c(as.character(entrezid)), cols, keytype="GENEID")
  start<-min(transcript.list$TXSTART)
  stop<-max(transcript.list$TXEND)
  chr<-transcript.list[1,"TXCHROM"]
  begin<-start-flank
  end<-stop+flank
  annEPICsub<-annEPIC[annEPIC$chr==chr&(annEPIC$pos>=begin&annEPIC$pos<=end),]
  # rs<-c(as.character(chr), as.character(start), as.character(stop))
  return (annEPICsub)
}

ResultIntersection <- function (set, outfile) {
  #find intersection
  ecad<-intersect(sigCpGs.entrez$sig.eg, set)
gene.list<-select(org.Hs.eg.db, keys=ecad, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
result<-ADvsND.STG[ADvsND.STG$adj.P.Val<0.05&ADvsND.STG$V1 %in% gene.list$SYMBOL,]
print(dim(result))
write.table(result, file=outfile, sep = "\t", row.names = F)
}

splitGencodeBasicV12_NAME <- function (Banner.ADvsND.STG.robust) {
  library("stringr")
  a<-str_split_fixed(Banner.ADvsND.STG.robust$GencodeBasicV12_NAME, ";", 7)
  a<-as.data.frame(a)
  a[1:2,]
  class(a$V1) #factor
  Banner.ADvsND.STG.robust<-cbind(Banner.ADvsND.STG.robust,a) 
}  

MHP.plot <- function(data1, title) {
  library("gap")
  data1$chr<-as.integer(gsub("chr","",data1$chr))
  data1<-data1[order(data1[,1], data1[,2]), ]
  colors<-rep(colors()[c(276,566)],12)
  labels<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
  pdf(paste(title,".pdf",sep="")) #save as image (change setting to highest quality)
  par(las="1", xpd=FALSE, cex.axis=1, cex=0.5, xaxs="i")
  ops <- mht.control(colors=colors,logscale=TRUE)
  mhtplot(data1,ops,pch=19)
  GWS.threshold<- -log10(0.05/dim(data1)[1])
  print(GWS.threshold)
  abline(h=GWS.threshold)
  axis(2, pos=2, at=1:20, cex=0.5)
  axis(3, col = "dark red", pos=GWS.threshold, lty = 3, lwd.ticks=0, labels = FALSE )
  box()
  title(title,cex.main=1.8)
  dev.off()  
}

getCoordinate2 <- function(annEPIC, chr, pos, flank) {
  begin<-pos-flank
  end<-pos+flank
  annEPICsub<-annEPIC[annEPIC$chr==chr&(annEPIC$pos>=begin&annEPIC$pos<=end),]
  return (annEPICsub)
}

getCoordinate3 <- function(annEPIC, chr, begin, end) {
  annEPICsub<-annEPIC[annEPIC$chr==chr&(annEPIC$pos>=begin&annEPIC$pos<=end),]
  return (annEPICsub)
}

enrichment <- function (DMP.df, all, dir.name, threshold) {
  #before adding the threshold argument, it was set to 0.0001
  library("missMethyl")
  citation("missMethyl")
  library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  library("limma")
  sigCpGs <- DMP.df[DMP.df$adj.P.Val<0.05,"Name"]
  print(length(sigCpGs))
  if (length(sigCpGs) < 500) {
    sigCpGs <- DMP.df[DMP.df$P.Value<threshold,"Name"] 
  }
  print(length(sigCpGs))
  
  PathToRes<-paste(dir.name, threshold, sep = "_")
  dir.create(PathToRes)
  n.top.geneset<-200
  png(paste(PathToRes, "topGo_bin.png", sep="/"))
  par(mfrow=c(1,1))
  gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, array.type = "EPIC", plot.bias=TRUE)
  dev.off()
  
  # Top GO categories
  res.df<-topGO(gst,ontology=c("BP"), number=n.top.geneset)
  # dotplot(res.df)
  # emapplot(res.df)
  write.csv(res.df,paste(PathToRes, "topGO.csv", sep="/")) 
  
  kst <- gometh(sig.cpg=sigCpGs, all.cpg=all,array.type = "EPIC", collection="KEGG")
  
  # Top KEGG categories
  res.df<-topKEGG(kst, number=n.top.geneset)
  write.csv(res.df,paste(PathToRes, "topKEGG.csv", sep="/")) 
  
  #Broad C2CP
  # biocLite("qusage")
  library("qusage")
  db.version<-"7.1"
  db.path<-paste("~/Banner_epigenetics/references/msigdb_v", db.version, "_files_to_download_locally/msigdb_v", db.version, "_GMTs/", sep = "")
  c2cp <- read.gmt(file = paste(db.path, "c2.cp.v", db.version, ".entrez.gmt", sep = ""))
  hsac2cp <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=c2cp, array.type = "EPIC")
  
  # sigCpGs.entrz<-getMappedEntrezIDs(sig.cpg=sigCpGs, all.cpg=all, array.type = "EPIC")
  
  res.df<-topGSA(hsac2cp, number=n.top.geneset)
  write.csv(res.df,paste(PathToRes, "hsac2cp.csv", sep="/")) 
  
  c5 <- read.gmt(file = paste(db.path, "c5.all.v", db.version, ".entrez.gmt", sep = ""))
  c5bp <- read.gmt(file = paste(db.path, "c5.bp.v", db.version, ".entrez.gmt", sep = ""))
  hsac5 <- gsameth(sig.cpg=sigCpGs, all.cpg=all, array.type=c("EPIC"), collection=c5)
  hsac5bp <- gsameth(sig.cpg=sigCpGs, all.cpg=all, array.type=c("EPIC"), collection=c5bp)
  res.df<-topGSA(hsac5, number=n.top.geneset)
  write.csv(res.df,paste(PathToRes, "hsac5.csv", sep="/")) 
  topGSA(hsac5bp, number=n.top.geneset)
  write.csv(res.df,paste(PathToRes, "hsac5bp.csv", sep="/")) 
  
  #GSEA type of pathway analysis
  library("methylGSA")
  citation("methylGSA")
  print(dim(DMP.df))
  cpg.p <- structure(as.numeric(DMP.df$P.Value), names = as.character(DMP.df$Name))
  kegg_gsea_ora1 = methylRRA(cpg.pval = cpg.p, method = "GSEA", 
                             minsize = 20, maxsize = 500, array.type = "EPIC",GS.type = "KEGG")
  head(kegg_gsea_ora1, 20)
  class(kegg_gsea_ora1) #data.frame
  write.csv(kegg_gsea_ora1,paste(PathToRes, "gsea_kegg.csv", sep="/"))
  
  # png(paste(PathToRes, "kegg_gsea_barplot.png", sep="/"), width = 480, height = 960)
  pdf(paste(PathToRes, "kegg_gsea_barplot.pdf", sep="/"))
  par(mfrow=c(1,1))
  # barplot(kegg_gsea_ora1[,c("ID","Description", "Size","pvalue",  "padj")], num = 55, colorby = "padj")
  barplot(kegg_gsea_ora1[,c("ID","Description", "Size","pvalue",  "padj")], num = 30, colorby = "padj")
  dev.off()
  
  react_gsea_ora1 = methylRRA(cpg.pval = cpg.p, method = "GSEA",minsize = 20, 
                              maxsize = 500,array.type = "EPIC", GS.type = "Reactome")
  head(react_gsea_ora1, 20)
  write.csv(react_gsea_ora1,paste(PathToRes, "gsea_reactome.csv", sep="/"))
  
  # png(paste(PathToRes, "react_gsea_barplot.png", sep="/"), width = 480, height = 1440)
  # https://stackoverflow.com/questions/3595582/saving-plot-to-tiff-with-high-resolution-for-publication-in-r
  pdf(paste(PathToRes, "react_gsea_barplot.pdf", sep="/"))
  par(mfrow=c(1,1))
  barplot(react_gsea_ora1, num = 30, colorby = "padj")
  
  library("enrichplot")
  # dotplot(react_gsea_ora1)
  dev.off()

  go_gsea1 = methylRRA(cpg.pval = cpg.p, method = "GSEA", 
                             minsize = 20, maxsize = 500, array.type = "EPIC", GS.type = "GO")
  head(go_gsea1, 20)
  class(go_gsea1) #data.frame
  write.csv(go_gsea1,paste(PathToRes, "gsea_go_methylRRA.csv", sep="/"))
  
  go_gsea1 = methylRRA(cpg.pval = cpg.p, method = "GSEA", 
                       minsize = 4, maxsize = 500, array.type = "EPIC", GS.type = "GO")
  write.csv(go_gsea1,paste(PathToRes, "gsea_go_methylRRA_4_500.csv", sep="/"))
  
  go_gsea2<-methylglm(cpg.pval = cpg.p, array.type = "EPIC",  
                        GS.type = "GO", minsize = 100, maxsize = 500)
  write.csv(go_gsea2,paste(PathToRes, "gsea_go_methylglm.csv", sep="/"))
  
  go_gsea2<-methylglm(cpg.pval = cpg.p, array.type = "EPIC",  
                      GS.type = "GO", minsize = 4, maxsize = 500)
  write.csv(go_gsea2,paste(PathToRes, "gsea_go_methylglm_4_500.csv", sep="/"))

  kegg_gsea2<-methylglm(cpg.pval = cpg.p, array.type = "EPIC",  
                      GS.type = "KEGG", minsize = 20, maxsize = 500)
  write.csv(kegg_gsea2,paste(PathToRes, "gsea_kegg_methylglm.csv", sep="/"))
  
  Reactome_gsea2<-methylglm(cpg.pval = cpg.p, array.type = "EPIC",  
                      GS.type = "Reactome", minsize = 20, maxsize = 500)
  write.csv(Reactome_gsea2,paste(PathToRes, "gsea_Reactome_methylglm.csv", sep="/"))
  
  Reactome_gsea3 <- methylgometh(cpg.pval = cpg.p, array.type = "EPIC",  sig.cut = 0.0001,  
                      GS.type = "Reactome", minsize = 20, maxsize = 500)
  write.csv(Reactome_gsea3,paste(PathToRes, "gsea_Reactome_methylgometh_0001.csv", sep="/"))
  
  go_gsea3<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC", sig.cut = 0.0001,  
                      GS.type = "GO", minsize = 20, maxsize = 500) 
  write.csv(go_gsea3,paste(PathToRes, "gsea_go_methylgometh_0001.csv", sep="/"))
  
  go_gsea3<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC", sig.cut = 0.0001,  
                         GS.type = "GO", minsize = 4, maxsize = 500) # 13519 gene sets are being tested
  write.csv(go_gsea3,paste(PathToRes, "gsea_go_methylgometh_0001_4_500.csv", sep="/"))
  
  kegg_gsea3<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC", sig.cut = 0.0001,   
                           GS.type = "KEGG", minsize = 20, maxsize = 500)
  write.csv(kegg_gsea3, paste(PathToRes, "gsea_kegg_methylgometh_0001.csv", sep="/"))
  
  Reactome_gsea4 <- methylgometh(cpg.pval = cpg.p, array.type = "EPIC",   
                                 GS.type = "Reactome", minsize = 20, maxsize = 500)
  write.csv(Reactome_gsea4,paste(PathToRes, "gsea_Reactome_methylgometh_001.csv", sep="/"))
  
  go_gsea4<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC",   
                         GS.type = "GO", minsize = 20, maxsize = 500)
  write.csv(go_gsea4,paste(PathToRes, "gsea_go_methylgometh_001.csv", sep="/"))
  
  go_gsea4<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC",   
                         GS.type = "GO", minsize = 4, maxsize = 500)
  write.csv(go_gsea4,paste(PathToRes, "gsea_go_methylgometh_001_4_500.csv", sep="/"))
  
  kegg_gsea4<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC", sig.cut = 0.001,   
                           GS.type = "KEGG", minsize = 20, maxsize = 500)
  write.csv(kegg_gsea4, paste(PathToRes, "gsea_kegg_methylgometh_001.csv", sep="/"))
  
  Reactome_gsea5 <- methylgometh(cpg.pval = cpg.p, array.type = "EPIC", sig.cut = 0.01,  
                                 GS.type = "Reactome", minsize = 20, maxsize = 500)
  write.csv(Reactome_gsea5, paste(PathToRes, "gsea_Reactome_methylgometh_01.csv", sep="/"))
  
  go_gsea5<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC", sig.cut = 0.01,   
                         GS.type = "GO", minsize = 20, maxsize = 500)
  write.csv(go_gsea5, paste(PathToRes, "gsea_go_methylgometh_01.csv", sep="/"))
  
  go_gsea5<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC", sig.cut = 0.01,   
                         GS.type = "GO", minsize = 4, maxsize = 500)
  write.csv(go_gsea5, paste(PathToRes, "gsea_go_methylgometh_01_4_500.csv", sep="/"))
  
  kegg_gsea5<-methylgometh(cpg.pval = cpg.p, array.type = "EPIC", sig.cut = 0.01,   
                         GS.type = "KEGG", minsize = 20, maxsize = 500)
  write.csv(kegg_gsea5, paste(PathToRes, "gsea_kegg_methylgometh_01.csv", sep="/"))
}

createM.norm<-function(pheno, M.norm){
  #below has to be character, as factor causing wield results
  epi.group<-as.character(pheno$barcodes)
  class(epi.group) #character
  
  o<-order(epi.group)
  print(o)
  M.norm.epi<-M.norm[, epi.group]
  dim(M.norm.epi) #863718    536
  M.norm.epi<-M.norm.epi[,o]
  pheno<-pheno[o,]
  dim(pheno) #536 421
  print(all(rownames(pheno)==colnames(M.norm.epi)))
  print(all(rownames(pheno) %in% colnames(M.norm.epi)))
  return (list(pheno, M.norm.epi)) #this may be the reason for the problem
}

getDMRGenes<-function (table.name) {
  meta.DMR<-read.xlsx("data/media-2.xlsx", sheetName = table.name, header = TRUE, startRow = 2)  #Cross Cortex
  dim(meta.DMR) #231   5
  meta.DMR<-meta.DMR[is.na(meta.DMR$Position)==F,]
  print(dim(meta.DMR)) #221  5
  meta.DMR$gene1<-gsub("\\s.*","", meta.DMR$Genes.with.TSS.within.5Kb.upstream.........GREAT.annotation., perl = T)
  meta.DMR$gene<-gsub("\\s.*","", meta.DMR$Genes.with.TSS.within.1Kb.downstream..GREAT.annotation., perl = T)
  return (meta.DMR)
}

DMR.plot<-function(dmr1, pdf_filename) {
  print(range(dmr1$P.Value)) #0.005777028 0.470149637
  dmr1$color<-as.factor(sign(dmr1$logFC))
  dmr1$y<- -log10(dmr1$P.Value)
  pdf(pdf_filename)
  class(dmr1)
  dmr1[1:2,]
  ggplot(as.data.frame(dmr1), aes(x=pos, y=y, color=color)) +
    scale_y_continuous(name = "-log10(p-value)") +
    geom_hline(yintercept = -log10(5.79e-8), linetype="dotted", color = "blue", size=1) +
    scale_x_continuous(name ="pos") +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45))
  dev.off()
}
