
 extra_evi<-function(gene)
 {

  #id_conversion<-read.delim("/home/wangq6/cgg/Breast_Cancer_causal_gene/ensembl/Homo_sapiens.GRCh38.85.ena.tsv.gz",as.is=T)
  gene$Name<-substr(gene$Name,1,15)
  #gene$ts_id<-id_conversion[match(gene$Name,id_conversion$gene_stable_id),]$transcript_stable_id
  
 
  ###======================================= distance to TSS =========================###
 
  #cat("Collecting extra evidence: ")
 
  compute_dist<-function(gene)
  {
   index<-gene$strand=="+"
   tss_dist<-0
   tss_dist[index]<-abs(gene[index,]$start_hg19-gene[index,]$SNP_pos_hg19)
   tss_dist[!index]<-abs(gene[!index,]$end_hg19-gene[!index,]$SNP_pos_hg19)
   gene$dist<-tss_dist
   gene 
  }

  gene<-compute_dist(gene) 


  ### ============================================   adding capture Hi-C   =========================================== ####
 
  cat("Collecting extra evidence: capture Hi-C...\n")
 
  adding_caphic<-function(gene)
  {
   path<-"./supporting_files/capHiC/"
   cap4<-paste(path,"GM12878_DRE_number",sep="")
   cap4<-read.delim(cap4,as.is=T)

   gene$cap4_enhancer_no<-cap4[match(gene$Name,cap4$Name),]$cap4_enhancer_no
   gene
  }

  gene<-adding_caphic(gene)

 ###===============================  end of adding capture Hi-C  ====================================================###
 

 ###======================================= adding FANTOM5  =========================================================###

 cat("Collecting extra evidence: FANTOM5...\n")

 adding_fantom5<-function(gene)
 {
  path<-"./supporting_files/Fantom5/" 
  eep<-paste(path,"enhancer_tss_associations.bed.txt",sep="")
  eep<-read.delim(eep,as.is=T)

  eep<-eep[unlist(lapply(strsplit(eep$name,";"),function(x) length(x)))==5,]
  pos<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[1]]))
  chr<-unlist(lapply(strsplit(pos,":"),function(x)x[[1]]))
  start<-as.numeric(lapply(strsplit(unlist(lapply(strsplit(pos,":"),function(x)x[[2]])),"-"),function(x) x[[1]]))
  end<-as.numeric(lapply(strsplit(unlist(lapply(strsplit(pos,":"),function(x)x[[2]])),"-"),function(x) x[[2]]))
  fdr<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[length(x)]]))
  fdr<-as.numeric(lapply(strsplit(fdr,":"),function(x) x[[2]]))
  genes<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[length(x)-2]]))

  eep$chr<-chr; eep$start<-start; eep$end<-end; eep$fdr<-fdr; eep$gene<-genes
  eep<-eep[,13:17];  eep<-eep[eep$fdr<1,]
  eep$enhancer<-paste(eep$chr,eep$start,eep$end,sep=":")
 
  gene$fantom5_enhancer_no<-0 
  for(i in 1:nrow(gene))
  {
   temp<-eep[eep$gene==gene[i,]$official_name,]

   if(nrow(temp)>0)
    gene[i,]$fantom5_enhancer_no<-length(unique(temp$enhancer))
  }
  gene
 } 

 gene<-adding_fantom5(gene)

 ###=================================== end of adding FANTOM5  ========================================================###


 ###======================================= adding de novo mutation ===================================================###

 cat("Collecting extra evidence: De novo mutations...\n")

 adding_mutation<-function()
 {
  path<-"./supporting_files/SCZ_DNM/"

  dnm<-read.delim(paste(path,"2014_Feb_De_novo_mutation_in_SCZ_Nature_623_trios_s3.txt",sep=""),as.is=T)
  dnm<-dnm[,c("Child.phenotype","Study","Genes","Gene.annotations")]

  fromer<-read.delim(paste(path,"2014_Feb_De_novo_mutation_in_SCZ_Nature_623_trios_s2.txt",sep=""),as.is=T)
  fromer$Child.phenotype<-"SZ"; fromer$Study<-"Fromer"
  fromer<-fromer[,c("Child.phenotype","Study","Genes","Gene.annotations")]

  gir<-dnm[dnm$Study=="Girard",]
  xu<-dnm[dnm$Study=="Xu",]
  gul<-dnm[dnm$Study=="Gulsuner",]

  case<-rbind(gir,xu[xu[,1]=="SZ",],gul[gul[,1]=="SZ",],fromer)
  index2<-is.element(case$Gene.annotations,c("esplice","frameshift","nonsense","missense","codon-deletion","code-insertion"))
  case<-case[index2,]
  case
  #length(unique(case$Genes))/20000
 }

 case<-adding_mutation()

 ###======================================= end of adding de novo mutation ===================================================###

 ###======================================= adding brain HiC ===================================================###

 cat("Collecting extra evidence: Brain Hi-C...\n")

 adding_brainhic<-function(gene)
 {
  path1<-"./supporting_files/BrainHiC/S22_TSS_CP.txt"
  path2<-"./supporting_files/BrainHiC/S23_TSS_GZ.txt"

  cp<-read.delim(path1,as.is=T)
  gz<-read.delim(path2,as.is=T)
  cp$enhancer<-paste(cp$chr,cp$interacting_bin_start,cp$interacting_bin_end,sep=":")
  gz$enhancer<-paste(gz$chr,gz$interacting_bin_start,gz$interacting_bin_end,sep=":")

  gene$brain_cp<-0; gene$brain_gz<-0 
 
  for(i in 1:nrow(gene))
  {
   temp<-cp[cp$ENSGID_for_TSS==gene[i,]$Name,]
   if(nrow(temp)>0)
    gene[i,]$brain_cp<-length(unique(temp$enhancer))
  }

  for(i in 1:nrow(gene))
  {
   temp<-gz[gz$ENSGID_for_TSS==gene[i,]$Name,]
   if(nrow(temp)>0)
    gene[i,]$brain_gz<-length(unique(temp$enhancer))
  }
  gene
 }

 gene<-adding_brainhic(gene)

 ###======================================= end of adding brain HiC ===================================================###

 ###======================================= adding differential expression ================================================###

 cat("Collecting extra evidence: Differential expression...\n")

 path<-paste("./supporting_files/SCZ_DE/",
             "CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-differentialExpression-includeAncestry-DxSCZ-DE.tsv",
              sep="")
 dexpr<-read.delim(path)
 
 ###==================================== end of adding differential expression =============================================###


 ###====================================   combining p-values   =============================================###
  
 #extra<-gene[,c("dist","cap4_enhancer_no","fantom5_enhancer_no","tis_expr","brain_cp","brain_gz")]
 extra<-gene[,c("dist","cap4_enhancer_no","fantom5_enhancer_no","brain_cp","brain_gz")]
 #extra<-gene[,c("cap4_enhancer_no","fantom5_enhancer_no")]
 
 for(i in 1:ncol(extra))
  extra[,i][is.na(extra[,i])]<-median(extra[,i][!is.na(extra[,i])]) 

 sigma<-cov(extra)
 s<-svd(sigma)
 s1<- s$u %*% diag((s$d)^0.5) %*% t(s$v)
 s2<- s$u %*% diag((s$d)^(-0.5)) %*% t(s$v)

 mu<-apply(extra,2,mean)
 extra_tran<-t(s2%*%(apply(extra,1,function(x) x-mu)))

 lower_tail<-c(T,F,F,F,F)
 for(i in 1:ncol(extra_tran))
  if(!lower_tail[i]) 
   extra_tran[,i]<-(-extra_tran[,i])

 extra_p<-apply(extra_tran,2,pnorm)

 p1<-length(unique(case$Genes))/20000; p_denovo<-numeric()                  ### adding de novo p-value
 p_denovo[is.element(gene$official_name,case$Genes)]<-p1
 p_denovo[!is.element(gene$official_name,case$Genes)]<-1-p1
 extra_p<-cbind(extra_p,p_denovo)
 
 extra_p<-cbind(extra_p,dexpr[match(gene$Name,dexpr$genes),]$P.Value)       ### adding differential expression p-value

 for(i in 1:ncol(extra_p))
  extra_p[,i][is.na(extra_p[,i])]<-median(extra_p[,i][!is.na(extra_p[,i])])

 #extra_weight1<-apply(extra_p,1,function(x)  (-2*sum(log(x))))
 #extra_weight2<-apply(extra_p,1,function(x) { pchisq(-2*sum(log(x)),df=2*length(x))})
 gene$extra_weight<-(-log(apply(extra_p,1,function(x) { pchisq(-2*sum(log(x)),df=2*length(x),lower.tail=F)})))

 extra_weight<-gene[,c("official_name","extra_weight")]
 extra_weight<-extra_weight[match( unique(extra_weight$official_name),extra_weight$official_name),]
 colnames(extra_weight)[1]<-"gene"
 extra_weight
 }

