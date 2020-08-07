 
 extract_candidate_genes<-function(args)
 {
  gwas_file<-args[1]; flanking<-as.numeric(args[2])
 
  gwas<-read.delim(gwas_file,as.is=T)
  if( all(is.element(c("SNP","Chr","Pos_hg19"),colnames(gwas))) )                ### check column names
  gwas<-gwas[,c("SNP","Chr","Pos_hg19")]  else stop("Wrong column names!\n")
  if(substr(gwas$Chr[1],1,3)!="chr")  gwas$Chr<-paste("chr",gwas$Chr,sep="")     ### add prefix "chr" in column "Chr"

  if( flanking<0 ) stop("Please verify the flank region number!\n")              ### check flanking region number

  gene_path<-"./supporting_files/All_human_genes"
  gene<-read.delim(gene_path,as.is=T)
 
  output<-list()
  for(i in 1:nrow(gwas))
  {
   start<-max(gwas[i,]$Pos_hg19-flanking,0)
   end<-gwas[i,]$Pos_hg19+flanking
   temp<-gene[gene$chrom==gwas[i,]$Chr,]
   index1<-temp$start_hg19<end & temp$start_hg19>start
   index2<-temp$end_hg19<end & temp$end_hg19>start
   temp<-temp[index1|index2,]  

   if(nrow(temp)>0)
   {
    temp$SNP<-gwas[i,]$SNP
    temp$SNP_chr<-gwas[i,]$Chr
    temp$SNP_pos_hg19<-gwas[i,]$Pos_hg19
    output<-rbind(output,temp)
   }
  } 
  output
 } 


