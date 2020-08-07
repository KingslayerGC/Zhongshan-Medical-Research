
 sampling_procedure<-function(args)
 {
  burnin_round<-3000          ### maximum sampling round allowed at burn in step
  after_burnin_round<-3000    ### maximum sampling round allowed at post burn in step
  exclude_samegene<-T

  ###======================= specify different network =============================###
  transp<-"./supporting_files/go_propogation_probality_rp_0.3.RData" 
  cat("Loading propagation probabilities...\n")
  load(transp); nodes<-colnames(pro_p)

  gene<-extract_candidate_genes(args)
  gene<-gene[(!is.na(gene$official_name)),]
  gene<-gene[is.element(gene$official_name,nodes),]
  region<-split(gene$official_name,gene$SNP)
  cat(paste(length(unique(gene$official_name)),"genes from",length(region),"loci were found with propagation probability...\n"))

  pro_p<-pro_p[,(is.element(nodes,unique(gene$official_name)))]
  extra_weight<-extra_evi(gene)

  ####================= end of loading data ===============####

  ####================= burn in step ======================####
 
  #t0<-proc.time() 

  thres<-0.01; pickup<-0; 
  num_region<-length(region); circle<-1; chosen<-NULL
  remaining<-unlist(lapply(region,function(x) sample(x,1)))
  num0<-rep(0,sum(unlist(lapply(region,length))))

  dif<-thres+1; dif_record<-NULL 
  #while(dif>thres && circle<100)
  while(dif>thres && circle<(burnin_round+1))
  {
   pickup<-pickup%%num_region+1  
   if(pickup==1)
    if(!is.null(chosen))
    {
     num1<-NULL
     for(j in 1:length(region))
      num1<-c(num1,unlist(lapply(region[[j]],function(x) sum(chosen[,j]==x)))) 
     num1<-num1+num0
     if(circle>1)
     {
      freq0<-num0/(num_region*(circle-1))
      freq1<-num1/(num_region*circle)
      dif<-(sum((freq0-freq1)^2))^0.5
      if( circle%%50==0 )
      {
       cat("Burnin sampling, sampling circle:",circle,"\n")
      }
      #dif_record<-c(dif_record,dif)
     }
     num0<-num1; chosen<-NULL; circle<-circle+1
    }

  pickup_p<-pro_p[,is.element(colnames(pro_p),remaining[-pickup])]
  pickup_p<-pickup_p[is.element(rownames(pickup_p),unlist(region[pickup])),]  

  if(exclude_samegene && !is.null(dim(pickup_p)))                              ### if there is overlap between candidate genes and
  {                                                                            ### conditional genes, exclude the same genes 
   pickup_p<-pickup_p[,!is.element(colnames(pickup_p),rownames(pickup_p))]
   if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
  }

  if(is.null(dim(pickup_p))) { pickup_p<-1; names(pickup_p)<-region[[pickup]] } else    ### when there is only one candiate gene  
   { pickup_p<-apply(pickup_p,1,sum); pickup_p<-extra_weight[match(names(pickup_p),extra_weight$gene),]$extra_weight*pickup_p } 

  if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i]<-1/length(pickup_p)     ### avoid probabilit=0 

  remaining[pickup]<-sample(names(pickup_p),1,replace=T,prob=pickup_p)
  chosen<-rbind(chosen,remaining)  
  }
  
  #proc.time()-t0 

  ###===================== end of burn in step ===============================###

  ###======================= post-burnin step ===================================###
  #t0<-proc.time()

  pickup<-0; num_region<-length(region); circle<-1; chosen<-NULL
  num0<-rep(0,sum(unlist(lapply(region,length))))
  
  joi_dis<-matrix(0,nrow=nrow(gene),ncol=nrow(gene))
  temp<-NULL;
  for(j in 1:length(region))
   temp<-c(temp,paste(names(region[j]),region[[j]],sep="_"))
  colnames(joi_dis)<-temp; rownames(joi_dis)<-temp

  thres<-0.01; dif<-thres+1
  while(dif>thres && circle<(after_burnin_round+1) )
  {
   pickup<-pickup%%num_region+1
   if(pickup==1)
    if(!is.null(chosen))  
    {
     ###================================= calculate frequency =========================###
     num1<-NULL
     for(j in 1:length(region))
      num1<-c(num1,unlist(lapply(region[[j]],function(x) sum(chosen[,j]==x))))
     num1<-num1+num0
     if(circle>1)
     {
      freq0<-num0/(num_region*(circle-1))
      freq1<-num1/(num_region*circle)
      dif<-(sum((freq0-freq1)^2))^0.5
      if( circle%%50==0 )
      {
       cat("Post-burnin sampling, sampling circle:",circle,"\n")
      }
      #dif_record<-c(dif_record,dif)
     }
     num0<-num1; circle<-circle+1; chosen<-NULL
     ###============================= end of calculating frequency =======================###
    }

   pickup_p<-pro_p[,is.element(colnames(pro_p),remaining[-pickup])]
   pickup_p<-pickup_p[is.element(rownames(pickup_p),unlist(region[pickup])),]

   if(exclude_samegene && !is.null(dim(pickup_p)))                              ### if there is overlap between candidate genes and
   {                                                                            ### conditional genes, exclude the same genes
    pickup_p<-pickup_p[,!is.element(colnames(pickup_p),rownames(pickup_p))]
    if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
   }

   if(is.null(dim(pickup_p))) { pickup_p<-1; names(pickup_p)<-region[[pickup]] } else    ### when there is only one candiate gene
    { pickup_p<-apply(pickup_p,1,sum); pickup_p<-extra_weight[match(names(pickup_p),extra_weight$gene),]$extra_weight*pickup_p } 

   if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i]<-1/length(pickup_p)     ### avoid probabilit=0

   remaining[pickup]<-sample(names(pickup_p),1,replace=T,prob=pickup_p)
   chosen<-rbind(chosen,remaining)

   ###=============================== calculating joint distribution  =====================###
   index_col<-match(paste(names(remaining[-pickup]),remaining[-pickup],sep="_"),colnames(joi_dis))  
   index_row<-match(paste(names(remaining[pickup]),remaining[pickup],sep="_"),colnames(joi_dis))  
   joi_dis[index_row,index_col]<-joi_dis[index_row,index_col]+1 
  }

  #proc.time()-t0

  ###=====================  end of post-burnin step  =======================================###

  ###===================== summarize and record the results ================================###
  freq<-cbind(unlist(region),freq1)

  region_indicator<-NULL
  gene_num<-as.numeric(lapply(region,length))
  for(i in 1:length(gene_num))
   region_indicator<-c(region_indicator,rep(names(region[i]),gene_num[i]))

  freq<-cbind(freq,region_indicator)
  colnames(freq)<-c("gene","post_prob","region")
  freq<-as.data.frame(freq,stringsAsFactors=F)
  freq[,2]<-as.numeric(freq[,2])

  output<-NULL                                      #### sort according to posterior probability
  for(i in unique(freq$region))
  {
   temp<-freq[freq$region==i,]
   output<-rbind(output,temp[order(temp$post_prob,decreasing=T),])
  }
  freq<-output


  cat("Recording the results!\n")

  if(is.na(args[3])) res_path<-"./iRIGS_result/"     else res_path<-args[3]
  if(substr(res_path,nchar(res_path),nchar(res_path))!="/")  res_path<-paste(res_path,"/",sep="")
  if(!file.exists(res_path)) dir.create(res_path,showWarnings=T,recursive=T)

  if(is.na(args[4])) res_file<-strsplit(args[1],"/")[[1]][length(strsplit(args[1],"/")[[1]])]  else res_file<-args[4]
  res_file<-paste(res_path,res_file,"_risk_genes_predicted_by_iRIGS",sep="")
  write.table(freq,res_file,quote=F,row.names=F,sep="\t")
  output
 }

