 
 args<-commandArgs(TRUE)

 gibbs<-function(args)
 {
  if(length(args)==0)
   stop("No arguments!\nFor help, please type \"Rscript Gibbs.R -h\" or \"Rscript Gibbs.R --help\".\n")

  if(length(args)<2 & args[1]!="-h" & args[1]!="-help")
   stop("No enough arguments!\nFor help, please type \"Rscript Gibbs.R -h\" or \"Rscript Gibbs.R --help\".\n")

  if(args[1]=="-h" | args[1]=="--help")
  {
   cat("\nUsage: Rscript Gibbs.R SNP_file flank res_path res_pref\nPlease desinage at least the first two parameters!\n") 
   cat("\nExplanation of parameters:\n")
   cat("  --SNP_file   File name of GWAS SNPs. The file must include three columns: SNP, Chr, and Pos_hg19.\n")
   cat("  --flank      An integer indicating the flanking region (basepair) of a SNP when identifying candidate genes.\n")
   cat("               1000000 is our recommended number.\n")
   cat("  --res_path   Path of result files. Optional; \"./iRIGS_result/\" by default if not designated.\n")
   cat("  --res_pref   Prefix of the result file. Optional; SNP_file name by default if not designated.\n")
   cat("\nFor example Rscript Gibbs.R ./SNP_file/SCZ_108_loci 1000000 ./iRIGS_result_108_loci/ SCZ\n\n")
  } else
  {
   source("extract_candidate_genes.R")
   source("sampling_procedure.R")
   source("extra_evi.R") 
 
   #args<-c("./SNP_file/SCZ_20_loci_test",1000000,"./iRIGS_result_108_loci/","SCZ")
  
   output<-sampling_procedure(args)
   cat("All analysis finished!\n")
  }
 }

 gibbs(args) 
 q("no")


