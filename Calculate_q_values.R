####################################################
#
#
#
####################################################


# INSTALLING Q-VALUE PACKAGE
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue")

library(qvalue)

# read simulations description table

sim_table <- read.table("Simulations/Fc.simparams",header=T)

number_of_replicates <- 100
for (sim in sim_table$simID){
  for (replic in 1:number_of_replicates){
    success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
    
    if(class(success)=="try-error"){
      cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
    }else{
      SNP_list$FC_q_value[!is.na(SNP_list$FC_p_value)] <- qvalue(SNP_list$FC_p_value[!is.na(SNP_list$FC_p_value)])$qvalues 
      save(SNP_list,maf,file=paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) 
    }
  
  }
}


