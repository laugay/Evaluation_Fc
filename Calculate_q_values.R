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
  
  message(paste("Simulation",sim,"\n"))
  
  for (replic in 1:number_of_replicates){
    success <- try( load(paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) , silent=T)
    
    if(class(success)=="try-error"){
      cat(paste0("Simulation ",sim," replicate ",replic," file does not exist\n"))
    }else{
      if(replic==41 && sim=="NM7"){}
      else if(replic==58 && sim=="NM7"){}
      else if(replic==72 && sim=="NM7"){}
      else if(replic==4 && sim=="NM14"){}
      else if(replic==10 && sim=="NM14"){}
      else if(replic==35 && sim=="NM14"){}
      else if(replic==47 && sim=="NM14"){}
      else if(replic==48 && sim=="NM14"){}
      else if(replic==51 && sim=="NM14"){}
      else if(replic==4 && sim=="NM21"){}
      else if(replic==18 && sim=="NM21"){}
      else if(replic==54 && sim=="NM21"){}
      else if(replic==64 && sim=="NM21"){}
      else if(replic==82 && sim=="NM21"){}
      else if(replic==72 && sim=="SV7"){}
      else if(replic==80 && sim=="SV7"){}
      else if(replic==89 && sim=="SV7"){}
      else if(replic==97 && sim=="SV7"){}
      else if(replic==8 && sim=="SV14"){}
      else if(replic==21 && sim=="SV14"){}
      else if(replic==51 && sim=="SV14"){}
      else if(replic==63 && sim=="SV14"){}
      else if(replic==64 && sim=="SV14"){}
      else if(replic==71 && sim=="SV14"){}
      else if(replic==77 && sim=="SV14"){}
      else if(replic==84 && sim=="SV14"){}
      else if(replic==100 && sim=="SV14"){}
      else if(replic==42 && sim=="SV21"){}
      else if(replic==45 && sim=="SV21"){}
      else if(replic==50 && sim=="SV21"){}
      else if(replic==51 && sim=="SV21"){}
      else if(replic==86 && sim=="SV21"){}
      else{
        if(length(which(!is.na(SNP_list$FC_q_value)))==0)
          SNP_list$FC_q_value[!is.na(SNP_list$FC_p_value)] <- qvalue(SNP_list$FC_p_value[!is.na(SNP_list$FC_p_value)])$qvalues 
        save(SNP_list,maf,Fc,file=paste0("Simulations/",sim,"/",sim,"_",replic,"_Fstat.RData")) 
        
      }
    }
  
  }
}

replic

