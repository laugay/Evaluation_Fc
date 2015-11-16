#####################################################################################
#
# Script to test the "Fc outlier method" modified from Goldringer and Bataillon 2002
#
# by A. Becheler, R. Vitalis & M. Navascués
# in collaboration with L. Gay and J. Ronfort
# 
#####################################################################################

# This scripts uses SLiM (Messer 2013) for simulating the data
# http://messerlab.org/software/

# Run from command line examples:
# R --no-save --args seed4random  1234                 < Main.R 
# R --no-save --args simID        "SelfAdapt_project"  < Main.R 
# R --no-save --args Ne_only      "c(T)"               < Main.R 
# R --no-save                                          < Main.R

#---------------------------------
# FIXED SETTINGS AND DEPENDENCIES
#---------------------------------

# This script uses the following packages:
require(batch)

# high penality for scientific notation
# (necessary to input large numbers in SLiM, which does not take scientific notation as input)
options("scipen"=999)

# tools to write/read SLiM input/output
source("slim_tools.R")
source("manipulate_output.R")
# tools to calculate F statistics
source("F_stats_tools.R")

#----------------------------------------------------------------------
# PARAMETERS: DEFAULT VALUES + VALUES FROM COMMAND LINE CALL ARGUMENTS
#----------------------------------------------------------------------
source("parameter_values.R")

#------------------------
# WRITE HEAD OF LOG FILE
#------------------------
source("write_head_log.R")

#-----------
# SIMULATION
#-----------
# write SLiM input files and run SLiM for drift period
source("simulation_drift.R")

# read SLiM output files from drift period
# write SLiM input files and run SLiM for selection period
source("simulation_selection.R")

# read SLiM output files from selection period
# sample individuals and loci
source("simulation_sample_output.R")


Fst <- compute_F_ST(genotype_data)




fastPHASE_in <-character()
for (chr in 1:chr_num){
  fastPHASE_in <- c(fastPHASE_in,paste0(simID,"_",replic, "_chr_",chr,"_fastPHASE_in.txt"))
}
fastPHASE_pop           <- paste0(simID,"_",replic, "_fastPHASE_pop.txt")
mutable_name            <- paste0(simID,"_",replic, "_muttable.RData") 
Fstat_name              <- paste0(simID,"_",replic, "_Fstat.RData")
haplotypes_sampled_file <- paste0(simID,"_",replic, "_haplotypes.RData")
genotypes_sampled_file  <- paste0(simID,"_",replic, "_genotypes.RData")
    
      
      
      
      

              pop1Fis <- FIS.compute(H1=t(pop1hap1),
                                     H2=t(pop1hap2),
                                     sample_size=sample_size)
              pop2Fis <- FIS.compute(H1=t(pop2hap1),
                                     H2=t(pop2hap2),
                                     sample_size=sample_size)
              
              n1 <- nrow(pop1hap1)
              n2 <- nrow(pop2hap1)
              hap1 <- rbind(pop1hap1,pop2hap1)
              hap2 <- rbind(pop1hap2,pop2hap2)
              n     <- nrow(hap1)
              nloci <- ncol(hap1)
              M <- array( c(hap1,hap2), dim=c(n,nloci,2) )
              genotype_data <- apply(M, c(1,2), sum)
              indID1 <- seq_len(n1)
              indID2 <- seq_len(n2)
              genotype_data <- cbind( seq_len(n), c(rep(1,each=n1),rep(2,each=n2)), genotype_data )

              
              
              genotypes_sampled_file
              save(genotype_data,Fst,file=genotypes_sampled_file)
              
              

              sigma_hat_1 <- 2*pop1Fis$WC_Fis/(1+pop1Fis$WC_Fis)
              sigma_hat_2 <- 2*pop2Fis$WC_Fis/(1+pop2Fis$WC_Fis)
              
              meanFis <-  mean(pop1Fis$WC_Fis,pop2Fis$WC_Fis)
              sigma_hat <- 2*meanFis/(1+meanFis)

              if (sigma_hat<0) sigma_hat<-0
              if (sigma_hat>1) sigma_hat<-1
            
            
            # MAKES INPUT FOR fastPHASE FOR HAPLOTYPE-BASED ANALYSES (1 FILE PER CHORMOSOME)

            #for (chr in 1:length(fastPHASE_in)){
            #  if (chr==1){
            #    limit1 <- 1
            #  }else{
            #    limit1 <- chromosome_structure[(chr*2)-2,1]
            #  }
            #  limit2 <- chromosome_structure[(chr*2)-1,1]
            #  
            #  loci <- setdiff(intersect(which(SNP_table$x>=limit1),which(SNP_table$x<=limit2)),m2$x)
            #  
            #  loci_position_and_name <- matrix(loci,nrow=length(loci),ncol=1,dimnames=list(SNP_table$x[loci],"loci"))
            #  
            #  write(sample_size*2,fastPHASE_in[chr],ncolumns=1)
            #  write(length(loci),fastPHASE_in[chr],ncolumns=1,append=T)
            #  write( c("P",SNP_table$x[loci]),fastPHASE_in[chr],ncolumns=length(loci)+1,append=T )
            #  
            #  for (i in 1:sample_size){
            #    write(paste0("# ind",i,"pop1"),fastPHASE_in[chr],append=T)
    #  
     #           haplo <- array(0,length(loci))
    #            derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1_list[[(i*2)-1]],SNP_table$x[loci])),]
     #           haplo[ derived_alleles ] <- 1
      #          write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
       #         pop1hap1[i,] <- haplo
                
             #   haplo <- array(0,length(loci))
            #    derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1_list[[(i*2)]],SNP_table$x[loci])),]
            #    haplo[ derived_alleles ] <- 1
            #    write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
            #    pop1hap2[i,] <- haplo
                
                
                
            #  }
            #  for (i in 1:sample_size){
            #    write(paste0("# ind",i,"pop2"),fastPHASE_in[chr],append=T)
                
            #    haplo <- array(0,length(loci))
            #    derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2_list[[(i*2)-1]],SNP_table$x[loci])),]
            #    haplo[ derived_alleles ] <- 1
            #    write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
            #    pop2hap1[i,] <- haplo
                
            #    haplo <- array(0,length(loci))
            #    derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2_list[[(i*2)]],SNP_table$x[loci])),]
            #    haplo[ derived_alleles ] <- 1
            #    write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
            #    pop2hap2[i,] <- haplo
            #  }
            #}
            #write(c(rep(1,sample_size),rep(2,sample_size)),ncolumns =sample_size*2,fastPHASE_pop)
            
             
            
            
            # ./fastPHASE -oTEST -utest_fastPHASE_pop.txt -i -Z -H200 -T1 -K8 -Pzp test_fastPHASE_in.txt
            
            
            
            
            # define sample of SNP loci that fulfil MAF criterion
            nx <- SNP_table[,"derived.x"] + SNP_table[,"ancestral.x"]
            px <- SNP_table[,"derived.x"] / nx
            ny <- SNP_table[,"derived.y"] + SNP_table[,"ancestral.y"]
            py <- SNP_table[,"derived.y"] / ny
            maf <- ( ((px + py) / 2 >= MAF_threshold) & ((px + py) / 2 <= 1-MAF_threshold) )
            
            if (!quiet) cat(paste(Sys.time(),nrow(SNP_table)-length(which(maf)),"of",nrow(SNP_table),"loci do not fulfil MAF criterion in sample for simulation",simID,"replicate",replic, "\n"))
          
            save(SNP_table,maf,file=mutable_name)
            
          
            # END ANALYSIS
            # analyse with Fc
            Fstat <- FstatFun_from_dataframe(SNP_table=SNP_table,
                                           maf=maf,
                                           dT=selection_period_duration,
                                           nbsimul=num_of_sim_test,
                                           MAF_threshold=MAF_threshold)
            
      
            #head(Fstat$Fc_list)
            
            #m2
            #head(SNP_table)
            
            chromosome_under_selection          <- ceiling( m2$x / (genome_length / chr_num) )
            first_bp_chromosome_under_selection <- (chromosome_under_selection-1) * (genome_length / chr_num)+1
            last_bp_chromosome_under_selection  <- (chromosome_under_selection) * (genome_length / chr_num)
            
            distance_bp <- array(NA,nrow(SNP_table))
          
            for (locus in intersect(which(SNP_table$x>=first_bp_chromosome_under_selection),which(SNP_table$x<=last_bp_chromosome_under_selection)) ){
              distance_bp[locus] <- abs(m2$x - SNP_table[locus,"x"])
            }
            distance_cM <- distance_bp * r * 100
            
            SNP_table <- cbind(SNP_table,distance_bp,distance_cM,Fstat$Fc_list,pop1Fis$Fis_all,pop2Fis$Fis_all)
            names(SNP_table) <- c(names(SNP_table)[1:19],"Fis_t0","Fis_t1")
            Fc <- Fstat$Fc_global
            Fis <- list(Fis_t0       <- pop1Fis$WC_Fis,
                        sigma_hat_t0 <- pop1Fis$sigma_hat,
                        Fis_t1       <- pop2Fis$WC_Fis,
                        sigma_hat_t1 <- pop2Fis$sigma_hat)
            
            save(SNP_table,maf,Fc,Fis,file=Fstat_name)
            save(pop1hap1,pop1hap2,pop2hap1,pop2hap2,file=haplotypes_sampled_file)
    
            files2save <- c(Fstat_name,
                            haplotypes_sampled_file,
                            genotypes_sampled_file,
                            #fastPHASE_in,
                            slim_in_drift,
                            slim_in_selection,
                            slim_log_drift,
                            slim_log_selection)
            
            files2delete <- c(slim_out_drift, 
                              slim_out_selection, 
                              slim_init_selection,
                              mutable_name)
            
            
            for (item in files2save){
              system( paste("mv",item,simID) )
            }
            for (item in files2delete){
              system( paste("rm",item) )
            }
            
            write(paste("Total number of polymorphic loci in the sample:",total_S), file=log_file,append=T)
            write(paste("Number of loci that do not fulfill MAF criterion:",nrow(SNP_table)-length(which(maf))), file=log_file,append=T)
            write(paste("Chormosome under selection:", chromosome_under_selection), file=log_file,append=T)
            
            

