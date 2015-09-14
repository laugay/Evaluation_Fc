###################################################################################
#
# Script to test the Fc outlier method modified from Goldringer and Bataillon 2002
#
# by A. Becheler, R. Vitalis & M. Navascués
# 
###################################################################################

# This scripts uses SLiM (Messer 2013) for simulating the data
# http://messerlab.org/software/

# Run from command line examples:
# R --no-save --args seed4random 1234  < Main.R 
# R --no-save --args simID "SelfAdapt_project"  < Main.R 
# R --no-save < Main.R

# This script uses the following packages:
require(batch)

#---------------
# MAIN SETTINGS
#---------------

# Set seed for random number generation
seed4random <- 9238121 

# set working directory
#working_directory <- "/home/miguel/Work/Research/2012.SelfAdapt/Evaluation_Fc"


# high penality for scientific notation
# (necessary to input large numbers in SLiM, which does not take scientific notation as input)
options("scipen"=999)

# tools to write/read SLiM input/output
source("slim_tools.R")

# Simulation ID for file identification
simID <- "test"


#----------------------------
# PARAMETERS: DEFAULT VALUES
#----------------------------

# selfing rate
sigma <- 0.0 
# mutation rate per bp
u <- 1e-8
# recombination rate per bp
r <- 1e-8
# genome size
genome_length <- 500000000
# number of chromosomes
chr_num <- 2
# effective population size (Ne)
N <- 500
# selection coeficient
sel_coef <- 0.005
# dominance cofficient
dominance_coef <- 1
# length of pure drift period (number of times the population size)
number_of_times       <- 10 # see below
# Adaptation mode: "NM"=new mutation; "SV"=standing variation
selection_mode <- "SV" #  
# number of generations between samples
selection_period_duration <- 20
# sample size
sample_size      <- 50      # number of individuals sampled
sample_size_loci <- 10000    # number of loci sampled for demographic inference)
#threshold for mimimum allele frequency
MAF_threshold <- 0.01
# number of simulations for testing Fc
num_of_sim_Fc <- 10000
# number of replicates for each scenario
num_of_replicates <- 100
# number of attemps to get the selected allele present in tha last generation
# AND in the sample
# AND fulfilling MAF criterion
max_number_of_trials <- 200

#  Do not output progress messages
quiet <- F 

# sample diploid individuals or chromosomes
sample_diploid <- T

# gets parameter and setting values from command line (using package 'batch')
# example:
#   R --vanilla --args seed4random 104415 dominance_coef 0.5 < Main.R > Main.out
parseCommandArgs()

# Final setup lines
set.seed(seed4random)
#setwd(working_directory)

log_file <- paste0(simID,"_log.txt")

write("Evaluation of Fc statistics for detection of selection. Log file", file=log_file)


# Some population genetics quantities of interest:
theta         <- N*4*u*genome_length
Ns            <- N*sel_coef
expected_S    <- u*genome_length*4*N*sum(1/1:(sample_size-1))
Ne_corrected  <- N*(2-sigma)

# length of pure drift period
drift_period_duration <- number_of_times*N

# check that sample sizes (individuals and loci) are OK with parameters
if (N<sample_size) warning(paste(Sys.time(),"/!\\ Sample size larger than population size in simulation",simID,"\n"))
if (expected_S<sample_size_loci) warning(paste(Sys.time(),"/!\\ Number of loci to sample lower than the expected number of polymorphic sites in simulation",simID,"\n"))
# check that selection strength makes scenario not be nearly-neutral
if (Ns<=1) warning(paste(Sys.time(),"/!\\ Selection strength is too low, neutral or nearly neutral scenario in simulation",simID,"\n"))

write("DEMOGRAPHY", file=log_file,append=T)
write(paste("Population size:",N), file=log_file,append=T)
write(paste("Selfing rate:",sigma), file=log_file,append=T)
write(paste("Effective population size:",Ne_corrected), file=log_file,append=T)

write("MUTATION", file=log_file,append=T)
write(paste("Mutation rate per bp:",u), file=log_file,append=T)
write(paste("Genome length:",genome_length), file=log_file,append=T)
write(paste("Theta (genome wide):",theta), file=log_file,append=T)
write(paste("Expected number of polymorphisms:",expected_S), file=log_file,append=T)

write("RECOMBINATION", file=log_file,append=T)
write(paste("Recombination rate per bp:",r), file=log_file,append=T)
write(paste("Number of chormosomes:",chr_num), file=log_file,append=T)

write("SELECTION", file=log_file,append=T)
write(paste("Selection coefficient:",sel_coef), file=log_file,append=T)
write(paste("Dominance coefficient:",dominance_coef), file=log_file,append=T)
write(paste("Ns:",Ns), file=log_file,append=T)
if(selection_mode=="NM"){
  write("Selection acting on a new mutation", file=log_file,append=T)
}else if(selection_mode=="SV"){
  write("Selection acting on standing variation", file=log_file,append=T)  
}else{
  write("Selection mode undefined, using selection on standing variation", file=log_file,append=T)    
}

write("SAMPLE", file=log_file,append=T)
write(paste("Pre-simulation for mutatuion-drift equilibrium (neutral):",drift_period_duration), file=log_file,append=T)
write(paste("Time period between samples (selection):",selection_period_duration), file=log_file,append=T)
write(paste("Sample size (individuals):",sample_size), file=log_file,append=T)
write(paste("Sample size (loci):",sample_size_loci), file=log_file,append=T)
write(paste("Minimum allele frequency threshold:",MAF_threshold), file=log_file,append=T)
write(paste("Simulations for Fc null hypothesis:",num_of_sim_Fc), file=log_file,append=T)
write(paste("Replicates of scenario:",num_of_replicates), file=log_file,append=T)
write(paste("Seed for random number generation:",seed4random), file=log_file,append=T)

write("REPLICATES", file=log_file,append=T)




#-----------
# SIMULATION
#-----------

system( paste("mkdir",simID) )
for (replic in 1:num_of_replicates){

  write(paste("Replicate:",replic), file=log_file,append=T)
  
  if (selection_mode=="NM"){
    advantageous_allele <- "derived"
  }else{
    if (selection_mode!="SV"){
      selection_mode<-"SV"
      warning("Adaptation mode undefined, using standing variation") 
    }
    # choose advantageous allele between derived and ancestral state 
    advantageous_allele <- sample(c("derived","ancestral"),size=1)
    if (advantageous_allele=="derived"){
      sel_coef <- sel_coef
      dominance_coef <- dominance_coef
      if (!quiet) cat(paste(Sys.time(),"Derived allele is advantageous in simulation",simID,"replicate",replic,"\n"))
      write("Derived allele is advantageous", file=log_file,append=T)
    }else if (advantageous_allele=="ancestral"){
      sel_coef <- -(sel_coef/(1+sel_coef))
      dominance_coef <- 1 - dominance_coef
      if (!quiet) cat(paste(Sys.time(),"Ancestral allele is advantageous in simulation",simID,"replicate",replic,"\n"))
      write("Ancestral allele is advantageous", file=log_file,append=T)
    }
  }
  
  
  # 1. SIMULATION OF A PURE DRIFT PERIOD
  
  # Set file names
  slim_in_drift  <- paste0(simID,"_",replic, "_drift.txt")
  slim_out_drift <- paste0(simID,"_",replic, "_drift.out")
  slim_log_drift <- paste0(simID,"_",replic, "_drift.log")
  
  # Write slim input file (functions from slim_tools_2.R)
  writeMutation          (file=slim_in_drift, number_of_types=2, h=c(0.5,dominance_coef), DFE=c("f","f"), s=c(0,sel_coef), append=F, append_mutation=F)
  writeMutationRate      (file=slim_in_drift, u=u)  
  writeGenomicElement    (file=slim_in_drift, number_of_types=1, mut_type=list(c("m1","m2")), prop=list(c(1,0)))
  writeChromosome        (file=slim_in_drift, element_type="g1", start=1, end=genome_length)
  chromosome_structure <- writeRecombinationChrom(file=slim_in_drift, chr_num=chr_num, genome_length=genome_length, r=r, append=T)
  writeGenerations       (file=slim_in_drift, t=drift_period_duration, append=T)
  writeDemography        (file=slim_in_drift, type="P", time=1, pop="p1", N=N)
  writeDemography        (file=slim_in_drift, type="S", time=1, pop="p1", sigma=sigma, append_demography=T)
  writeOutput            (file=slim_in_drift, type="A", time=drift_period_duration, filename=slim_out_drift)
  writeSeed              (file=slim_in_drift, seed=round(runif(1,-2^31,2^31)) )
  
  if (selection_mode=="NM"){
    writePredeterminedMutations(file=slim_in_drift, time=drift_period_duration, mut_type="m2", x=sample(genome_length,1) )
  }
  
  if (!quiet) cat(paste(Sys.time(),"SIMULATION OF A PURE DRIFT PERIOD STARTS. Simulation:",simID,"replicate",replic, "\n"))
  
  # Run SLiM
  system(paste("./slim",slim_in_drift,">",slim_log_drift))
  # NB: SLiM executable must be in the same folder as SLiM input file (this is a requirement from SLiM)
  
  if (!quiet) cat(paste(Sys.time(),"END OF SIMULATION OF A PURE DRIFT PERIOD. Simulation:",simID,"replicate",replic, "\n"))

  advantageous_allele_not_lost <- FALSE
  counter <- 0
  while(!advantageous_allele_not_lost && counter<max_number_of_trials){
    counter <- counter+1
    
    # 2. SIMULATION OF THE PERIOD WITH SELECTION (in between samples)
    if (!quiet) cat(paste(Sys.time(),"Selection period for simulation:",simID,"replicate",replic, "attempt number",counter,"\n"))
    
    
    # Set file names
    slim_in_selection   <- paste0(simID,"_",replic, "_selection.txt")
    slim_out_selection  <- paste0(simID,"_",replic, "_selection.out")
    slim_log_selection  <- paste0(simID,"_",replic, "_selection.log")
    slim_init_selection <- paste0(simID,"_",replic, "_selection_init.txt")
    fastPHASE_in <-character()
    for (chr in 1:chr_num){
      fastPHASE_in <- c(fastPHASE_in,paste0(simID,"_",replic, "_chr_",chr,"_fastPHASE_in.txt"))
    }
    fastPHASE_pop       <- paste0(simID,"_",replic, "_fastPHASE_pop.txt")
    mutable_name        <- paste0(simID,"_",replic, "_muttable.RData") 
    Fstat_name          <- paste0(simID,"_",replic, "_Fstat.RData")
    haplotypes_sampled_file <- paste0(simID,"_",replic, "_haplotypes_sampled.RData")
    
    
    
    # Read slim output from DRIFT period
    out_drift_lines      <- readLines(con=slim_out_drift)
    out_drift_pop_line   <- which(out_drift_lines=="Populations:")
    out_drift_mut_line   <- which(out_drift_lines=="Mutations:")
    out_drift_gen_line   <- which(out_drift_lines=="Genomes:")
    out_drift_num_of_mut <- out_drift_gen_line-out_drift_mut_line-1
    
    # Get table of mutations from end of drift period
    mutation_table <- read.table(file=slim_out_drift, skip=out_drift_mut_line, nrows=out_drift_num_of_mut)
    
    if (selection_mode!="NM") {
      if (selection_mode!="SV"){
        selection_mode<-"SV"
        warning("Adaptation mode undefined, using standing variation") 
      }
      
      # SELECTION ON STANDING VARIATION (changes selection coefficient of on a random locus)
      
      # chose a random locus
      sampled_mut                   <- sample(x=out_drift_num_of_mut, size=1)
      # change selection coefficient of locus
      levels(mutation_table[,2])    <- c("m1","m2") 
      mutation_table[sampled_mut,2] <- "m2"
      mutation_table[sampled_mut,4] <- sel_coef
      mutation_table[sampled_mut,5] <- dominance_coef
    }
    
    # write initialzation file for slim (state of population at starting point of selection period)
    write(out_drift_lines[out_drift_pop_line:out_drift_mut_line], slim_init_selection)
    write.table(mutation_table, slim_init_selection, append=T, quote=F, col.names=F, row.names=F)
    write(out_drift_lines[-(1:(out_drift_gen_line-1))], slim_init_selection, append=T)
    
    # if (!quiet) cat(paste(Sys.time(),"Episode de SV slim commence, répétition n° ",sim))
    
    # write input file for slim
    writeMutation          (file=slim_in_selection, number_of_types=2, h=c(0.5,dominance_coef), DFE=c("f","f"), s=c(0,sel_coef), append=F, append_mutation=F)
    writeMutationRate      (file=slim_in_selection, u=u)  
    writeGenomicElement    (file=slim_in_selection, number_of_types=1, mut_type=list(c("m1","m2")), prop=list(c(1,0)))
    writeChromosome        (file=slim_in_selection, element_type="g1", start=1, end=genome_length)
    writeRecombinationChrom(file=slim_in_selection, chr_num=chr_num, genome_length=genome_length, r=r, append=T)  
    writeGenerations       (file=slim_in_selection, t=selection_period_duration, append=T)
    writeDemography        (file=slim_in_selection, type="S", time=1, pop="p1", sigma=sigma)
    writeOutput            (file=slim_in_selection, type="A", time=selection_period_duration, filename=slim_out_selection)
    writeOutput            (file=slim_in_selection, type="F", time=selection_period_duration, append_output=T)
    writeSeed              (file=slim_in_selection, seed=round(runif(1,-2^31,2^31)) )
    
    write("#INITIALIZATION",file=slim_in_selection,ncolumns=1,append=TRUE)
    write(slim_init_selection,file=slim_in_selection,ncolumns=1,append=TRUE)
    
    if (!quiet) cat(paste(Sys.time(),"SIMULATION OF THE PERIOD WITH SELECTION STARTS. Simulation:",simID,"replicate",replic, "\n"))
    system(paste("./slim",slim_in_selection,">",slim_log_selection))
    if (!quiet) cat(paste(Sys.time(),"END OF SIMULATION OF THE PERIOD WITH SELECTION. Simulation:",simID,"replicate",replic, "\n"))
    
    # Read slim output and log from SELECTION period
    out_selection_lines    <- readLines(con=slim_out_selection)
    out_selection_mut_line <- which(out_selection_lines=="Mutations:")
    out_selection_gen_line <- which(out_selection_lines=="Genomes:")
    
    log_selection_lines    <- readLines(con=slim_log_selection)
    log_selection_mut_line <- which(log_selection_lines=="Mutations:")
    if (length(log_selection_mut_line)>0){
      log_selection_num_mut <- length(log_selection_lines)-log_selection_mut_line
    }else{
      log_selection_num_mut <- 0
    }
    
    # VERIFYING THE SUCCESS OF SELECTION (was the advantageous allele lost by drift?)
    
    if (advantageous_allele=="derived"){
      if (log_selection_num_mut > 0) {
        if (length(grep(pattern="m2",x=out_selection_lines[(out_selection_mut_line+1):(out_selection_gen_line-1)]))>0 || length(grep(pattern="m2",x=log_selection_lines[(log_selection_mut_line+1):length(log_selection_lines)]))>0){
          advantageous_allele_not_lost <- TRUE
        }else{ advantageous_allele_not_lost <- FALSE }  
      }else if(log_selection_num_mut==0){
        if (length(grep(pattern="m2",x=out_selection_lines[(out_selection_mut_line+1):(out_selection_gen_line-1)]))>0){
          advantageous_allele_not_lost <- TRUE
        }else{ advantageous_allele_not_lost <- FALSE }
      }
    }else if (advantageous_allele=="ancestral"){
      if (log_selection_num_mut > 0) {
        if (length(grep(pattern="m2",x=log_selection_lines[(log_selection_mut_line+1):length(log_selection_lines)]))>0){
          advantageous_allele_not_lost <- FALSE
        }else{ advantageous_allele_not_lost <- TRUE }  
      }else if(log_selection_num_mut==0){
        advantageous_allele_not_lost <- TRUE
      }
    }
    
    if(!advantageous_allele_not_lost){
      if (!quiet) cat(paste(Sys.time(),"Advantageous allele was lost in simulation",simID,"replicate",replic, ": REPEATING SELECTION STEP\n"))  
    }else{
      if (!quiet) cat(paste(Sys.time(),"Advantageous allele was NOT lost in simulation",simID,"replicate",replic, "\n"))
      
      # Make a list of ALL polymorphic sites
      
      SNP_list <- Make.SNP.list(file_in_1=slim_init_selection,
                                file_in_2=slim_out_selection,
                                file_fixed=slim_log_selection,
                                populations=2,
                                sample_size=sample_size,
                                N=N,
                                sample_diploid=sample_diploid)
      
      if (!quiet) cat(paste(Sys.time(),"List of all polymorphic SNP completed for simulation",simID,"replicate",replic, "\n"))     
      
      sampled_haplotypes_1_list <- SNP_list$sampled_haplotypes_1_list
      sampled_haplotypes_2_list <- SNP_list$sampled_haplotypes_2_list
      m2                        <- SNP_list$m2
      SNP_list                  <- SNP_list$SNP_list
      
      total_S <- nrow(SNP_list)
      
      if (nrow(SNP_list)<sample_size_loci){
        if (!quiet) cat(paste(Sys.time(),"Less than",sample_size_loci,"polymorphic loci (",nrow(SNP_list),") in sample for simulation",simID,"replicate",replic, ": REPEATING SELECTION STEP\n"))  
        advantageous_allele_not_lost <-F
      }else{
        if (!quiet) cat(paste(Sys.time(),"The number of polymorphic loci in the sample (",dim(SNP_list)[1],") allows to sample",
                  sample_size_loci,"loci for demographic inference in simulation",simID,"replicate",replic, "\n"))  
    
        # reduce table of polymoprhic SNP to sample_size_loci
        
        if (length(which(SNP_list[,"type"] == "m2")) == 0) {
          if (!quiet) cat(paste(Sys.time(),"SNP under selection is absente in sample from simulation",simID,"replicate",replic, ": REPEATING SELECTION STEP\n"))
          # Sample all loci randomly
          advantageous_allele_not_lost <-F
        }else{
          if (!quiet) cat(paste(Sys.time(),"SNP under selection is present in sample from simulation",simID,"replicate",replic, "\n"))
          # Sample locus under selection and the rest randomly 
          SNP_list <- SNP_list[sort(c(sample(which(SNP_list[,"type"]=="m1"),sample_size_loci-1),which(SNP_list[,"type"]=="m2"))),]
          m2_MAF <- m2.gbmaf(SNP_list,MAF_threshold)      
          if (!m2_MAF){
            if (!quiet) cat(paste(Sys.time(),"SNP under selection does NOT fulfil MAF criterion in sample from simulation",simID,"replicate",replic, "\n"))
            advantageous_allele_not_lost <-F
          }else{
            if (!quiet) cat(paste(Sys.time(),"SNP under selection fulfils MAF criterion in sample from simulation",simID,"replicate",replic, "\n"))

            
            
            pop1hap1 <- matrix(NA,nrow=sample_size,ncol=sample_size_loci)
            pop1hap2 <- matrix(NA,nrow=sample_size,ncol=sample_size_loci)
            pop2hap1 <- matrix(NA,nrow=sample_size,ncol=sample_size_loci)
            pop2hap2 <- matrix(NA,nrow=sample_size,ncol=sample_size_loci)
            for (i in 1:(2*sample_size) ){
              sampled_haplotypes_1_list[[i]]<-  setdiff(intersect(sampled_haplotypes_1_list[[i]],SNP_list$x),m2$x)
              sampled_haplotypes_2_list[[i]]<-  setdiff(union(intersect(sampled_haplotypes_2_list[[i]],SNP_list$x),SNP_list[which(SNP_list$n_pop.y==2*N),"x"]),m2$x)
            }
            loci_position_and_name <- matrix(1:sample_size_loci,nrow=sample_size_loci,ncol=1,dimnames=list(SNP_list$x,"loci"))
            for (i in 1:sample_size){
              haplo <- array(0,sample_size_loci)
              derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1_list[[(i*2)-1]],SNP_list$x)),]
              haplo[ derived_alleles ] <- 1
              pop1hap1[i,] <- haplo
              
              haplo <- array(0,sample_size_loci)
              derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1_list[[(i*2)]],SNP_list$x)),]
              haplo[ derived_alleles ] <- 1
              pop1hap2[i,] <- haplo
              
              haplo <- array(0,sample_size_loci)
              derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2_list[[(i*2)-1]],SNP_list$x)),]
              haplo[ derived_alleles ] <- 1
              pop2hap1[i,] <- haplo
              
              haplo <- array(0,sample_size_loci)
              derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2_list[[(i*2)]],SNP_list$x)),]
              haplo[ derived_alleles ] <- 1
              pop2hap2[i,] <- haplo
            }
            
            pop1Fis <- FIS.compute(H1=t(pop1hap1),
                                   H2=t(pop1hap2),
                                   sample_size=sample_size)
            pop2Fis <- FIS.compute(H1=t(pop2hap1),
                                   H2=t(pop2hap2),
                                   sample_size=sample_size)
            
            # MAKES INPUT FOR fastPHASE FOR HAPLOTYPE-BASED ANALYSES (1 FILE PER CHORMOSOME)

            #for (chr in 1:length(fastPHASE_in)){
            #  if (chr==1){
            #    limit1 <- 1
            #  }else{
            #    limit1 <- chromosome_structure[(chr*2)-2,1]
            #  }
            #  limit2 <- chromosome_structure[(chr*2)-1,1]
            #  
            #  loci <- setdiff(intersect(which(SNP_list$x>=limit1),which(SNP_list$x<=limit2)),m2$x)
            #  
            #  loci_position_and_name <- matrix(loci,nrow=length(loci),ncol=1,dimnames=list(SNP_list$x[loci],"loci"))
            #  
            #  write(sample_size*2,fastPHASE_in[chr],ncolumns=1)
            #  write(length(loci),fastPHASE_in[chr],ncolumns=1,append=T)
            #  write( c("P",SNP_list$x[loci]),fastPHASE_in[chr],ncolumns=length(loci)+1,append=T )
            #  
            #  for (i in 1:sample_size){
            #    write(paste0("# ind",i,"pop1"),fastPHASE_in[chr],append=T)
    #  
     #           haplo <- array(0,length(loci))
    #            derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1_list[[(i*2)-1]],SNP_list$x[loci])),]
     #           haplo[ derived_alleles ] <- 1
      #          write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
       #         pop1hap1[i,] <- haplo
                
             #   haplo <- array(0,length(loci))
            #    derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1_list[[(i*2)]],SNP_list$x[loci])),]
            #    haplo[ derived_alleles ] <- 1
            #    write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
            #    pop1hap2[i,] <- haplo
                
                
                
            #  }
            #  for (i in 1:sample_size){
            #    write(paste0("# ind",i,"pop2"),fastPHASE_in[chr],append=T)
                
            #    haplo <- array(0,length(loci))
            #    derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2_list[[(i*2)-1]],SNP_list$x[loci])),]
            #    haplo[ derived_alleles ] <- 1
            #    write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
            #    pop2hap1[i,] <- haplo
                
            #    haplo <- array(0,length(loci))
            #    derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2_list[[(i*2)]],SNP_list$x[loci])),]
            #    haplo[ derived_alleles ] <- 1
            #    write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
            #    pop2hap2[i,] <- haplo
            #  }
            #}
            #write(c(rep(1,sample_size),rep(2,sample_size)),ncolumns =sample_size*2,fastPHASE_pop)
            
             
            
            
            # ./fastPHASE -oTEST -utest_fastPHASE_pop.txt -i -Z -H200 -T1 -K8 -Pzp test_fastPHASE_in.txt
            
            
            
            
            # define sample of SNP loci that fulfil MAF criterion
            nx <- SNP_list[,"derived.x"] + SNP_list[,"ancestral.x"]
            px <- SNP_list[,"derived.x"] / nx
            ny <- SNP_list[,"derived.y"] + SNP_list[,"ancestral.y"]
            py <- SNP_list[,"derived.y"] / ny
            maf <- ( ((px + py) / 2 >= MAF_threshold) & ((px + py) / 2 <= 1-MAF_threshold) )
            
            if (!quiet) cat(paste(Sys.time(),nrow(SNP_list)-length(which(maf)),"of",nrow(SNP_list),"loci do not fulfil MAF criterion in sample for simulation",simID,"replicate",replic, "\n"))
          
            save(SNP_list,maf,file=mutable_name)
            
          
            # END ANALYSIS
            # analyse with Fc
            Fstat <- FstatFun_from_dataframe(SNP_list=SNP_list,
                                           maf=maf,
                                           dT=selection_period_duration,
                                           nbsimul=num_of_sim_Fc,
                                           MAF_threshold=MAF_threshold)
            
      
            #head(Fstat$Fc_list)
            
            #m2
            #head(SNP_list)
            
            chromosome_under_selection          <- ceiling( m2$x / (genome_length / chr_num) )
            first_bp_chromosome_under_selection <- (chromosome_under_selection-1) * (genome_length / chr_num)+1
            last_bp_chromosome_under_selection  <- (chromosome_under_selection) * (genome_length / chr_num)
            
            distance_bp <- array(NA,nrow(SNP_list))
          
            for (locus in intersect(which(SNP_list$x>=first_bp_chromosome_under_selection),which(SNP_list$x<=last_bp_chromosome_under_selection)) ){
              distance_bp[locus] <- abs(m2$x - SNP_list[locus,"x"])
            }
            distance_cM <- distance_bp * r * 100
            
            SNP_list <- cbind(SNP_list,distance_bp,distance_cM,Fstat$Fc_list,pop1Fis$Fis_all,pop2Fis$Fis_all)
            names(SNP_list) <- c(names(SNP_list)[1:19],"Fis_t0","Fis_t1")
            Fc <- Fstat$Fc_global
            Fis <- list(Fis_t0       <- pop1Fis$WC_Fis,
                        sigma_hat_t0 <- pop1Fis$sigma_hat,
                        Fis_t1       <- pop2Fis$WC_Fis,
                        sigma_hat_t1 <- pop2Fis$sigma_hat)
            
            save(SNP_list,maf,Fc,Fis,file=Fstat_name)
            save(pop1hap1,pop1hap2,pop2hap1,pop2hap2,file=haplotypes_sampled_file)
    
            files2save <- c(Fstat_name,
                            haplotypes_sampled_file,
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
            write(paste("Number of loci that do not fulfill MAF criterion:",nrow(SNP_list)-length(which(maf))), file=log_file,append=T)
            write(paste("Chormosome under selection:", chromosome_under_selection), file=log_file,append=T)
            
            
            
          }  
        }
      }
    }
  }
} # END for (replic in 1:num_of_replicates)

