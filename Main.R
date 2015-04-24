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
seed4random <- 3653698 

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
sigma <- 0.5 
# mutation rate per bp
u <- 1e-9
# recombination rate per bp
r <- 1e-8
# genome size
genome_length <- 2000000000
# number of chromosomes
chr_num <- 10
# effective population size (Ne)
N <- 500
# selection coeficient
sel_coef <- 0.5
# dominance cofficient
dominance_coef <- 1
# length of pure drift period (number of times the population size)
number_of_times       <- 10 # see below
# Adaptation mode: "NM"=new mutation; "SV"=standing variation
selection_mode <- "NM" #  
# number of generations between samples
selection_period_duration <- 25
# sample size
sample_size      <- 50      # number of individuals sampled
sample_size_loci <- 10000    # number of loci sampled for demographic inference)
#threshold for mimimum allele frequency
MAF_threshold <- 0.01
# number of simulations for testing Fc
num_of_sim_Fc <- 1000
# number of replicates for each scenario
num_of_replicates <- 10
starting_replicate <- 1
# number of attemps to get the selected allele present in tha last generation
# AND in the sample
# AND fulfilling MAF criterion
max_number_of_trials <- 100

#  Do not output progress messages
quiet <- T 

# gets parameter and setting values from command line (using package 'batch')
# example:
#   R --vanilla --args seed4random 104415 dominance_coef 0.5 < Main.R > Main.out
parseCommandArgs()

# Final setup lines
set.seed(seed4random)
#setwd(working_directory)

# Some population genetics quantities of interest:
theta         <- N*4*u*genome_length
Ns            <- N*sel_coef
expected_S    <- u*genome_length*4*N*sum(1/1:(sample_size-1))
Ne_corrected  <- N*(2-sigma)/2

# length of pure drift period
drift_period_duration <- number_of_times*N

# check that sample sizes (individuals and loci) are OK with parameters
if (N<sample_size) warning(paste(Sys.time(),"/!\\ Sample size larger than population size in simulation",simID,"\n"))
if (expected_S<sample_size_loci) warning(paste(Sys.time(),"/!\\ Number of loci to sample lower than the expected number of polymorphic sites in simulation",simID,"\n"))
# check that selection strength makes scenario not be nearly-neutral
if (Ns<=1) warning(paste(Sys.time(),"/!\\ Selection strength is too low, neutral or nearly neutral scenario in simulation",simID,"\n"))



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
    if (!quiet) cat(paste(Sys.time(),"Derived allele is advantageous in simulation",simID,"\n"))
  }else if (advantageous_allele=="ancestral"){
    sel_coef <- -(sel_coef/(1+sel_coef))
    dominance_coef <- 1 - dominance_coef
    if (!quiet) cat(paste(Sys.time(),"Ancestral allele is advantageous in simulation",simID,"\n"))
  }
}


#-----------
# SIMULATION
#-----------

system( paste("mkdir",simID) )
for (rep in starting_replicate:(starting_replicate+num_of_replicates)){

  # 1. SIMULATION OF A PURE DRIFT PERIOD
  
  # Set file names
  slim_in_drift  <- paste0(simID,"_",rep,"_drift.txt")
  slim_out_drift <- paste0(simID,"_",rep,"_drift.out")
  slim_log_drift <- paste0(simID,"_",rep,"_drift.log")
  
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
  
  if (!quiet) cat(paste(Sys.time(),"SIMULATION OF A PURE DRIFT PERIOD STARTS. Simulation:",simID,"replicate",rep,"\n"))
  
  # Run SLiM
  system(paste("./slim",slim_in_drift,">",slim_log_drift))
  # NB: SLiM executable must be in the same folder as SLiM input file (this is a requirement from SLiM)
  
  if (!quiet) cat(paste(Sys.time(),"END OF SIMULATION OF A PURE DRIFT PERIOD. Simulation:",simID,"replicate",rep,"\n"))

  advantageous_allele_not_lost <- FALSE
  counter <- 0
  while(!advantageous_allele_not_lost && counter<max_number_of_trials){
    counter <- counter+1
    
    # 2. SIMULATION OF THE PERIOD WITH SELECTION (in between samples)
    if (!quiet) cat(paste(Sys.time(),"Selection period for simulation:",simID,"replicate",rep,"attempt number",counter,"\n"))
    
    
    # Set file names
    slim_in_selection   <- paste0(simID,"_",rep,"_selection.txt")
    slim_out_selection  <- paste0(simID,"_",rep,"_selection.out")
    slim_log_selection  <- paste0(simID,"_",rep,"_selection.log")
    slim_init_selection <- paste0(simID,"_",rep,"_selection_init.txt")
    fastPHASE_in <-character()
    for (chr in 1:chr_num){
      fastPHASE_in <- c(fastPHASE_in,paste0(simID,"_",rep,"_chr_",chr,"_fastPHASE_in.txt"))
    }
    mutable_name        <- paste0(simID,"_",rep,"_muttable.RData") 
    Fstat_name          <- paste0(simID,"_",rep,"_Fstat.RData")
    sampled_ind_name    <- paste0(simID,"_",rep,"_sampled_ind.RData")
    
    
    
    
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
    writeOutput            (file=slim_in_selection, type="A", time=selection_period_duration, filename=slim_out_selection)
    writeOutput            (file=slim_in_selection, type="F", time=selection_period_duration, append_output=T)
    writeSeed              (file=slim_in_selection, seed=round(runif(1,-2^31,2^31)) )
    
    write("#INITIALIZATION",file=slim_in_selection,ncolumns=1,append=TRUE)
    write(slim_init_selection,file=slim_in_selection,ncolumns=1,append=TRUE)
    
    if (!quiet) cat(paste(Sys.time(),"SIMULATION OF THE PERIOD WITH SELECTION STARTS. Simulation:",simID,"replicate",rep,"\n"))
    system(paste("./slim",slim_in_selection,">",slim_log_selection))
    if (!quiet) cat(paste(Sys.time(),"END OF SIMULATION OF THE PERIOD WITH SELECTION. Simulation:",simID,"replicate",rep,"\n"))
    
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
      if (!quiet) cat(paste(Sys.time(),"Advantageous allele was lost in simulation",simID,"replicate",rep,": REPEATING SELECTION STEP\n"))  
    }else{
      if (!quiet) cat(paste(Sys.time(),"Advantageous allele was NOT lost in simulation",simID,"replicate",rep,"\n"))
      
      # Make a list of ALL polymorphic sites
      
      SNP_list <- Make.SNP.list(file_in_1=slim_init_selection,
                                file_in_2=slim_out_selection,
                                file_fixed=slim_log_selection,
                                populations=2,
                                sample_size=sample_size,
                                N=N)
      
      if (!quiet) cat(paste(Sys.time(),"List of all polymorphic SNP completed for simulation",simID,"replicate",rep,"\n"))     
      
      sampled_haplotypes_1_list <- SNP_list$sampled_haplotypes_1_list
      sampled_haplotypes_2_list <- SNP_list$sampled_haplotypes_2_list
      m2                        <- SNP_list$m2
      SNP_list                  <- SNP_list$SNP_list
      
      if (nrow(SNP_list)<sample_size_loci){
        if (!quiet) cat(paste(Sys.time(),"Less than",sample_size_loci,"polymorphic loci (",nrow(SNP_list),") in sample for simulation",simID,"replicate",rep,": REPEATING SELECTION STEP\n"))  
        advantageous_allele_not_lost <-F
      }else{
        if (!quiet) cat(paste(Sys.time(),"The number of polymorphic loci in the sample (",dim(SNP_list)[1],") allows to sample",
                  sample_size_loci,"loci for demographic inference in simulation",simID,"replicate",rep,"\n"))  
    
        # reduce table of polymoprhic SNP to sample_size_loci
        
        if (length(which(SNP_list[,"type"] == "m2")) == 0) {
          if (!quiet) cat(paste(Sys.time(),"SNP under selection is absente in sample from simulation",simID,"replicate",rep,": REPEATING SELECTION STEP\n"))
          # Sample all loci randomly
          advantageous_allele_not_lost <-F
        }else{
          if (!quiet) cat(paste(Sys.time(),"SNP under selection is present in sample from simulation",simID,"replicate",rep,"\n"))
          # Sample locus under selection and the rest randomly 
          SNP_list <- SNP_list[sort(c(sample(which(SNP_list[,"type"]=="m1"),sample_size_loci-1),which(SNP_list[,"type"]=="m2"))),]
          m2_MAF <- m2.gbmaf(SNP_list,MAF_threshold)      
          if (!m2_MAF){
            if (!quiet) cat(paste(Sys.time(),"SNP under selection does NOT fulfil MAF criterion in sample from simulation",simID,"replicate",rep,"\n"))
            advantageous_allele_not_lost <-F
          }else{
            if (!quiet) cat(paste(Sys.time(),"SNP under selection fulfils MAF criterion in sample from simulation",simID,"replicate",rep,"\n"))
 
            
            
            # MAKES INPUT FOR fastPHASE FOR HAPLOTYPE-BASED ANALYSES (1 FILE PER CHORMOSOME)
            for (i in 1:(2*sample_size) ){
              sampled_haplotypes_1_list[[i]]<-  setdiff(intersect(sampled_haplotypes_1_list[[i]],SNP_list$x),m2$x)
              sampled_haplotypes_2_list[[i]]<-  setdiff(union(intersect(sampled_haplotypes_2_list[[i]],SNP_list$x),SNP_list[which(SNP_list$n_pop.y==2*N),"x"]),m2$x)
            }
            for (chr in 1:length(fastPHASE_in)){
              if (chr==1){
                limit1 <- 1
              }else{
                limit1 <- chromosome_structure[(chr*2)-2,1]
              }
              limit2 <- chromosome_structure[(chr*2)-1,1]
              
              loci <- setdiff(intersect(which(SNP_list$x>=limit1),which(SNP_list$x<=limit2)),m2$x)
              
              loci_position_and_name <- matrix(loci,nrow=length(loci),ncol=1,dimnames=list(SNP_list$x[loci],"loci"))
              
              write(sample_size*2,fastPHASE_in[chr],ncolumns=1)
              write(length(loci),fastPHASE_in[chr],ncolumns=1,append=T)
              write( c("P",SNP_list$x[loci]),fastPHASE_in[chr],ncolumns=length(loci)+1,append=T )
              
              for (i in 1:sample_size){
                write(paste0("# ind",i,"pop1"),fastPHASE_in[chr],append=T)
      
                haplo <- array(0,length(loci))
                derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1_list[[(i*2)-1]],SNP_list$x[loci])),]
                haplo[ derived_alleles ] <- 1
                write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
                
                haplo <- array(0,length(loci))
                derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1_list[[(i*2)]],SNP_list$x[loci])),]
                haplo[ derived_alleles ] <- 1
                write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
              }
              for (i in 1:sample_size){
                write(paste0("# ind",i,"pop2"),fastPHASE_in[chr],append=T)
                
                haplo <- array(0,length(loci))
                derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2_list[[(i*2)-1]],SNP_list$x[loci])),]
                haplo[ derived_alleles ] <- 1
                write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
                
                haplo <- array(0,length(loci))
                derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2_list[[(i*2)]],SNP_list$x[loci])),]
                haplo[ derived_alleles ] <- 1
                write(paste(haplo,collapse=""),fastPHASE_in[chr],append=T)
              }
            }
            
            
            
            
            
            
            
            # define sample of SNP loci that fulfil MAF criterion
            nx <- SNP_list[,"derived.x"] + SNP_list[,"ancestral.x"]
            px <- SNP_list[,"derived.x"] / nx
            ny <- SNP_list[,"derived.y"] + SNP_list[,"ancestral.y"]
            py <- SNP_list[,"derived.y"] / ny
            maf <- ( ((px + py) / 2 >= MAF_threshold) & ((px + py) / 2 <= 1-MAF_threshold) )
            
            if (!quiet) cat(paste(Sys.time(),nrow(SNP_list)-length(which(maf)),"of",nrow(SNP_list),"loci do not fulfil MAF criterion in sample for simulation",simID,"replicate",rep,"\n"))
          
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
            
            SNP_list <- cbind(SNP_list,distance_bp,distance_cM,Fstat$Fc_list)
            Fc <- Fstat$Fc_global
            
            save(SNP_list,maf,Fc,file=Fstat_name)
        
            files2save <- c(Fstat_name,
                            fastPHASE_in,
                            slim_in_drift,
                            slim_init_selection,
                            slim_in_selection)
            
            files2delete <- c(slim_out_drift, 
                              slim_log_drift, 
                              slim_out_selection, 
                              slim_log_selection,  
                              mutable_name)
            
            
            for (item in files2save){
              system( paste("mv",item,simID) )
            }
            for (item in files2delete){
              system( paste("rm",item) )
            }
          }  
        }
      }
    }
  }
} # END for (rep in 1:num_of_replicates)