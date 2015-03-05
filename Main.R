###################################################################################
#
# Script to test the Fc outlier method modified from Goldringer and Bataillon 2002
#
# by A. Becheler, R. Vitalis & M. Navascués
# 
###################################################################################

# This scripts uses SLiM (Messer 2013) for simulating the data
# http://messerlab.org/software/

# This script uses the following packages:
require(batch)
require(adegenet)
require(pegas)
require(seqinr)


#---------------
# MAIN SETTINGS
#---------------

# Set seed for random number generation
seed4random <- 3698 
set.seed(seed4random)

# set working directory
working_directory <- "/home/miguel/Work/Research/2012.SelfAdapt/Evaluation_Fc"
setwd(working_directory)

# high penality for scientific notation (necessary to input large numbers in SLiM)
options("scipen"=999)

# tools to write/read SLiM input/output
source("slim_tools.R")
source("slim_tools_2.R")
source("slim_tools_RV_2014-05-06.R")

# Simulation ID for file identification
simID <- "test"


#----------------------------
# PARAMETERS: DEFAULT VALUES
#----------------------------

# selfing rate
sigma <- 0.9 
# mutation rate per bp
u <- 2e-6
# recombination rate per bp
r <- 2e-6
# genome size
genome_length <- 25000
# number of chromosomes
chr_num <- 10
# theta= 4Neu (for the genome)
theta <- 50
# selection coeficient
sel_coef <- 0.6
# dominance cofficient
dominance_coef <- 1
# length of pure drift period (number of times the population size)
number_of_times       <- 30 # see below
# Adaptation mode: "standing variation" or "new mutation" ()
mode <- "standing variation" 
# number of generations between samples
selection_period_duration <- 20
# sample size
sample_size <- 50

# gets parameter and setting values from command line (using package 'batch')
# example:
#   R --vanilla --args seed4random 104415 dominance_coef 0.5 < Main.R > Main.out
parseCommandArgs()

# population size
N <- theta/(4*u*genome_length)
if (N<sample_size) message(paste(Sys.time(),"/!\\ Sample size larger than population size in simulation",simID))
# length of pure drift period
drift_period_duration <- number_of_times*N
# Simulation scenario ID for file identification
simID <- paste0("scenario_",simID)



#-----------
# SIMULATION
#-----------

# 1. SIMULATION OF A PURE DRIFT PERIOD

# Set file names
slim_in_drift  <- paste0(simID,"_drift.txt")
slim_out_drift <- paste0(simID,"_drift.out")
slim_log_drift <- paste0(simID,"_drift.log")

# Write slim input file (functions from slim_tools_2.R)
writeMutation          (file=slim_in_drift, number_of_types=1, h=0.5, DFE="f", s=0)
writeMutationRate      (file=slim_in_drift, u=u)  
writeGenomicElement    (file=slim_in_drift, number_of_types=1, mut_type=list("m1"), prop=list(1))
writeChromosome        (file=slim_in_drift, element_type="g1", start=1, end=genome_length)
writeRecombinationChrom(file=slim_in_drift, chr_num=chr_num, genome_length=genome_length, r=r, append=T)
writeGenerations       (file=slim_in_drift, t=drift_period_duration, append=T)
writeDemography        (file=slim_in_drift, type="P", time=1, pop="p1", N=N)
writeDemography        (file=slim_in_drift, type="S", time=1, pop="p1", sigma=sigma, append_demography=T)
writeOutput            (file=slim_in_drift, type="A", time=drift_period_duration, filename=slim_out_drift)
writeSeed              (file=slim_in_drift, seed=round(runif(1,-2^31,2^31)) )

message(paste(Sys.time(),"SIMULATION OF A PURE DRIFT PERIOD STARTS. Simulation:",simID))

# Run SLiM
system(paste("./slim",slim_in_drift,">",slim_log_drift))
# NB: SLiM executable must be in the same folder as SLiM input file (this is a requirement from SLiM)

message(paste(Sys.time(),"END OF SIMULATION OF A PURE DRIFT PERIOD. Simulation:",simID))


# 2. SIMULATION OF THE PERIOD WITH SELECTION (in between samples)

# Set file names
slim_in_selection   <- paste0(simID,"_selection.txt")
slim_out_selection  <- paste0(simID,"_selection.out")
slim_log_selection  <- paste0(simID,"_selection.log")
slim_init_selection <- paste0(simID,"_selection_init.txt")
selestim_in         <- paste0(simID,"_selestim_in")
selestim_out        <- paste0(simID,"_selestim_out")
mutable_name        <- paste0(simID,"_muttable.RData") 
Fstat_name          <- paste0(simID,"_Fstat.RData")
sampled_ind_name    <- paste0(simID,"_sampled_ind.RData")




# Read slim output from DRIFT period
out_drift_lines      <- readLines(con=slim_out_drift)
out_drift_pop_line   <- which(out_drift_lines=="Populations:")
out_drift_mut_line   <- which(out_drift_lines=="Mutations:")
out_drift_gen_line   <- which(out_drift_lines=="Genomes:")
out_drift_num_of_mut <- out_drift_gen_line-out_drift_mut_line-1

if (mode=="new mutation")       cat("Adaptation though new mutation not implemented yet, using standing variation instead")
if (mode!="standing variation") cat("Adaptation mode undefined, using standing variation") 
mode <- "standing variation" # TO DO: implement new mutation scenario

# SELECTION ON STANDING VARIATION (changes selection coefficient of on a random locus)

# chose a random locus
sampled_mut                   <- sample(x=out_drift_num_of_mut, size=1)
# choose advantageous allele between derived and ancestral state 
advantageous_allele <- sample(c("derived","ancestral"),size=1)
if (advantageous_allele=="derived"){
  sel_coef <- sel_coef
  dominance_coef <- dominance_coef
  message(paste(Sys.time(),"Derived allele is advantageous in simulation",simID))
}else if (advantageous_allele=="ancestral"){
  sel_coef <- -(sel_coef/(1+sel_coef))
  dominance_coef <- 1 - dominance_coef
  message(paste(Sys.time(),"Ancestral allele is advantageous in simulation",simID))
}
# change selection coefficient of locus
mutation_table                <- read.table(file=slim_out_drift, skip=out_drift_mut_line, nrows=out_drift_num_of_mut)
levels(mutation_table[,2])    <- c("m1","m2") 
mutation_table[sampled_mut,2] <- "m2"
mutation_table[sampled_mut,4] <- sel_coef
mutation_table[sampled_mut,5] <- dominance_coef

# write initialzation file for slim (state of population at starting point of selection period)
write(out_drift_lines[out_drift_pop_line:out_drift_mut_line], slim_init_selection)
write.table(mutation_table, slim_init_selection, append=T, quote=F, col.names=F, row.names=F)
write(out_drift_lines[-(1:(out_drift_gen_line-1))], slim_init_selection, append=T)

# message(paste(Sys.time(),"Episode de SV slim commence, répétition n° ",sim))

# write input file for slim
writeMutation          (file=slim_in_selection, number_of_types=2, h=c(0.5,dominance_coef), DFE="f", s=c(0,sel_coef), append=F, append_mutation=F)
writeMutationRate      (file=slim_in_selection, u=u)  
writeGenomicElement    (file=slim_in_selection, number_of_types=1, mut_type=list(c("m1","m2")), prop=list(c(1,0)))
writeChromosome        (file=slim_in_selection, element_type="g1", start=1, end=genome_length)
writeRecombinationChrom(file=slim_in_selection, chr_num=chr_num, genome_length=genome_length, r=r, append=T)  
writeGenerations       (file=slim_in_selection, t=selection_period_duration, append=T)
writeOutput            (file=slim_in_selection, type="A", time=selection_period_duration, filename=slim_out_selection)
writeSeed              (file=slim_in_selection, seed=round(runif(1,-2^31,2^31)) )

write("#INITIALIZATION",file=slim_in_selection,ncolumns=1,append=TRUE)
write(slim_init_selection,file=slim_in_selection,ncolumns=1,append=TRUE)

message(paste(Sys.time(),"SIMULATION OF THE PERIOD WITH SELECTION STARTS. Simulation:",simID))
system(paste("./slim",slim_in_selection,">",slim_log_selection))
message(paste(Sys.time(),"END OF SIMULATION OF THE PERIOD WITH SELECTION. Simulation:",simID))

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

if (advantageous_allele_not_lost) {
  message(paste(Sys.time(),"Advantageous allele was NOT lost in simulation",simID))
}else{
  message(paste(Sys.time(),"Advantageous allele was lost in simulation",simID))  
}

# Make a list of plymorphic sites

SNP_list <- Make.SNP.list(file_in_1=slim_init_selection,
                          file_in_2=slim_out_selection,
                          file_fixed=slim_log_selection,
                          populations=2, sample_size=sample_size, samp_ind_name=sampled_ind_name)















######################### CONTINUE OR NOT THE ANALISIS OF THE SIMULATION #######################################

if(! m2succeed){
  # STOP ANALYSIS
  sim <- sim
  message(paste(Sys.time(),"m2 était absente de la simul slim, répétition n° ",sim))
  
}else if(m2succeed){
  message(paste(Sys.time(),"m2 était présente dans la simul slim, répétition n° ",sim))
  
  # TRY SAMPLINGS TO BEST EXPLOIT THE SIMULATION
  samp_try <- 1
  m2_MAF <- FALSE
  
  while(samp_try < 5){
    print(paste("sample try n°",samp_try))
    # individuals are sampled and SNP_list is made
    SNP_list <- Make.SNP.list(file_in_1=slim_init_selection,
                                 file_in_2=slim_out_selection,
                                 file_fixed=slim_log_selection,
                                 populations=2, sample_size=50, samp_ind_name=sampled_ind_name)
    message(paste(Sys.time(),"On vient d'achever la fonction Make.SNP, répétition n° ",sim))     
    if (length(which(SNP_list[,"type"] == "m2")) == 0) {
      message(paste(Sys.time(),"m2 absente..., essai ",samp_try, "répétition n° ",sim))        
      m2_MAF <- FALSE
      samp_try <- samp_try +1
    } else {
      m2_MAF <- m2.gbmaf(SNP_list)      
      if (m2_MAF == FALSE){
        # RE-TRY SAMPLINGS
        samp_try <- samp_try +1
        message(paste(Sys.time(),"Retente un échantillonnage, essai ",samp_try, "répétition n° ",sim))        
      } else if(m2_MAF == TRUE) {
        # STOP SAMPLINGS
        samp_try <- 5
        print("m2 sampled ! ")
        message(paste(Sys.time(),"On a réussi à échantillonner m2, répétition n° ",sim))
      }
    }
  } # end of while(samp_try < 5)
  
  if(m2_MAF==FALSE){
    print("m2 did not fil MAF recquirements")
    message(paste(Sys.time(),"m2 ne remplissait pas la MAF -> autre simul slim, répétition n° ",sim))
    
    # STOP ANALYSIS
    sim <- sim
  }else if(m2_MAF==TRUE){
    print("m2 has been sampled AND fills MAF recquirements")
    message(paste(Sys.time(),"m2 remplissait la MAF, on va échantillonner les SNP, répétition n° ",sim))
    
    # ensure that m2 will be in the analisis... and sample 999 other loci
    
    m2 <- which(SNP_list[,"type"]=="m2")
    SNP_list_m2 <- SNP_list[m2,]
    SNP_list_without_m2 <- SNP_list[-m2,]
    
    nx <- SNP_list_without_m2[,"derived.x"] + SNP_list_without_m2[,"ancestral.x"]
    px <- SNP_list_without_m2[,"derived.x"] / nx
    ny <- SNP_list_without_m2[,"derived.y"] + SNP_list_without_m2[,"ancestral.y"]
    py <- SNP_list_without_m2[,"derived.y"] / ny
    maf <- ((px + py) / 2 >= 0.01 & (px + py) / 2 <= 0.99)
    #		sample <- sample(nrow(SNP_list_without_m2[maf,]),size = 999)
    sample <- sample(which(maf),size = 999)
    SNP_list_m1 <- SNP_list_without_m2[sample,]
    SNP_list_m1_m2 <- rbind(SNP_list_m2,SNP_list_m1)
    SNP_final_list <- SNP_list_m1_m2[order(SNP_list_m1_m2 $x,decreasing=FALSE),]
    
    #save informations
    save(SNP_final_list,file=mutable_name)
    SNP_final_list <- SNP_final_list[,c("derived.x","ancestral.x","derived.y","ancestral.y")]
    
    # write selestim_file
    if (nrow(SNP_final_list)>0 ) {
      con <- file(selestim_in, open = "w")
      vect <- c(as.character(2),as.character(nrow(SNP_final_list)))
      writeLines(vect, con = con)
      write.table(x=SNP_final_list,file=con,row.names=FALSE,col.names=FALSE)
      close(con)
    }else{print("Your SNP list is empty : no polymorphism detected ")}
    
    message(paste(Sys.time(),"Fichier selestim écrit, Analyse Fstat commence, répétition n° ",sim))
    
    # END ANALYSIS
    # analyse with Fstatistics
    Fstat <- FstatFun_RV(data_file=selestim_in,dT=dT,nbsimul=simF)
    
    save(Fstat,file=Fstat_name)
    message(paste(Sys.time(),"Analyse Fstat terminée, simulation terminée, répétition n° ",sim))
    sim<- sim + 1
  } # end of else if(m2_MAF==TRUE)
} # end of else if(m2succeed){












