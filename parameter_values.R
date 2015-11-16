# DEFAULT VALUES

# Set seed for random number generation
seed4random <- 123456
# Simulation ID for file identification
simID <- "test"
# selfing rate
sigma <- 0.0 
# mutation rate per bp
u <- 1e-8
# recombination rate per bp
r <- 1e-8
# genome size
genome_length <- 5e5
# number of chromosomes
chr_num <- 2
# effective population size (Ne)
N <- 500
# selection coeficient
sel_coef <- 0.0 
# dominance cofficient
dominance_coef <- 1
# length of pure drift period (number of times the population size)
number_of_times       <- 5 # see below
# Adaptation mode: "NM"=new mutation; "SV"=standing variation
selection_mode <- "SV"   
# number of generations between samples
selection_period_duration <- 20
# sample size
sample_size      <- 50   # number of diploid individuals sampled
sample_size_loci <- 5    # number of loci sampled for demographic inference)
#threshold for mimimum allele frequency
MAF_threshold   <- 0.01
# number of simulations for neutrality test 
num_of_sim_test <- 1000000
# replicate ID
replic <- 0
# if TRUE: estimate Ne from simulations only, do not perform neutrality tests
Ne_only <- T
# Do not output progress messages
quiet <- F
# Do not output whole population
no_whole_pop_out <- F


# VALUES FROM COMMAND LINE CALL ARGUMENTS

# gets parameter and setting values from command line (using package 'batch')
parseCommandArgs()

# set compound parameter values that may have change through command line call
#------------------------------------------------------------------------------

# set random seed
set.seed(seed4random)

if (Ne_only)             num_of_sim_Fc    <- 0

# length of pure drift period
drift_period_duration <- number_of_times*N

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
  }else if (advantageous_allele=="ancestral"){
    sel_coef <- -(sel_coef/(1+sel_coef))
    dominance_coef <- 1 - dominance_coef
  }
}

# FILE NAMING

# directory to file output
working_dir <- paste0(simID) 
system( paste("mkdir",working_dir) )
# name for log file for each replicate
log_file       <- paste0(working_dir,"/",simID,"_",replic,"_log.txt")
whole_pop_file <- paste0(working_dir,"/",simID,"_",replic,"_whole_pop.RData") 
data_file      <- paste0(working_dir,"/",simID,"_",replic,"_data.RData") 
F_stats_file   <- paste0(working_dir,"/",simID,"_",replic,"_F_stats.RData") 

