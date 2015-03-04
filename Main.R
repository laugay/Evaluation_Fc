##################################################################################
#
# Script to test the Fc outlier method modified from Goldringer and Bataillon 2002
# 
##################################################################################

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
seed4random <- 368898 
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
chr_num <- 8
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

# gets parameter and setting values from command line (using package 'batch')
# example:
#   R --vanilla --args seed4random 104415 dominance_coef 0.5 < Main.R > Main.out
parseCommandArgs()

# population size
N <- theta/(4*u*genome_length)
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
advantageus_allele <- sample(c("derived","ancestral"),size=1)
if (advantageus_allele=="derived"){
  sel_coef <- sel_coef
  dominance_coef <- dominance_coef
}else if (advantageus_allele=="derived"){
  sel_coef <- -(sel_coef/(1+sel_coef))
  dominance_coef <- 1 - dominance_coef
}
# change selection coefficient of locus
mutation_table                <- read.table(file=slim_out_drift, skip=out_drift_mut_line, nrows=out_drift_num_of_mut)
levels(mutation_table[,2])    <- c("m1","m2") 
mutation_table[sampled_mut,2] <- "m2"
mutation_table[sampled_mut,4] <- sel_coef
mutation_table[sampled_mut,5] <- dominance_coef

# write initialzation file for slim (state of population at starting point of selection period)
write(out_drift_lines[out_drift_pop_line:out_drift_mut_line], slim_init_selection)
write.table(mutation_table, slim_init_selection, append=T, quote=F, col.names=F)
write(out_drift_lines[-(1:(out_drift_gen_line-1))], slim_init_selection, append=T)

# message(paste(Sys.time(),"Episode de SV slim commence, répétition n° ",sim))

# write input file for slim
writeMutation          (file=slim_in_selection, number_of_types=2, h=0.5, DFE="f", s=c(0,sel_coef), append=F, append_mutation=F)
writeMutationRate      (file=slim_in_selection, u=u)  
writeGenomicElement    (file=slim_in_selection, number_of_types=1, mut_type=list("m1","m2"), prop=list(1,0))
writeChromosome        (file=slim_in_selection, element_type="g1", start=1, end=genome_length)
writeRecombinationChrom(file=slim_in_selection, chr_num=chr_num, genome_length=genome_length, r=r, append=T)  
writeGenerations       (file=slim_in_selection, t=selection_period_duration, append=T)
writeOutput            (file=slim_in_selection, type="A", time=selection_period_duration, filename=slim_out_selection)
writeSeed              (file=slim_in_selection, seed=round(runif(1,-2^31,2^31)) )

write("#INITIALIZATION",file=slim_in_selection,ncolumns=1,append=TRUE)
write(slim_init_selection,file=slim_in_selection,ncolumns=1,append=TRUE)

system(paste("./slim",slim_in_selection,">",slim_log_selection))

# message(paste(Sys.time(),"Episode de SV slim terminé, répétition n° ",sim))

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

if (log_selection_num_mut > 0) {
  if (length(grep(pattern="m2",x=out_selection_lines[(out_selection_mut_line+1):(out_selection_gen_line-1)]))>0 || length(grep(pattern="m2",x=log_selection_lines[(log_selection_mut_line+1):length(log_selection_lines)]))>0){
    m2succeed <- TRUE
  }else{ m2succeed <- FALSE }  
}else if(log_selection_num_mut==0){
  if (length(grep(pattern="m2",x=out_selection_lines[(out_selection_mut_line+1):(out_selection_gen_line-1)]))>0){
    m2succeed <- TRUE
  }else{ m2succeed <- FALSE }
}














