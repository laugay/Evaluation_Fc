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

# high penality for scientific notation (necessary to input large number in SLiM)
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
# genome size
genome_length <- 25000
# theta= 4Neu (for the genome)
theta <- 50
# selection coeficient
sel_coef <- 0.6
# dominance cofficient
dominance_coef <- 1
# length of pure drift period
number_of_times       <- 30 # see below
# Adaptation mode: "standing variation" or "new mutation" ()
mode <- "standing variation" 
# number of generations between samples
selection_period_duration <- 20

# gets parameter and setting values from command line
# example:
#   R --vanilla --args seed4random 104415 dominance_coef 0.5 < Main.R > Main.out
parseCommandArgs()

# population size
N <- theta/(4*u*genome_length)
# length of pure drift period
drift_period_duration <- number_of_times*N




#-----------
# SIMULATION
#-----------

# 1. SIMULATION OF A PURE DRIFT PERIOD

# Set file names
slim_in_drift  <- paste0(simID,"_drift.txt")
slim_out_drift <- paste0(simID,"_drift.out")
slim_log_drift <- paste0(simID,"_drift.log")

# Write slim input file
writeMutation      (file=slim_in_drift, number_of_types=1, h=0.5, DFE="f", s=0)
writeMutationRate  (file=slim_in_drift, u=u)  
writeGenomicElement(file=slim_in_drift, number_of_types=1, mut_type=list("m1"), prop=list(1))
writeChromosome    (file=slim_in_drift, element_type="g1", start=1, end=genome_length)
writeRecombination (file=slim_in_drift, interval_end=genome_length, r=0)
writeGenerations   (file=slim_in_drift, t=drift_period_duration, append=T)
writeDemography    (file=slim_in_drift, type="P", time=1, pop="p1", N=N)
writeDemography    (file=slim_in_drift, type="S", time=1, pop="p1", sigma=sigma, append_demography=T)
writeOutput        (file=slim_in_drift, type="A", time=drift_period_duration, filename=slim_out_drift)
writeSeed          (file=slim_in_drift, seed=round(runif(1,-2^31,2^31)) )

# Run slim
system(paste("./slim",slim_in_drift,">",slim_log_drift))
# NB: slim executable must be in the same folder as slim input file

# 2. SIMULATION OF THE PERIOD WITH SELECTION (in between samples)

# Set file names
slim_in_selection   <- paste0(simID,"_selection.txt")
slim_out_selection  <- paste0(simID,"_selection.out")
slim_log_selection  <- paste0(simID,"_selection.log")
slim_init_selection <- paste0(simID,"_selection_init.txt")

# Read slim output from DRIFT period
out_drift_lines <- readLines(con=slim_out_drift)
pop_line        <- which(out_drift_lines=="Populations:")
mut_line        <- which(out_drift_lines=="Mutations:")
gen_line        <- which(out_drift_lines=="Genomes:")
num_of_mut      <- gen_line-mut_line-1

if (mode=="new mutation")       cat("Adaptation though new mutation not implemented yet, using standing variation instead")
if (mode!="standing variation") cat("Adaptation mode undefined, using standing variation") 

# SELECTION ON STANDING VARIATION (changes selection coefficient of on a random locus)

# chose a random locus
sampled_mut                   <- sample(x=num_of_mut, size=1)
# change selection coefficient of locus (derived allele being advantageous) # TODO: choose randomly the advanageuous allele
mutation_table                <- read.table(file=slim_out_drift, skip=mut_line, nrows=num_of_mut)
levels(mutation_table[,2])    <- c("m1","m2") 
mutation_table[sampled_mut,2] <- "m2"
mutation_table[sampled_mut,4] <- sel_coef
mutation_table[sampled_mut,5] <- dominance_coef

# write initialzation file for slim (state of populaion at starting point of selection period)
write(out_drift_lines[pop_line:mut_line], slim_init_selection)
write.table(mutation_table, slim_init_selection, append=T, quote=F, col.names=F)
write(out_drift_lines[-(1:(gen_line-1))], slim_init_selection, append=T)

# write input file for slim
writeMutation      (file=slim_in_selection, number_of_types=2, h=0.5, DFE="f", s=c(0,sel_coef), append=F, append_mutation=F)
writeMutationRate  (file=slim_in_selection, u=u)  
writeGenomicElement(file=slim_in_selection, number_of_types=1, mut_type=list("m1","m2"), prop=list(1,0))
writeChromosome    (file=slim_in_selection, element_type="g1", start=1, end=genome_length)
writeRecombination (file=slim_in_selection, interval_end=genome_length, r=0)
writeGenerations   (file=slim_in_selection, t=selection_period_duration, append=T)
writeOutput        (file=slim_in_selection, type="A", time=selection_period_duration, filename=slim_out_selection)
writeSeed          (file=slim_in_selection, seed=round(runif(1,-2^31,2^31)) )

write("#INITIALIZATION",file=slim_in_selection,ncolumns=1,append=TRUE)
write(slim_init_selection,file=slim_in_selection,ncolumns=1,append=TRUE)

system(paste("./slim",slim_in_selection,">",slim_log_selection))


