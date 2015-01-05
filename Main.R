##################################################################################
#
# Script to test the Fc outlier method modified from Goldringer and Bataillon 2002
# 
##################################################################################

# This scripts uses SLiM (Messer 2013) for simulating the data
# http://messerlab.org/software/

# This script uses the following packages:
require(adegenet)
require(pegas)
require(seqinr)


#---------------
# MAIN SETTINGS
#---------------

# Set seed for random number generation
set.seed(368898)

# set working directory
setwd("/home/miguel/Work/Research/2012.SelfAdapt/Evaluation_Fc")

# high penality for scientific notation (necessary to input large number in SLiM)
options("scipen"=999)

# tools to write/read SLiM input/output
source("slim_tools.R")


#------------
# PARAMETERS
#------------

# selfing rate
sigma <- 0.9 
# mutation rate per bp
u <- 2e-6
# genome size
genome_length <- 25000
# theta= 4Neu (for the genome)
theta <- 50
# population size
N <- theta/(4*u*genome_length)
# selection coeficient
new_scoef <- 0.6
# length of pure drift period
drift_period_duration <- 30*N


#-----------
# SIMULATION
#-----------

# 1. SIMULATION OF A PURE DRIFT PERIOD

######################## SIMULATION OF A PURE DRIFT PERIOD ########################################################

simID    <- "drift"
slim_in  <- paste(simID,".txt",sep="")
slim_out <- paste(simID,".out",sep="")
slim_log <- paste(simID,".log",sep="")

writeMutation(file=slim_in, number_of_types=1, h=0.5, DFE="f", s=0)
writeMutationRate(file=slim_in, u=u)  
writeGenomicElement(file=slim_in, number_of_types=1, mut_type=list("m1"), prop=list(1))
writeChromosome(file=slim_in, element_type="g1", start=1, end=genome_length)
writeRecombination(file=slim_in, interval_end=genome_length, r=0)
writeGenerations(file=slim_in, t=drift_period_duration, append=T)
writeDemography(file=slim_in, type="P", time=1, pop="p1", N=N)
writeDemography(file=slim_in, type="S", time=1, pop="p1", sigma=sigma, append_demography=T)
writeOutput(file=slim_in, type="A", time=drift_period_duration, filename=slim_out)
writeSeed(file=slim_in, seed=round(runif(1,-2^31,2^31)) )

system(paste("./slim",slim_in,">",slim_log))
# NB: slim executable must be in the same folder as slim input file






