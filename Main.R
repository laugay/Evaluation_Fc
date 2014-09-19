##################################################################################
#
# Script to test the Fc outlier method modified from Goldringer and Bataillon 2002
# 
##################################################################################

# This scripts uses SLiM (Messer 2013) for simulating the data
# http://messerlab.org/software/

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
N <- theta/(4*u*genome_length)
# selection coeficient
new_scoef <- 0.6

set.seed(368898)

#

