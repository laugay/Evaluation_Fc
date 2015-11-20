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

# functions to write/read SLiM input/output
source("slim_tools.R")
source("manipulate_output.R")
# functions to calculate F statistics and analyse data
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
source("sample_output.R")

# Analyse data
# 1. estimate Fst, Ne (all loci)
# 2. test for outliers (locus by locus)
results <- FC_outlier_test(SNP_data$genotype_data,
                           MAF_threshold=MAF_threshold,
                           delta_T=selection_period_duration,
                           num_of_sim_test=num_of_sim_test)

(results)




# Ne_hat_FC   <- EstimateNe.F_C (Fstats$F_C,selection_period_duration,sample_size,sample_size)
# selfing_hat <- 2*Fis_hat/(1+Fis_hat)
# write(paste("Estimated Fst between time samples:"       ,Fst_hat,    "; Expected value:", E_Fst), file=log_file,append=T)
# write(paste("Estimated effective population size (Fst):",Ne_hat_FST, "; True value:",     E_Ne) , file=log_file,append=T)
# write(paste("Estimated effective population size (Fc):" ,Ne_hat_FC,  "; True value:",     E_Ne) , file=log_file,append=T)







  
#


files2delete <- c(slim_out_drift, 
                  slim_out_selection, 
                  slim_init_selection,
                  slim_in_drift,
                  slim_in_selection,
                  slim_log_drift,
                  slim_log_selection)
for (item in files2delete){
  system( paste("rm",item) )
}





# test
# R --no-save --args sigma 0.0 sel_coef 0.0 simID "test" selection_mode "SV" selection_period_duration 20 seed4random 4444 replic 0 N 500 chr_num 2 genome_length 5e5 sample_size 50 sample_size_loci 5 < Main.R
# R --no-save --args sigma 0.0 sel_coef 0.0 simID "test" selection_mode "SV" selection_period_duration 20 seed4random 4444 replic 0 N 500 chr_num 2 genome_length 5e8 sample_size 50 sample_size_loci 1000 < Main.R
