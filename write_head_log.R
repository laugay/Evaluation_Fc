
write("Evaluation of Fc statistics for detection of selection. Log file", file=log_file)
if(Ne_only){
  write("Estimation of Ne only. No test for selection performed.", file=log_file,append=T)
}

# Some population genetics quantities of interest:
theta         <- N*4*u*genome_length
Ns            <- N*sel_coef
expected_S    <- u*genome_length*4*N*sum(1/1:(sample_size[2]-1))
E_Ne          <- (2-sigma)*N  # number of gene copies
E_Fst         <- selection_period_duration/(selection_period_duration+2*E_Ne)

# check that sample sizes (individuals and loci) are OK with parameters
if (any(N<sample_size)) warning(paste(Sys.time(),"/!\\ Sample size larger than population size in scenario",simID,"\n"))
if (expected_S<sample_size_loci) warning(paste(Sys.time(),"/!\\ Number of loci to sample lower than the expected number of polymorphic sites in scenario",simID,"\n"))
# check that selection strength makes scenario be nearly-neutral
if (Ns>=-1 & Ns<=1) warning(paste(Sys.time(),"/!\\ Selection strength is very low: neutral or nearly neutral model in scenario",simID,"\n"))

write("", file=log_file,append=T)
write("DEMOGRAPHY", file=log_file,append=T)
write(paste("Population size:",N), file=log_file,append=T)
write(paste("Selfing rate:",sigma), file=log_file,append=T)
write(paste("Effective population size:",E_Ne), file=log_file,append=T)
write(paste("Expected Fst between time samples:",E_Fst), file=log_file,append=T)

write("", file=log_file,append=T)
write("MUTATION", file=log_file,append=T)
write(paste("Mutation rate per bp:",u), file=log_file,append=T)
write(paste("Genome length:",genome_length), file=log_file,append=T)
write(paste("Theta (genome wide):",theta), file=log_file,append=T)

write("", file=log_file,append=T)
write("RECOMBINATION", file=log_file,append=T)
write(paste("Recombination rate per bp:",r), file=log_file,append=T)
write(paste("Number of chormosomes:",chr_num), file=log_file,append=T)

write("", file=log_file,append=T)
write("SELECTION", file=log_file,append=T)
write(paste("Selection coefficient:",sel_coef), file=log_file,append=T)
write(paste("Dominance coefficient:",dominance_coef), file=log_file,append=T)
write(paste("Ns:",Ns), file=log_file,append=T)
if (sel_coef!=0){
  if(selection_mode=="NM"){
    write("Selection acting on a new mutation", file=log_file,append=T)
  }else{
    write("Selection acting on standing variation", file=log_file,append=T)  
    if (advantageous_allele=="derived"){
      write("Derived allele is advantageous", file=log_file,append=T)
    }else{
      write("Ancestral allele is advantageous", file=log_file,append=T)
    }
  }
}
write("", file=log_file,append=T)
write("SAMPLE", file=log_file,append=T)
write(paste("Pre-simulation for mutation-drift equilibrium (neutral):",drift_period_duration), file=log_file,append=T)
write(paste("Time period between samples (selection):",selection_period_duration), file=log_file,append=T)
write(paste("Sample size (individuals) time",1:2,":",sample_size), file=log_file,append=T)
write(paste("Sample size (loci):",sample_size_loci), file=log_file,append=T)
write(paste("Expected number of polymorphisms (in a one generation sample):",expected_S), file=log_file,append=T)
write(paste("Minimum allele frequency threshold:",MAF_threshold), file=log_file,append=T)
if(!Ne_only){
  write(paste("Simulations for hypothesis test:",num_of_sim_test), file=log_file,append=T)
}
write(paste("Seed for random number generation:",seed4random), file=log_file,append=T)

write("", file=log_file,append=T)
write("SIMULATION", file=log_file,append=T)
write(paste("Replicate:",replic), file=log_file,append=T)



