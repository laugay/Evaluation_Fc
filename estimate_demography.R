
# compute Fst, Fis, Fc
Fstats   <- compute_Fstats(SNP_data$genotype_data)
Fst_hat  <- Fstats$F_ST
Fis_hat  <- Fstats$F_IS

# compute Ne estimates
if (Fst_hat>0){
  Ne_hat <- selection_period_duration * (1 - Fst_hat) / (2 * Fst_hat) 
}else{
  Ne_hat <- NA
}
selfing_hat <- 2*Fis_hat/(1+Fis_hat)

write(paste("Estimated effective population size:",Ne_hat,  "; True value:",     E_Ne) , file=log_file,append=T)
write(paste("Estimated Fst between time samples:" ,Fst_hat, "; Expected value:", E_Fst), file=log_file,append=T)


save(Fstats,Ne_hat,selfing_hat,file=F_stats_file)