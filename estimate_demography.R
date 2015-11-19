
# compute Fst, Fis, Fc
Fstats   <- compute_Fstats(SNP_data$genotype_data)
Fst_hat  <- Fstats$F_ST
Fis_hat  <- Fstats$F_IS
FC       <- Fstats$F_C
#remove(Fstats)

# compute Ne estimates
Ne_hat_FC   <- EstimateNe.F_C (Fstats$F_C,selection_period_duration,sample_size,sample_size)
Ne_hat_FST  <- EstimateNe.F_ST(Fstats$F_ST,selection_period_duration)
selfing_hat <- 2*Fis_hat/(1+Fis_hat)

write(paste("Estimated Fst between time samples:"       ,Fst_hat,    "; Expected value:", E_Fst), file=log_file,append=T)
write(paste("Estimated effective population size (Fst):",Ne_hat_FST, "; True value:",     E_Ne) , file=log_file,append=T)
write(paste("Estimated effective population size (Fc):" ,Ne_hat_FC,  "; True value:",     E_Ne) , file=log_file,append=T)
#write(paste("Estimated Fst between time samples:" ,FC,    "; Expected value:", E_FC), file=log_file,append=T)


save(Fstats,Ne_hat,selfing_hat,file=F_stats_file)