#	The input data file should be formatted as follows:
#
#	1	1	0	1	2	2	1	0
#	2	1	0	1	0	0	1	2
#	1	2	1	1	0	0	1	2
#	2	2	2	1	1	2	0	2
#	3	2	2	2	1	1	1	1
#
#	Each line corresponds to an individual
# 	The first column contains IDs for individual
# 	The second column contains IDs for sampled demes
# 	The next columns correspond to loci
# 	`0' corresponds to homozygotes for type-1 alleles
# 	`1' corresponds to heterozygotes
# 	`2' corresponds to homozygotes for type-2 alleles

# data as an R matrix with appropriate format or infile as a text file in the appropriate format
# MAF_threshold:   threshold to filter loci with a mi9numim maf (minor allele frequency)
# delta_T:         number of generations between time samples
# num_of_sim_test: number of replicates of drift simulations to build null FC distribution
FC_outlier_test <- function(data,infile=NA,MAF_threshold,delta_T,num_of_sim_test){
  # read the data
  if (!is.na(infile)) data <- read.table(infile) 
  
  nbr.pops <- length(unique(data[,2]))																# compute the number of samples
  nbr.loci <- ncol(data)-2																						# compute the number of loci
  stopifnot(nbr.pops==2) 
  
  
  # compute Fst, Fis, Fc
  Fstats   <- compute_Fstats(data,MAF_threshold=MAF_threshold)
  # compute Ne estimates
  Ne_hat_FST  <- EstimateNe.F_ST(Fstats$F_ST,selection_period_duration)
  Ne_hat      <- round(Ne_hat_FST)
  
  # set of unique initial allele frequencies (minor allele)
  starting_freq <- Fstats$p[1,] 
  starting_freq[starting_freq>0.5] <- 1-starting_freq[starting_freq>0.5]
  

  parameter_combinations <- matrix(nrow=0,ncol=3)
  p_value <- array(NA,dim=nbr.loci)
  for (locus in seq_len(nbr.loci)) {
    
    if (Fstats$maf[locus]){
      position <- intersect(intersect(which(parameter_combinations[,1]==starting_freq[locus]),
                                      which(parameter_combinations[,2]==Fstats$n[1,locus])),
                                      which(parameter_combinations[,3]==Fstats$n[2,locus]))
      if (length(position)==1){
        FC_distribution <- get(paste0("FC_distribution_",position))
      }else{
        parameter_combinations <- rbind(parameter_combinations,c(starting_freq[locus],Fstats$n[1,locus],Fstats$n[2,locus]))
        FC_distribution <- Drift.simulation.4.test(Ne=Ne_hat,
                                                   dT=delta_T,
                                                   starting_freq=starting_freq[locus],
                                                   sample_size=Fstats$n[,locus],
                                                   num_of_sim_test=num_of_sim_test,
                                                   MAF_threshold=MAF_threshold)
        assign( paste0("FC_distribution_",nrow(parameter_combinations)) , FC_distribution )
      }
      p_value[locus] <- sum(FC_distribution[as.numeric(names(FC_distribution)) >= Fstats$F_C_locus[locus] ])/num_of_sim_test
      remove(FC_distribution)
    }
  }
  
  results_total    <- list(F_ST  =Fstats$F_ST,
                           F_IS  =Fstats$F_IS,
                           F_C   =Fstats$F_C,
                           F_IS  =Fstats$F_IS,
                           Ne_hat=Ne_hat_FST)
  
  results_by_locus <- data.frame(F_ST=Fstats$F_ST_locus,
                            F_IS=Fstats$F_IS_locus,
                            F_C =Fstats$F_C_locus,
                            p_value,
                            maf=Fstats$maf)

  return(list(results_total=results_total,results_by_locus=results_by_locus)) 
}
























Drift.simulation.4.test <- function(Ne,dT,starting_freq,sample_size,num_of_sim_test,MAF_threshold){
  if (starting_freq==0){
    px <- rep(1/(2*sample_size[1]) ,num_of_sim_test)
  }else{
    px <- rep(starting_freq,num_of_sim_test)
  }
  starting_freq <- rep(starting_freq,num_of_sim_test)
  py <- rep(NA,num_of_sim_test)
  sim <- 1
  while (sim <= num_of_sim_test) {
    f.drift <- px[sim]
    for (g in 1:dT) {
      f.drift <- (rbinom(1,Ne,f.drift)) / Ne
    }
    py[sim] <- rbinom(1, 2*sample_size[2],f.drift) / (2 * sample_size[2])
    if (((px[sim] + py[sim]) / 2) >= MAF_threshold && ((px[sim] + py[sim]) / 2) <= (1-MAF_threshold) ) {
      sim <- sim + 1
    }
  }
  F_C.numerator   <- (starting_freq - py)^2
  F_C.denominator <- ((starting_freq + py) / 2 - starting_freq * py)
  FC_sim <- F_C.numerator/F_C.denominator
  FC_distribution <- table(sort(FC_sim))
  return(FC_distribution)
}




####### Create filter for loci using a maf threshold
maf_filter <- function(p,MAF_threshold){
  test <- rbind( (p[1,]+p[2,])/2 >= MAF_threshold , (p[1,]+p[2,])/2 <= (1-MAF_threshold)) 
  test <- apply(test,2,all)
  return(test)
}





######## Estimating Effective Population Size (Ne)
EstimateNe.F_C <- function(FC,dT,S1,S2) {
  Ne_hat <- 2*((dT -2)/ (2*(FC - (1/(2*S1)) - (1/(2*S2)))))
  if (Ne_hat<0) Ne_hat <- NA
  return(Ne_hat)
}
EstimateNe.F_ST <- function(Fst_hat,dT) {
  Ne_hat <- dT * (1 - Fst_hat) / (2 * Fst_hat) 
  if (Ne_hat<0) Ne_hat <- NA
  return(Ne_hat)
}


FstatFun_from_dataframe <- function (SNP_list,maf,dT,nbsimul,MAF_threshold) {
  
  S1 <- (SNP_list[1,"derived.x"]+SNP_list[1,"ancestral.x"])/2
  S2 <- (SNP_list[1,"derived.y"]+SNP_list[1,"ancestral.y"])/2
  init_freq_4_test <- SNP_list[,"derived.x"]/(2*S1)
  
  
  
  Fc_list <- matrix(NA,nrow(SNP_list),5)
  colnames(Fc_list)<-c("FC_num","FC_denom","FC_obs","FC_p_value","FC_q_value")
  
  Fc_global <- list(FC_sum_num=NA,FC_sum_denom=NA,FC_multi=NA,Ne_FC_multi=NA)
  
  Fc_list[,"FC_num"]   <- Compute.locus.F_c.num(SNP_list)
  Fc_list[,"FC_denom"] <- Compute.locus.F_c.denom(SNP_list)
  Fc_list[,"FC_obs"]   <- Fc_list[,"FC_num"]/Fc_list[,"FC_denom"]
  
  Fc_global$FC_sum_num   <- sum(Fc_list[maf,"FC_num"])
  Fc_global$FC_sum_denom <- sum(Fc_list[maf,"FC_denom"])
  Fc_global$FC_multi     <- Fc_global$FC_sum_num/Fc_global$FC_sum_denom
  Fc_global$Ne_FC_multi  <- Compute.F_c.N_e(Fc_global$FC_multi,dT,S1,S2)
  
  
  
  max_row_FC <- length(unique(init_freq_4_test)) 
  
  new_list <- list(n_FC=0, Fmat=Fc_list, sim_FC=matrix(NA, nrow=max_row_FC, ncol=nbsimul+2) )
  
  for (locus in 1:nrow(Fc_list) ) {
    new_list <- Drift.simulation.FC(Ne             = round(Fc_global$Ne_FC_multi),
                                    locus          = locus,
                                    new_list       = new_list,
                                    freq           = init_freq_4_test,
                                    nbsimul        = nbsimul,
                                    dT             = dT,
                                    S1             = S1,
                                    S2             = S2,
                                    MAF_threshold  = MAF_threshold)
  }
  
  new_list$Fmat[which(maf==F),4:5] <- NA
  
  return(list(Fc_list=new_list$Fmat,Fc_global=Fc_global))
}














































#	The input data file should be formatted as follows:
#
#	1	1	0	1	2	2	1	0
#	2	1	0	1	0	0	1	2
#	1	2	1	1	0	0	1	2
#	2	2	2	1	1	2	0	2
#	3	2	2	2	1	1	1	1
#
#	Each line corresponds to an individual
# 	The first column contains IDs for individual
# 	The second column contains IDs for sampled demes
# 	The next columns correspond to loci
# 	`0' corresponds to homozygotes for type-1 alleles
# 	`1' corresponds to heterozygotes
# 	`2' corresponds to homozygotes for type-2 alleles
#
#	In the above example, there are two sampled demes made of 2 and 3 diploid individuals, respectively
#	For the above example, Genepop gives the following estimates of F-statistics:
#
#	Multilocus estimates for diploid data
#	Locus           Fwc(is)     Fwc(st)     Fwc(it)
#	------------    -------     -------     -------
#	loc_1            0.0526      0.7543      0.7672
#	loc_2           -0.5652     -0.0376     -0.6241
#	loc_3            0.3793     -0.3232      0.1787
#	loc_4            0.7391     -0.5682      0.5909
#	loc_5           -0.5652     -0.0376     -0.6241
#	loc_6            0.6327     -0.1575      0.5748
#	           All:  0.1847      0.0310      0.2099
#	-----------------------------------------------

compute_Fstats <- function (data,infile=NA,MAF_threshold=0.0) {   # compute the multi-locus estimate of F_ST for diploid data, following Weir (1996)
  
  if (!is.na(infile)) data <- read.table(infile)                  # read the data
  rownames(data) <- NULL																				  # rename the rows
  colnames(data) <- c("ind","pop",paste("loc_",seq(1,ncol(data) - 2),sep = ""))	# rename the columns
  
  gen <- subset(data,select = -c(ind,pop))												# remove the first two columns (ind. ID and sample ID)
  pop <- subset(data,select = pop)																# take the vector of sample ID
  
  lst.pops <- unique(pop)																					# get the list of sample ID
  nbr.pops <- length(lst.pops)																		# compute the number of samples
  nbr.loci <- ncol(gen)																						# compute the number of loci
  counts <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci)						# define the matrix of allele counts
  nbr.hmzgtes.p.1 <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci) 	# define the matrix of type-1 homozygotes
  nbr.hmzgtes.p.2 <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci) 	# define the matrix of type-2 homozygotes
  n <- matrix(NA,nrow = nbr.pops,ncol = nbr.loci)									# define the matrix of sample sizes
  
  cpt <- 1																										
  for (i in lst.pops) {																						          # loop over samples
    counts[cpt,] <- colSums(as.matrix(gen[which(pop == i),]),na.rm=TRUE)		# compute the allele counts from the dataset (per sample and per locus)
    nbr.hmzgtes.p.1[cpt,] <- colSums(as.matrix(gen[which(pop == i),] == 0),na.rm=TRUE)	# compute the number of type-1 homozygotes
    nbr.hmzgtes.p.2[cpt,] <- colSums(as.matrix(gen[which(pop == i),] == 2),na.rm=TRUE)	# compute the number of type-2 homozygotes
    n[cpt,] <- colSums(as.matrix(gen[which(pop == i),] != 9),na.rm=TRUE)								# compute the NUMBER OF INDIVIDUALS (per sample and per locus)
    cpt <- cpt + 1
  }
  
  n. <- colSums(n)																							  # compute the total NUMBER OF INDIVIDUALS (per locus)
  n2 <- colSums(n^2)																							# compute the sum of squared sample sizes (per locus)
  r <- nrow(unique(pop))																					# compute the number of sampled demes
  nc <- (n. - n2 / n.) / (r - 1.0)																# compute the n_c term in Weir (1996)
  
  p.1 <- (2 * n - counts) / (2 * n)																# compute the allele frequency for type-1 alleles (per sample and per locus)
  p.2 <- counts / (2 * n)																					# compute the allele frequency for type-2 alleles (per sample and per locus)
  
  
  vec.p.1.bar <- colSums(2 * n - counts) / colSums(2 * n)					# compute the overall allele frequency for type-1 alleles (per locus)
  vec.p.2.bar <- colSums(counts) / colSums(2 * n)									# compute the overall allele frequency for type-2 alleles (per locus)
  
  p.1.bar <- replicate(r,vec.p.1.bar)													# this is required to perform matrix operations afterwards
  p.2.bar <- replicate(r,vec.p.2.bar)													# this is required to perform matrix operations afterwards
  if (nbr.loci>1){
    p.1.bar <- t(p.1.bar)													# this is required to perform matrix operations afterwards
    p.2.bar <- t(p.2.bar)													# this is required to perform matrix operations afterwards
  }
  
  
  frq.hmzgtes.p.1 <- nbr.hmzgtes.p.1 / n													# compute the frequency of homozygotes of type 1
  frq.hmzgtes.p.2 <- nbr.hmzgtes.p.2 / n													# compute the frequency of homozygotes of type 1
  
  SSG <- colSums(n * (p.1 - frq.hmzgtes.p.1)) + colSums(n * (p.2 - frq.hmzgtes.p.2)) # compute the sum of squares for genes within individuals (Table 5.4, p. 177 in Weir 1996)
  SSI <- colSums(n * (p.1 + frq.hmzgtes.p.1 - 2 * p.1^2)) + colSums(n * (p.2 + frq.hmzgtes.p.2 - 2 * p.2^2)) # compute the sum of squares for individuals within populations (Table 5.4, p. 177 in Weir 1996)
  SSP <- 2 * colSums(n * (p.1 - p.1.bar)^2) + 2 * colSums(n * (p.2 - p.2.bar)^2) # compute the sum of squares for populations (Table 5.4, p. 177 in Weir 1996)
  
  MSG <- SSG / n.																								  # compute the observed mean squares for genes within individuals (Table 5.4, p. 177 in Weir 1996)
  MSI <- SSI / (n. - r)																						# compute the observed mean squares for individuals within populations (Table 5.4, p. 177 in Weir 1996)
  MSP <- SSP / (r - 1.0)																					# compute the observed mean squares for populations (Table 5.4, p. 177 in Weir 1996)
  
  F_ST_locus <-        (MSP - MSI) / (MSP + (nc - 1) * MSI + nc * MSG)   # compute the F_ST at each locus 
  F_IS_locus <-        (MSI - MSG) / (MSI + MSG)                         # compute the F_IS at each locus
  #F_IT_locus <- 1 - (2 * nc * MSG) / (MSP + (nc - 1) * MSI + nc * MSG)   # compute the F_IT at each locus
  

  maf <- maf_filter(p.1,MAF_threshold)
  
  F_ST <- sum(MSP[maf] - MSI[maf]) / sum(MSP[maf] + (nc[maf] - 1) * MSI[maf] + nc[maf] * MSG[maf])		# compute the multilocus F_ST (see Weir 1996, p. 178)
  F_IS <- sum(MSI[maf] - MSG[maf]) / sum(MSI[maf] + MSG[maf])                         # compute the multilocus F_IS

  
  if (length(lst.pops)==2){
    F_C.numerator   <- (p.1[1,] - p.1[2,])^2
    F_C.denominator <- ((p.1[1,] + p.1[2,]) / 2 - p.1[1,] * p.1[2,])
    F_C_locus       <- F_C.numerator/F_C.denominator
    F_C             <- sum(F_C.numerator[maf])/sum(F_C.denominator[maf])

    Fstats <- list(F_ST = F_ST,
                   F_IS = F_IS,
                   F_C  = F_C,
                   F_ST_locus = F_ST_locus,
                   F_IS_locus = F_IS_locus,
                   F_C_locus  = F_C_locus,
                   p = p.1,
                   n = n,
                   maf=maf)
  }else{
    Fstats <- list(F_ST = F_ST,
                   F_IS = F_IS,
                   F_ST_locus = F_ST_locus,
                   F_IS_locus = F_IS_locus,
                   p = p.1,
                   n = n,
                   maf=maf)
  }
  return (Fstats)
}


