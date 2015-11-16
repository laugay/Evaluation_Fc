
############################################# FUNCTIONS FOR FC ANALISIS ######################################

######## FC (numerator)
Compute.locus.F_c.num <- function(SNP_from_slim) {
  counts <- SNP_from_slim[,c("derived.x", "ancestral.x",  "derived.y", "ancestral.y")]
  
  r <- ncol(counts) / 2
  l <- seq(1,(2 * r),2)
  ss <- counts[,l] + counts[,(l + 1)]  
  n <- rowSums(ss)
  p <- counts[,l] / ss
  q <- counts[,(l + 1)] / ss
  FC.numerator <- (p[,1] - p[,2])^2 + (q[,1] - q[,2])^2
  return(FC.numerator)
}
######## FC (denominator)
Compute.locus.F_c.denom <- function(SNP_from_slim) {
  counts <- SNP_from_slim[,c("derived.x", "ancestral.x",  "derived.y", "ancestral.y")]
  
  r <- ncol(counts) / 2
  l <- seq(1,(2 * r),2)
  ss <- counts[,l] + counts[,(l + 1)]  
  n <- rowSums(ss)
  p <- counts[,l] / ss
  q <- counts[,(l + 1)] / ss
  FC.denominator <- ((p[,1] + p[,2]) / 2 - p[,1] * p[,2]) + ((q[,1] + q[,2]) / 2 - q[,1] * q[,2])
  return(FC.denominator)
}

######## Computing Effective Size

Compute.F_c.N_e <- function(mean_FC,dT,S1,S2) {
  N_e      <- 2*((dT -2)/ (2*(mean_FC - (1/(2*S1)) - (1/(2*S2)))))
  return(N_e)
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



FIS.compute <- function(H1,H2,sample_size,num_of_pop=1){
  #H1 <- H_list_pop[[1]]
  #H2 <- H_list_pop[[2]]
  
  # matrix of TRUE/FALSE
  homo <- H1==H2
  
  # compute h_tilda, the heterozygosity per locus 
  h_tilda <-  (apply(!homo, MARGIN=1, FUN=sum))/sample_size
  
  # Estimation de FIS = f 
  #n <- c(sample_size, sample_size) 
  n <- sample_size
  r <- num_of_pop
  #  n_bar <- sum(n)/r # the average sample size
  #n_c = (r*mean(n)-sum(n*n/r/mean(n)))/(r-1)
  p_tilda <- ( apply(X=(H1==1), MARGIN=1, FUN=sum) + apply(X=(H2==1), MARGIN=1, FUN=sum) ) / 2/sample_size
  SSI <- 2* (n*p_tilda*(1-p_tilda)) - (1/2)*(n*h_tilda)
  MSI <- SSI/(n-r)
  MSG  <- h_tilda/2
  
  FIS_locus <- (MSI-MSG)/(MSI + MSG)
  
  # somme sur tous les loci
  #  FIS_hat <- (sum(n_c *MSI) - sum(n_c * MSG))/(sum(n_c*MSI)+sum(n_c*MSG))
  FIS_hat <- (sum(MSI) - sum(MSG))/(sum(MSI)+sum(MSG))
  selfing_hat <- 2*FIS_hat/(1+FIS_hat)
  return(list(Fis_all=FIS_locus,WC_Fis=FIS_hat,sigma_hat=selfing_hat))
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

compute_Fstats <- function (data,infile=NA) {                       # compute the multi-locus estimate of F_ST for diploid data, following Weir (1996)
  
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
  for (i in lst.pops) {																						# loop over samples
    counts[cpt,] <- colSums(as.matrix(gen[which(pop == i),]))								# compute the allele counts from the dataset (per sample and per locus)
    nbr.hmzgtes.p.1[cpt,] <- colSums(as.matrix(gen[which(pop == i),] == 0))	# compute the number of type-1 homozygotes
    nbr.hmzgtes.p.2[cpt,] <- colSums(as.matrix(gen[which(pop == i),] == 2))	# compute the number of type-2 homozygotes
    n[cpt,] <- colSums(as.matrix(gen[which(pop == i),] != 9))								# compute the NUMBER OF INDIVIDUALS (per sample and per locus)
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
  
  F_ST <- sum(MSP - MSI) / sum(MSP + (nc - 1) * MSI + nc * MSG)		# compute the multilocus F_ST (see Weir 1996, p. 178)
  F_IS <- sum(MSI - MSG) / sum(MSI + MSG)                         # compute the multilocus F_IS

  Fstats <- list(F_ST=F_ST,
                 F_IS=F_IS,
                 F_ST_locus=F_ST_locus,
                 F_IS_locus=F_IS_locus)
  
  return (Fstats)
}


