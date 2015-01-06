Make.SNP.list_RV <- function (file_in_1,file_in_2,file_fixed,populations,sample_size,samp_ind_name) {
	for (pop in 1:populations) {
    	con <- paste("file_in_",pop,sep="")
    	lines <- readLines(con=get(con))
    
    	pop_line <- which(lines=="Populations:")
    	mut_line <- which(lines=="Mutations:")
    	gen_line <- which(lines=="Genomes:")
    
    	num_of_pop <- mut_line-pop_line-1
    	num_of_mut <- gen_line-mut_line-1
    	pop_size <- as.numeric(strsplit(lines[pop_line+1],split=" ")[[1]][2])
    	assign(x = paste("pop_size_",pop,sep = ""),value = pop_size)
    
    	if (num_of_mut > 0) {
		    mut_table <- read.table(get(paste("file_in_",pop,sep="")),skip=mut_line, nrows=num_of_mut,row.names=1)
    		colnames(mut_table) <- c("type","x","s","h","n_pop")
		   	count <- matrix(0,num_of_mut,2)
		    colnames(count) <- c("derived","ancestral")
		    mut_table <- cbind(mut_table,count)

			nbr_columns <- count.fields(get(con),sep = " ")
			max_columns <- max(nbr_columns)
   			all_data <- read.table(get(con),skip = gen_line,fill = TRUE,col.names = 1:max_columns)
    
		    sampled_ind <- sample(x = seq(pop_size),size = sample_size)
		    assign(x = paste("sampled_ind_",pop,sep = ""),value = sampled_ind)
		    haplo_1 <- (2 * sampled_ind - 1)
    		haplo_2 <- haplo_1 + 1
			lines <- sort(c(haplo_1,haplo_2))
    
    		sampled_data <- all_data[lines,-1]
    		counts <- table(unlist(sampled_data),useNA = "no")
    		mut_table[names(counts),"derived"] <- counts
    		mut_table[,"ancestral"] <- 2 * sample_size - mut_table[,"derived"]
    	} else {
      		mut_table <- matrix(NA,1,7)
    	  	colnames(mut_table) <- c("type","x", "s", "h","n_pop","derived","ancestral")
    	}        
    	assign(x = paste("mut_table_",pop,sep=""),value = mut_table) 

		lines <- readLines(con = file_fixed)
		mut_line <- which(lines == "Mutations:")
		if (length(mut_line) > 0 ) {
    		num_of_mut <- length(lines) - mut_line
    		fmut_table <- matrix()
    		if (num_of_mut > 0) {
				fmut_table <- read.table(file_fixed,skip = mut_line,nrows = num_of_mut,row.names = 1)
				fmut_table[,"V6"] <- c()
				colnames(fmut_table) <- c("type","x","s","h")
				count <- matrix(NA,num_of_mut,3)
      			colnames(count) <- c("n_pop","derived","ancestral")
				count[,"derived"] <- 2 * sample_size
				count[,"ancestral"] <- 0
      			fmut_table <- cbind(fmut_table,count)
    		} else {
				fmut_table <- matrix(NA,nrow=0,length(c("type","x", "s", "h","n_pop","derived","ancestral")))
				colnames(fmut_table) <- c("type","x", "s", "h","n_pop","derived","ancestral")
			}
		} else {
   	 		fmut_table <- matrix(NA,nrow = 0,length(c("type","x", "s", "h","n_pop","derived","ancestral")))
    		colnames(fmut_table) <- c("type","x", "s", "h","n_pop","derived","ancestral")
    	}
    }
    save(sampled_ind_1,sampled_ind_2,mut_table_1,mut_table_2,fmut_table,file = samp_ind_name)

  ############# MERGE AND SORT TABLES ###################################################
  
	colnames(mut_table_1) <- c("type","x","s","h","n_pop.x","derived.x","ancestral.x")
	colnames(mut_table_2) <- colnames(fmut_table) <- c("type","x","s","h","n_pop.y","derived.y","ancestral.y")
  
  rownames(mut_table_1) <- rownames(mut_table_2) <- rownames(fmut_table) <- c()
  Sel_2 <- rbind(mut_table_2,fmut_table)
  
  # take out m2 lines where a conflict will occur with run.seq
  SelEstim_2 <- Sel_2[order(Sel_2$x,decreasing=FALSE),]
  SelEstim_1 <- mut_table_1[order(mut_table_1$x,decreasing=FALSE),]
  m2_1 <- SelEstim_1[which(SelEstim_1[,"type"]=="m2"),]
  SelEstim_1 <- SelEstim_1[-which(SelEstim_1[,"type"]=="m2"),]
  if(length(which(SelEstim_2[,"type"]=="m2")>0)){
    m2_2 <- SelEstim_2[which(SelEstim_2[,"type"]=="m2"),]
    SelEstim_2 <- SelEstim_2[-which(SelEstim_2[,"type"]=="m2"),]
  }
  m2_all <- merge(m2_1, m2_2, by=c("x","type","s","h"), all = TRUE)
  
  #merge Selestim_1 and Selestim_2 and save it
  #SNP_list <-merge(SelEstim_1, SelEstim_2, by=c("x","type","s","h"), all = TRUE)
  run.seq <- function(x) as.numeric(ave(paste(x), x, FUN = seq_along))
  L <- list(SelEstim_1, SelEstim_2)
  L2 <- lapply(L, function(A) cbind(A, run.seq = run.seq(A$x)))
  SNP_list <- Reduce(function(...) merge(..., all = TRUE), L2)[-5]
  #SNP_list2 <-merge(L2[[1]], L2[[2]], by=c("x","type","s","h","run.seq"), all = TRUE)
  
  SNP_list <- rbind(m2_all,SNP_list)
  
  keep_lines <- array(TRUE,nrow(SNP_list))
  #clear useless lines and fill NA according to sampling cases

	keep_lines <- !((is.na(SNP_list[,"derived.x"]) & SNP_list[,"derived.y"] == 0)
				| (is.na(SNP_list[,"derived.y"]) & SNP_list[,"derived.x"] == 0)
				| (is.na(SNP_list[,"derived.x"]) & SNP_list[,"derived.y"] == 0)
				| (SNP_list[,"derived.x"] == 2 * sample_size & SNP_list[,"derived.y"] == 2 * sample_size)
				| (SNP_list[,"derived.x"] == 0 & SNP_list[,"derived.y"] == 0))
				
	replace.x <- ((is.na(SNP_list[,"derived.x"]) & SNP_list[,"derived.y"] > 0))
    SNP_list[replace.x,"derived.x"] <- 0
    SNP_list[replace.x,"ancestral.x"] <- 2 * sample_size
	replace.y <- ((is.na(SNP_list[,"derived.y"]) & SNP_list[,"derived.x"] > 0))
    SNP_list[replace.y,"derived.y"] <- 0
    SNP_list[replace.y,"ancestral.y"] <- 2 * sample_size

	SNP_list <- SNP_list[keep_lines,]

  ######## TRYING TO CLEAN CO EXISTING MUTATIONS ON SAME NUCLEOTIDE ########################
  # Use CleanX
  
  SNP_list <- CleanX(SNPdata=SNP_list)

  return(SNP_list)
}

Drift.simulation.FC_RV <- function(locus,new_list,freq,nbsimul,dT,S1,S2) {
 	Ne <- as.integer(new_list$Fmat[locus,"Ne_FC_multi"])
 	if (is.finite(Ne) && (!is.na(Ne)) && (Ne > 0)) {
		dist <- which((freq[locus,1] == new_list$sim_FC[,1]) & (Ne == new_list$sim_FC[,2]))
		if (length(dist) > 0) {
			FC_sim <- new_list$sim_FC[dist,-c(1,2)]
  			FCobs <- new_list$Fmat[locus,"FC_obs"]
  			p_value <- length(which(FC_sim >= FCobs)) / nbsimul
			new_list$Fmat[locus,"FC_p_value"] <- p_value
		} else {
			px <- rep(freq[locus,1],nbsimul)
  			px[which(px == 0)] <- 1 / (2 * S1)
  			px[which(px == 1)] <- 1 - 1 / (2 * S1)
			FC_sim <- rep(0,nbsimul)
  			py <- rep(0,nbsimul)
  			sim <- 1
  			while (sim <= nbsimul) {
    			f.drift <- px[sim]
    			for (g in 1:dT) {
      				f.drift <- (rbinom(1,Ne,f.drift)) / Ne
    			}
   		   		py[sim] <- rbinom(1,2 * S2,f.drift) / (2 * S2)
    			if (((px[sim] + py[sim]) / 2) >= 0.01 && ((px[sim] + py[sim]) / 2) <= 0.99) {
      				sim <- sim + 1
    			}
  			}
  			qx <- 1 - px
  			qy <- 1 - py 
  			num1 <- (px - py) * (px - py)
  			num2 <- (qx - qy) * (qx - qy)
  			denum1 <- ((px + py) / 2) - (px * py)
  			denum2 <- ((qx + qy) / 2) - (qx * qy)
  			FC_sim <- (1 / 2) * ((num1 / denum1) + (num2 / denum2))
  			FCobs <- new_list$Fmat[locus,"FC_obs"]
	  		p_value <- length(which(FC_sim >= FCobs)) / nbsimul
			new_list$Fmat[locus,"FC_p_value"] <- p_value
			new_sim <- c(freq[locus,1],Ne,FC_sim)
			new_list$n_FC <- new_list$n_FC + 1
			new_list$sim_FC[new_list$n_FC,] <- new_sim
		}
	} else {
		new_list$Fmat[locus,"FC_p_value"] <- NA
	}
  	return(new_list)
}

Drift.simulation.FST_RV <- function(locus,new_list,freq,nbsimul,dT,S1,S2){
	Ne <- as.integer(new_list$Fmat[locus,"Ne_FST_multi"])
	if (is.finite(Ne) && (!is.na(Ne)) && (Ne > 0)) {

    	dist <- which((freq[locus,1] == new_list$sim_FST[,1]) & (Ne == new_list$sim_FST[,2]))
		if (length(dist) > 0) {
				FST_sim <- new_list$sim_FST[dist,-c(1,2)]
  			FSTobs <- new_list$Fmat[locus,"FST_obs"]
  			p_value <- length(which(FST_sim >= FSTobs)) / nbsimul
			new_list$Fmat[locus,"FST_p_value"] <- p_value
		} else {
  			mat <- matrix(NA,nbsimul,4)
  			px <- rep(freq[locus,1],nbsimul)
  			px[which(px == 0)] <- 1 / (2 * S1)
  			px[which(px == 1)] <- 1 - 1 / (2 * S1)
  			mat[,1]<- px
			sim <- 1
  			while (sim <= nbsimul) {
   		 		f.drift <- px[sim]
   	 			for (g in 1:dT) {
      				f.drift <- (rbinom(1,Ne,f.drift)) / Ne
    			}
    			mat[sim,3] <- rbinom(1,2 * S2,f.drift) / (2 * S2)
    			if (((px[sim] + mat[sim,3]) / 2) >= 0.01 && ((px[sim] + mat[sim,3]) / 2) <= 0.99){
      				sim <- sim +1
    			}
  			}
  			mat[,1] <- 2 * S1 * mat[,1]
  			mat[,3] <- 2 * S2 * mat[,3]
 	 		counts <- apply(X = mat,MARGIN = 2,FUN = as.integer)
  			counts[,2] <- 2 * S1 - counts[,1]
  			counts[,4] <- 2 * S2 - counts[,3]
  
			r <- ncol(counts) / 2
  			l <- seq(1,(2 * r),2)
  			ss <- counts[,l] + counts[,(l + 1)]
  			ss2 <- rowSums((counts[,l] + counts[,(l + 1)])^2)
  			n <- rowSums(ss)
  			nc <- (n - ss2 / n) / (r - 1.0);
  			p <- counts[,l] / ss
  			q <- counts[,(l + 1)] / ss
  			pbar <- rowSums(counts[,l]) / rowSums(ss)
  			qbar <- rowSums(counts[,(l + 1)]) / rowSums(ss)
  			SSI <- rowSums(ss * (p - p^2) + ss * (q - q^2))
  			SSP <- rowSums(ss * (p - pbar)^2 + ss * (q - qbar)^2)
  			MSI <- SSI / (n - r);
  			MSP <- SSP / (r - 1.0);
  			FST_sim <- (MSP - MSI) / (MSP + (nc - 1) * MSI)

	  		FSTobs   <- new_list$Fmat[locus,"FST_obs"]	
  			p_value <- length(which(FST_sim >= FSTobs)) / nbsimul
			new_list$Fmat[locus,"FST_p_value"] <- p_value
			new_sim <- c(freq[locus,1],Ne,FST_sim)
			new_list$n_FST <- new_list$n_FST + 1
			new_list$sim_FST[new_list$n_FST,] <- new_sim
  		}
  	} else {
		new_list$Fmat[locus,"FST_p_value"] <- NA
  	}
  	return(new_list)
}

FstatFun_RV <- function (data_file,dT,nbsimul) {
  remove <- seq(nrow(read.table(data_file, skip=2)))
  S1 <- sum(scan(file=data_file,skip=2,nlines=1)[1:2])/2
  S2 <- sum(scan(file=data_file,skip=2,nlines=1)[3:4])/2
  data_freq <- cbind(read.table(file=data_file,skip=2)[,1]/(2*S1),read.table(file=data_file,skip=2)[,3]/(2*S2))
  
  npop <- scan(file=data_file,nlines=1)
  if (npop!=2) print("Error, number of populations must be 2")
  Fmat <- matrix(NA,nrow(data_freq),18)
  colnames(Fmat)<-c("FC_num","FC_denom","FC_sum_num","FC_sum_denom","FC_obs","FC_multi","Ne_FC_multi","FC_p_value","FC_q_value",
                    "FST_num","FST_denom","FST_sum_num","FST_sum_denom","FST_obs","FST_multi","Ne_FST_multi","FST_p_value","FST_q_value")
  
  Fmat[,"FC_num"] <- Compute.locus.F_c.num(data_file)
  Fmat[,"FC_denom"] <- Compute.locus.F_c.denom(data_file)
  
  Fmat[,"FST_num"] <- Compute.numerator.F_ST(data_file)
  Fmat[,"FST_denom"] <- Compute.denominator.F_ST(data_file)
  
  Fmat[,"FC_sum_num"] <- sum(Fmat[,"FC_num"])
  Fmat[,"FC_sum_denom"] <- sum(Fmat[,"FC_denom"])
  
  Fmat[,"FST_sum_num"] <- sum(Fmat[,"FST_num"])
  Fmat[,"FST_sum_denom"] <- sum(Fmat[,"FST_denom"])
  
  Fmat[,"FC_obs"]<- Fmat[,"FC_num"]/Fmat[,"FC_denom"]
  Fmat[,"FST_obs"]<- Fmat[,"FST_num"]/Fmat[,"FST_denom"]
  
  Fmat[,"FC_multi"] <- (Fmat[,"FC_sum_num"]-Fmat[,"FC_num"])/(Fmat[,"FC_sum_denom"]-Fmat[,"FC_denom"])
  Fmat[,"FST_multi"] <- (Fmat[,"FST_sum_num"]-Fmat[,"FST_num"])/(Fmat[,"FST_sum_denom"]-Fmat[,"FST_denom"])
  
  
  
##  !!!! 
  Fmat[,"Ne_FC_multi"] <- mapply(Compute.F_c.N_e, Fmat[,"FC_multi"], dT,S1,S2)
  Fmat[,"Ne_FST_multi"] <-  mapply(Compute.F_ST.N_e,Fmat[,"FST_multi"], dT)
 ### !!!! 

	max_row_FC <- nrow(unique(cbind(data_freq[,1],as.integer(Fmat[,"Ne_FC_multi"]))))
	max_row_FST <- nrow(unique(cbind(data_freq[,1],as.integer(Fmat[,"Ne_FST_multi"]))))
	new_list <- list(n_FC = 0,n_FST = 0,Fmat = Fmat,sim_FC = matrix(-9,nrow = max_row_FC,ncol = (nbsimul + 2)),sim_FST = matrix(-9,nrow = max_row_FST,ncol = (nbsimul + 2)))
	for (i in remove) {
		new_list <- Drift.simulation.FC_RV(locus = i,new_list=new_list,freq=data_freq,nbsimul=nbsimul,dT=dT,S1=S1,S2=S2)
		new_list <- Drift.simulation.FST_RV(locus = i,new_list=new_list,freq=data_freq,nbsimul=nbsimul,dT=dT,S1=S1,S2=S2)
	}

  #   Fmat[,"FC_q_value"] <- qvalue(Fmat[,"FC_p_value"])$qvalues
  #   Fmat[,"FST_q_value"] <- qvalue(Fmat[,"FST_p_value"])$qvalues
	return(new_list$Fmat)
}
