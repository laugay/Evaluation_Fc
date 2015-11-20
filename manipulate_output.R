############# FUNCTION to identify multiple occurences on the same locus
duplicated2 <- function(x) duplicated(x) | duplicated(x, fromLast=TRUE)

# Make a table with the info of all polymorphic site at the population level
Make_SNP_table <- function (file_t1,
                            file_t2,
                            file_fixed) {
  # read files
  header <- c("type","x","s","h","pop","time","n_pop")
  for (pop in 1:2) {
    pop_file   <- get(paste("file_t",pop,sep=""))
    file_lines <- readLines( pop_file )
    
    pop_line <- which(file_lines=="Populations:")
    mut_line <- which(file_lines=="Mutations:")
    gen_line <- which(file_lines=="Genomes:")
    
    num_of_pop <- 1 # mut_line-pop_line-1
    num_of_mut <- gen_line-mut_line-1

    pop_size <- as.numeric(strsplit(file_lines[pop_line+1],split=" ")[[1]][2])
    assign(x = paste("pop_size_",pop,sep = ""),value = pop_size)
    remove(file_lines)
    
    if (num_of_mut > 0) {
      mut_table <- read.table( pop_file ,skip=mut_line, nrows=num_of_mut,row.names=1)
      colnames(mut_table) <- header

      max_columns <- max( count.fields( pop_file ,skip = gen_line, sep = " ") )
      all_data <- read.table( pop_file ,skip = gen_line,fill = TRUE, col.names = c("ID",seq_len(max_columns-1)) )

      assign(paste("haplotypes_",pop,sep = ""),all_data)
    } else {
      mut_table <- matrix(NA,1,length(header))
      colnames(mut_table) <- header
      assign(paste("haplotypes_",pop,sep = ""),NA)
    }        
    assign(x = paste("mut_table_",pop,sep=""),value = mut_table) 
  }

  # read files with fixed mutations
  file_lines <- readLines(file_fixed)
  mut_line <- which(file_lines == "Mutations:")
  if (length(mut_line) > 0 ) {
    num_of_mut <- length(file_lines) - mut_line
    fmut_table <- matrix()
    if (num_of_mut > 0) {
      fmut_table <- read.table(file_fixed,skip = mut_line,nrows = num_of_mut,row.names = 1)
      fmut_table[,"V8"] <- c()
      colnames(fmut_table) <- c("type","x","s","h","pop","time")
      count <- matrix(NA,num_of_mut,1)
      colnames(count) <- c("n_pop")
      count[,"n_pop"]     <- 2 * pop_size_2
      fmut_table <- cbind(fmut_table,count)
    } else {
      fmut_table <- matrix(NA,nrow=0,length(header))
      colnames(fmut_table) <- header
    }
  } else {
    fmut_table <- matrix(NA,nrow = 0,length(header))
    colnames(fmut_table) <- header
  }  
  remove(file_lines)
  
  ############# MERGE AND SORT TABLES ###################################################
  
  colnames(mut_table_1) <- c("type","x","s","h","pop","time","n_pop.t1")                         #,"derived.x","ancestral.x")
  colnames(mut_table_2) <- colnames(fmut_table) <- c("type","x","s","h","pop","time","n_pop.t2") #,"derived.y","ancestral.y")

  mut_table_1 <- mut_table_1[as.character(1:nrow(mut_table_1)),]
  mut_table_2 <- mut_table_2[as.character(1:nrow(mut_table_2)),]
  
  haplotypes_1<-as.matrix(haplotypes_1[,2:ncol(haplotypes_1)])
  haplotypes_2<-as.matrix(haplotypes_2[,2:ncol(haplotypes_2)])
  
  haplotypes_1x <- matrix(NA,nrow=nrow(haplotypes_1),ncol=ncol(haplotypes_1))
  for (i in 1:ncol(haplotypes_1) ){
    haplotypes_1x[,i] <-  mut_table_1[haplotypes_1[,i],"x"]
  }
  haplotypes_2x <- matrix(NA,nrow=nrow(haplotypes_2),ncol=ncol(haplotypes_2))
  for (i in 1:ncol(haplotypes_2) ){
    haplotypes_2x[,i] <-  mut_table_2[haplotypes_2[,i],"x"]
  }
  # add fixed derived mutation to haplotypes_2
  if(nrow(fmut_table)>0){
    haplotypes_2x <- cbind(haplotypes_2x,matrix(fmut_table$x,nrow=nrow(haplotypes_2x),ncol=nrow(fmut_table)))
  }

  rownames(mut_table_1) <- rownames(mut_table_2) <- rownames(fmut_table) <- c()
  mut_table_2 <- rbind(mut_table_2,fmut_table)
  
  mut_table_1 <- mut_table_1[order(mut_table_1$x,decreasing=FALSE),]
  mut_table_2 <- mut_table_2[order(mut_table_2$x,decreasing=FALSE),]

  mut_m2_1    <- mut_table_1[which(mut_table_1[,"type"]=="m2"),]
  mut_table_1 <- mut_table_1[-which(mut_table_1[,"type"]=="m2"),]
  
  if(length(which(mut_table_2[,"type"]=="m2"))>0){
    mut_m2_2    <- mut_table_2[which(mut_table_2[,"type"]=="m2"),]
    mut_table_2 <- mut_table_2[-which(mut_table_2[,"type"]=="m2"),]
    mut_m2_2$s  <- mut_m2_1$s
  }else{
    mut_m2_2 <- cbind( mut_m2_1[,c("type","x","s","h","pop","time")],
                       "n_pop.t2"=2 * pop_size_2 )      
  }
  mut_m2 <- merge(mut_m2_1, mut_m2_2, by=c("type","x","s","h","pop","time"), all = TRUE)
  
  mut_table <- merge(mut_table_1, mut_table_2, by=c("type","x","s","h","pop","time"), all = TRUE)
  
  removed_loci <- numeric()
  if( any(duplicated2(mut_table$x)) ){
    positions_to_remove <- which(duplicated2(mut_table$x))
    number_of_loci_multiple_hits <- length(levels(as.factor(mut_table[positions_to_remove,"x"])))
    write(paste("Number of loci removed due to multiple mutations:",number_of_loci_multiple_hits), file=log_file,append=T)
    removed_loci <- mut_table$x[positions_to_remove] 
    mut_table <- mut_table[-positions_to_remove,]
  }
  rownames(mut_table) <- rownames(mut_m2) <- c()
  
  SNP_table <- rbind(mut_m2,mut_table)
  
  return( list(SNP_table=SNP_table,
               m2=mut_m2,
               haplotypes_1=haplotypes_1x,
               haplotypes_2=haplotypes_2x,
               removed_loci=removed_loci) )
}

# makes a sample of individuals
Make_sample <- function (haplotypes_1,
                         haplotypes_2,
                         SNP_table,
                         removed_loci,
                         sample_size,
                         sample_size_loci,
                         Ne_only=F) {
 
  # sample individuals
  for (pop in 1:2) {
    pop_size <- nrow(get(paste0("haplotypes_",pop))) / 2  # in number of individuals 
    sampled_ind <- sample(x = seq(pop_size),size = sample_size[pop])
    haplo_1 <- (2 * sampled_ind - 1)
    haplo_2 <- haplo_1 + 1
    lines_2_sample <- sort(c(haplo_1,haplo_2))
    assign(paste0("sampled_haplotypes_",pop),as.matrix(get(paste0("haplotypes_",pop))[lines_2_sample,]))
  }

  # reduce table of SNP to those polymorhic in sample
  if(ncol(sampled_haplotypes_1)>ncol(sampled_haplotypes_2)){
    col_to_add <- ncol(sampled_haplotypes_1)-ncol(sampled_haplotypes_2)
    sampled_haplotypes_2 <- cbind(sampled_haplotypes_2,matrix(NA,nrow=nrow(sampled_haplotypes_2),ncol=col_to_add))
  }else if(ncol(sampled_haplotypes_1)<ncol(sampled_haplotypes_2)){
    col_to_add <- ncol(sampled_haplotypes_2)-ncol(sampled_haplotypes_1)
    sampled_haplotypes_1 <- cbind(sampled_haplotypes_1,matrix(NA,nrow=nrow(sampled_haplotypes_1),ncol=col_to_add))
  }
  loci_polymorphic_in_sample <- table(rbind(sampled_haplotypes_1,sampled_haplotypes_2))
  loci_polymorphic_in_sample <- loci_polymorphic_in_sample[which(loci_polymorphic_in_sample>0)]
  loci_polymorphic_in_sample <- loci_polymorphic_in_sample[which(loci_polymorphic_in_sample<(sum(sample_size*2)))]
  loci_polymorphic_in_sample <- as.numeric(names(loci_polymorphic_in_sample))
  loci_polymorphic_in_sample <- setdiff(loci_polymorphic_in_sample,removed_loci)
  write(paste("Total number of bi-allelic loci in the sample:",length(loci_polymorphic_in_sample)), file=log_file,append=T)

  SNP_table <- SNP_table[match(loci_polymorphic_in_sample,SNP_table$x ,nomatch=NULL),]

  if (nrow(SNP_table) < sample_size_loci){
    if (!quiet) cat(paste(Sys.time(),"Less than",sample_size_loci,"bi-allelic loci (",nrow(SNP_table),") in sample for simulation",simID,"replicate",replic, "\n"))  
    sample_size_loci <- nrow(SNP_table)
    write(paste("/!\\ Number of sampled loci changed to total number of bi-allelic loci:",nrow(SNP_table)), file=log_file,append=T)
  }else{
    if (!quiet) cat(paste(Sys.time(),"The number of bi-allelic loci in the sample (",nrow(SNP_table),") allows to sample",
                          sample_size_loci,"loci for demographic inference in simulation",simID,"replicate",replic, "\n"))
  }
  
  # sample loci
  if ( sum(SNP_table$type=="m2") == 0 | Ne_only) {
    if (!quiet) cat(paste(Sys.time(),"SNP under selection is absente in sample from simulation",simID,"replicate",replic, "\n"))
    # Sample all loci randomly
    SNP_table <- SNP_table[sort(sample(which(SNP_table[,"type"]=="m1"),sample_size_loci)),]
  }else{
    # sample locus under selection plus (sample_size_loci-1) random loci
    if (!quiet) cat(paste(Sys.time(),"SNP under selection is present in sample from simulation",simID,"replicate",replic, "\n"))
    SNP_table <- SNP_table[sort(c(sample(which(SNP_table[,"type"]=="m1"),sample_size_loci-1),which(SNP_table[,"type"]=="m2"))),]
  }

  loci_position_and_name <- matrix(1:sample_size_loci,nrow=sample_size_loci,ncol=1,dimnames=list(SNP_table$x,"loci"))
  
  pop1hap <- matrix(NA,nrow=sample_size[1]*2,ncol=sample_size_loci)
  pop2hap <- matrix(NA,nrow=sample_size[2]*2,ncol=sample_size_loci)
  for (i in 1:(2*sample_size[1])){
    haplo <- array(0,sample_size_loci)
    derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_1[i,],SNP_table$x)),]
    haplo[ derived_alleles ] <- 1
    pop1hap[i,] <- haplo
  }    
  for (i in 1:(2*sample_size[2])){
    haplo <- array(0,sample_size_loci)
    derived_alleles <- loci_position_and_name[as.character(intersect(sampled_haplotypes_2[i,],SNP_table$x)),]
    haplo[ derived_alleles ] <- 1
    pop2hap[i,] <- haplo
  }
  sampled_haplotypes_1 <- pop1hap
  sampled_haplotypes_2 <- pop2hap
  remove(haplo,pop1hap,pop2hap)
  hap1 <- rbind( as.matrix(sampled_haplotypes_1[((1:sample_size[1])*2)-1,]) , as.matrix(sampled_haplotypes_2[((1:sample_size[2])*2)-1,]) )
  hap2 <- rbind( as.matrix(sampled_haplotypes_1[((1:sample_size[1])*2)  ,]) , as.matrix(sampled_haplotypes_2[((1:sample_size[2])*2)  ,]) )
  M <- array( c(hap1,hap2), dim=c(sum(sample_size),sample_size_loci,2) )
  genotype_data <- apply(M, c(1,2), sum)
  genotype_data <- cbind( seq_len(sum(sample_size)), c(rep(1,each=sample_size[1]),rep(2,each=sample_size[2])), genotype_data )
  remove(M,hap1,hap2)

  return( list(SNP_table=SNP_table,
               genotype_data=genotype_data,
               sampled_haplotypes_1=sampled_haplotypes_1,
               sampled_haplotypes_2=sampled_haplotypes_2) )
}











Drift.simulation.FC <- function(Ne,locus,new_list,freq,nbsimul,dT,S1,S2,MAF_threshold) {
  if (is.finite(Ne) && (!is.na(Ne)) && (Ne > 0)) {
    dist <- which(freq[locus] == new_list$sim_FC[,1]) 
    if (length(dist) > 0) {
      FC_sim <- new_list$sim_FC[dist,-c(1,2)]
      FCobs <- new_list$Fmat[locus,"FC_obs"]
      p_value <- length(which(FC_sim >= FCobs)) / nbsimul
      new_list$Fmat[locus,"FC_p_value"] <- p_value
    } else {
      px <- rep(freq[locus],nbsimul)
      px[which(px == 0)] <- 1/(2*S1)
      px[which(px == 1)] <- 1 - 1/(2*S1)
      
      FC_sim <- rep(0,nbsimul)
      py <- rep(0,nbsimul)
      sim <- 1
      while (sim <= nbsimul) {
        f.drift <- px[sim]
        for (g in 1:dT) {
          f.drift <- (rbinom(1,Ne,f.drift)) / Ne
        }
        y[sim] <- rbinom(1,2 * S2,f.drift) / (2 * S2)
        if (((px[sim] + py[sim]) / 2) >= MAF_threshold && ((px[sim] + py[sim]) / 2) <= (1-MAF_threshold) ) {
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
      new_sim <- c(freq[locus],Ne,FC_sim)
      new_list$n_FC <- new_list$n_FC + 1
      new_list$sim_FC[new_list$n_FC,] <- new_sim
    }
  } else {
    new_list$Fmat[locus,"FC_p_value"] <- NA
  }
  return(new_list)
}

