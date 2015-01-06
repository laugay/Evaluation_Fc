##############################################################################################################
# function to write #MUTATION TYPES
##############################################################################################################
writeMutation <- function(file="slim_input.txt", number_of_types=1, h=0.5, DFE="f", s=0, shape_alpha=numeric(), 
                               append=F, append_mutation=F){
  
  # Check that values on parameters are OK
  if (any(is.na(match(DFE,c("f","g","e"))))) stop("parameter DFE can take only \"f\", \"e\" and \"g\" as values")
  if (any(is.na(c(h,DFE,s,shape_alpha))))    stop("parameters h, DFE, s and shape_alpha cannot take NA as values")
  
  # Check lengths of h, DFE, S and shape_alpha
  if (length(h)<number_of_types){
    warning("Number of dominance coefficient lower than number of mutation types, values will be reused")
    x<-numeric()
    for (i in 1:ceiling(number_of_types/length(h)))  x <- c(x,h)
    h<-x
  }
  if (length(DFE)<number_of_types){
    warning("Number of distribution of fitness effects (DFE) lower than number of mutation types, values will be reused")
    x<-character()
    for (i in 1:ceiling(number_of_types/length(DFE))) x <- c(x,DFE)
    DFE<-x
  }
  if (length(s)<number_of_types){
    warning("Number of (mean) selection coefficient lower than number of mutation types, values will be reused")
    x<-numeric()
    for (i in 1:ceiling(number_of_types/length(s))) x <- c(x,s)
    s<-x
  }
  if (length(shape_alpha)<length(which(DFE=="g"))){
    warning("Number of (mean) selection coefficient lower than number of mutation types with DFE=\"g\", values will be reused")
    x<-numeric()
    for (i in 1:ceiling(number_of_types/length(shape_alpha))) x <- c(x,shape_alpha)
    shape_alpha<-x
  }
  
  # Transforms 
  x<-array(NA,number_of_types)
  x[which(DFE=="g")]<-shape_alpha
  shape_alpha<-x
  
  if (!append_mutation){
    write("#MUTATION TYPES",file=file,ncolumns=1,append=append)
  }
  
  for (i in 1:number_of_types){
    if (DFE[i]=="f") write( c( paste("m",i,sep=""), h[i], "f", s[i]), file=file, ncolumns=4, append=T)
    if (DFE[i]=="e") write( c( paste("m",i,sep=""), h[i], "e", s[i]), file=file, ncolumns=4, append=T)
    if (DFE[i]=="g") write( c( paste("m",i,sep=""), h[i], "g", s[i], shape_alpha[i]), file=file, ncolumns=5, append=T)
  }
  
}

##############################################################################################################
# function to write #GENOMIC ELEMENT TYPES
##############################################################################################################
writeGenomicElement <- function(file="slim_input.txt", number_of_types=1, mut_type=list("m1"), prop=list(1), 
                                     append=T, append_genomic_element=F){
  
  # Check that values on parameters are OK
  if (length(mut_type)!=length(prop)) stop("parameters mut_type and prop must be of same length")
  for (i in 1:length(mut_type)){
    if (length(mut_type[[i]])!=length(prop[[i]])) stop("elements in parameters mut_type and prop must be of same length")
    if (any(is.na(match(mut_type[[i]], paste("m",1:100,sep=""))))) stop("parameter mut_type can take only \"m1\", \"m2\"... \"m100\" as values")
    if (any(is.na(c(mut_type[[i]],prop[[i]]))))    stop("parameters mut_type and prop cannot take NA as values")
  }
  
  # Check lengths parameters
  if (length(mut_type)<number_of_types){
    warning("List of type of mutations lower than number of mutation types, values will be reused")
    x<-list()
    for (i in 1:ceiling(number_of_types/length(mut_type)))  x <- c(x,mut_type)
    mut_type<-x
  }
  if (length(prop)<number_of_types){
    warning("List of relative proportion of mutations lower than number of mutation types, values will be reused")
    x<-list()
    for (i in 1:ceiling(number_of_types/length(prop)))  x <- c(x,prop)
    prop<-x
  }
  
  if (!append_genomic_element){
    write("#GENOMIC ELEMENT TYPES",file=file,ncolumns=1,append=append)
  }
  
  for (i in 1:number_of_types){
    x <- character()
    for (j in 1:length(mut_type[[i]])) x <- paste(x,mut_type[[i]][j],prop[[i]][j])
    write( c( paste("g",i,sep=""), x), file=file, ncolumns=2, append=T)
  }
  
}

##############################################################################################################
# function to write #CHROMOSOME ORGANIZATION
##############################################################################################################
writeChromosome <- function(file="slim_input.txt", element_type="g1", start=1, end=1000, 
                                        append=T, append_chromosome=F){
  
  # Check that values on parameters are OK
  
  # Check lengths parameters
  if (length(element_type)!=length(start)) stop("start must be of same length as element_type")
  if (length(element_type)!=length(end))   stop("end must be of same length as element_type")
  
  if (!append_chromosome){
    write("#CHROMOSOME ORGANIZATION",file=file,ncolumns=1,append=append)
  }
  
  for (i in 1:length(element_type)){
    write( c( element_type[i], start[i], end[i] ), file=file, ncolumns=3, append=T)
  }
  
}

##############################################################################################################
# function to write #RECOMBINATION RATE
##############################################################################################################
writeRecombination <- function(file="slim_input.txt", interval_end=1000, r=1e-8, 
                                        append=T, append_recombination=F){
  
  # Check that values on parameters are OK
  
  # Check lengths parameters
  if (length(interval_end)!=length(r)) stop("r must be of same length as interval_end")
  
  if (!append_recombination){
    write("#RECOMBINATION RATE",file=file,ncolumns=1,append=append)
  }
  
  for (i in 1:length(interval_end)){
    write( c( interval_end[i], r[i] ), file=file, ncolumns=2, append=T)
  }
  
}

##############################################################################################################
# function to write RECOMBINATION RATE AND CHROMOSOMES
##############################################################################################################
writeRecombinationChrom <- function(file="slim_input.txt", chr_num=8, genome_length=50000000, r=1e-8, 
                                    append=T){
  
  options(scipen=999)
  # Check that values on parameters are OK
  
  # construct nucleotide segments
  int_start <- seq(from=1, to=genome_length,by=floor(genome_length/chr_num))
  int_end <- seq(from=floor(genome_length/chr_num), to=genome_length, by=floor(genome_length/chr_num))
  nuc_seq <- rate_seq <- matrix(NA,nrow=2*chr_num-1,ncol=1)
  nuc_seq[seq(from=1,to=length(nuc_seq),by=2),]<-int_end
  nuc_seq[seq(from=2,to=length(nuc_seq)-1,by=2),]<-int_start[-1]
  rate_seq[seq(from=1,to=length(rate_seq),by=2),]<-r
  rate_seq[seq(from=2,to=length(rate_seq)-1,by=2),]<-0.5
  recomb <- cbind(nuc_seq,rate_seq)
  
  # Check lengths parameters
  #if (length(interval_end)!=length(r)) stop("r must be of same length as interval_end")
  
  write("#RECOMBINATION RATE",file=file,ncolumns=1,append=append)
  write.table(x=recomb, file=file,col.names=FALSE, row.names=FALSE,append=append)
  return(recomb)
}

##############################################################################################################
# function to write #GENE CONVERSION
##############################################################################################################
writeGenerations <- function(file="slim_input.txt", fraction=0, length=0, append=T){
  write("#GENE CONVERSION",file=file,ncolumns=1,append=append)
  write( c(fraction,length), file=file, ncolumns=2, append=T)
}

##############################################################################################################
# function to write #GENERATIONS
##############################################################################################################
writeGenerations <- function(file="slim_input.txt", t=1000, append=T){
  write("#GENERATIONS",file=file,ncolumns=1,append=append)
  write( t, file=file, ncolumns=1, append=T)
}

##############################################################################################################
# function to write #SEED
##############################################################################################################
writeSeed <- function(file="slim_input.txt", seed=123456, append=T){
  write("#SEED",file=file,ncolumns=1,append=append)
  write( seed, file=file, ncolumns=1, append=T)
}

##############################################################################################################
# function to write #MUTATION RATE
##############################################################################################################
writeMutationRate <- function(file="slim_input.txt", u=1e-8, append=T){
  write("#MUTATION RATE",file=file,ncolumns=1,append=append)
  write( u, file=file, ncolumns=1, append=T)
}

##############################################################################################################
# function to write #INITIALIZATION
##############################################################################################################
writeInitioalization <- function(file="slim_input.txt", filename="slim_output.txt", append=T){
  write("#INITIALIZATION",file=file,ncolumns=1,append=append)
  write( filename, file=file, ncolumns=1, append=T)
}

##############################################################################################################
# function to write #DEMOGRAPHY AND STRUCTURE
##############################################################################################################
writeDemography <- function(file="slim_input.txt", type="P", time=1, pop="p1", N=1000, source_pop=character(), target_pop, m, sigma,
                                        append=T, append_demography=F){
  
  # Check that values on parameters are OK
  if (any(is.na(match(type,c("S","P","M","N"))))) stop("parameter type can take only \"P\", \"N\", \"M\" or \"S\" as values")
  if (length(type)!=1) stop("length(type)!=1; only one demographic events can be set at a time")
  
  if (!append_demography){
    write("#DEMOGRAPHY AND STRUCTURE",file=file,ncolumns=1,append=append)
  }

  if (type=="P") write( c(time, "P", pop, N, source_pop), file=file, ncolumns=5, append=T) 
  if (type=="N") write( c(time, "N", pop, N), file=file, ncolumns=4, append=T) 
  if (type=="M") write( c(time, "M", target_pop, source_pop, m), file=file, ncolumns=5, append=T) 
  if (type=="S") write( c(time, "S", pop, sigma), file=file, ncolumns=4, append=T) 

  
  
}

##############################################################################################################
# function to write #OUTPUT
##############################################################################################################
writeOutput <- function(file="slim_input.txt", type="A", time=1000, filename="slim_output.txt",
                        pop, size, MS=F, mut_type,
                        append=T, append_output=F){
  
  # Check that values on parameters are OK
  if (any(is.na(match(type,c("A","R","F","T"))))) stop("parameter type can take only \"A\", \"R\", \"F\" or \"T\" as values")
  if (length(type)!=1) stop("length(type)!=1; only one output item can be set at a time")
  
  if (!append_output){
    write("#OUTPUT",file=file,ncolumns=1,append=append)
  }
  
  if (type=="A") write( c(time, "A", filename), file=file, ncolumns=3, append=T) 
  if (type=="R") write( c(time, "R", pop, size), file=file, ncolumns=4, append=T)  
  if (type=="F") write( c(time, "F"), file=file, ncolumns=2, append=T)  
  if (type=="T") write( c(time, "T", mut_type), file=file, ncolumns=3, append=T)  



}

##############################################################################################################
# function to write #PREDETERMINED MUTATIONS
##############################################################################################################
writePredeterminedMutations <- function(file="slim_input.txt", number_of_mutations=1, time=10, mut_type="m1", x=100,
                        pop="p1", nAA=0, nAa=1,
                        append=T, append_predetermined_mutations=F){
  
  # Check that values on parameters are OK

  if (!append_predetermined_mutations){
    write("#PREDETERMINED MUTATIONS",file=file,ncolumns=1,append=append)
  }
  
  for (i in 1:number_of_mutations){
    write( c(time[i],mut_type[i],x[i],pop[i],nAA[i],nAa[i]), file=file, ncolumns=6, append=T) 
  }

}

##############################################################################################################
# function to converty slim output into structure format
##############################################################################################################
slim2structure <- function(file_in="slim_output.txt",file_out="slim_structure.txt", dist=F,subsample="all",sample_size=NULL){

  lines <- readLines(con=file_in)

  pop_line <- which(lines=="Populations:")
  mut_line <- which(lines=="Mutations:")
  gen_line <- which(lines=="Genomes:")

  num_of_pop <- mut_line-pop_line-1
  num_of_mut <- gen_line-mut_line-1

  num_of_ind <- numeric()
  for (pop in 1:num_of_pop) num_of_ind <- c(num_of_ind,as.numeric(strsplit(lines[pop_line+pop],split=" ")[[1]][2]))
  num_of_gen <- num_of_ind*2
  
  mut_id <- mut_pos <- numeric()
  if (dist) mut_dist <- -1
  for (mut in 1:num_of_mut){
    mut_id   <- c(mut_id,as.numeric(strsplit(lines[mut_line+mut],split=" ")[[1]][1]))
    mut_pos  <- c(mut_pos,as.numeric(strsplit(lines[mut_line+mut],split=" ")[[1]][3]))
    if (mut>1 && dist) mut_dist <- c(mut_dist, mut_pos[mut]-mut_pos[mut-1])
  }

  if (subsample!="all"){
    if (is.null(sample_size)) {
      sample_size <- array(NA,num_of_pop)
      for (pop in 1:num_of_pop){
        if (subsample=="individuals") cat("\n How many individuals do you want to sample in population ",pop,"? ")
        if (subsample=="haplotypes") cat("\n How many haplotypes/gene-copies/haploid-genomes do you want to sample in population ",pop,"? ")
        sample_size[pop] <- as.integer(readLines(n = 1))
      }
    }
    if (subsample=="individuals"){
      sampled_gen<-list()
      for (pop in 1:num_of_pop){
        x <- sample(num_of_ind[pop],sample_size[pop])
        sampled_gen[[pop]] <- sort(c(x*2-1,x*2))
      }  
    }
    if (subsample=="haplotypes"){
      sampled_gen<-list()
      for (pop in 1:num_of_pop) sampled_gen[[pop]] <- sort(sample(num_of_gen[pop],sample_size[pop]))
    }
  }else{
    sampled_gen<-list()
    for (pop in 1:num_of_pop) sampled_gen[[pop]] <- 1:num_of_gen[pop]
  }
  
  
  
  
  
  
  
  write(mut_pos,file=file_out,ncolumns=num_of_mut)
  if (dist) write(mut_dist,file=file_out,ncolumns=num_of_mut,append=T)
  
  for (pop in 1:num_of_pop){
    x<-0
    for (genome in sampled_gen[[pop]]){
      genotype <- array(3,num_of_mut)
      genotype[match(as.numeric(strsplit(lines[gen_line+x+genome],split=" ")[[1]][-1]),mut_id)] <- 4 
      write( c(paste(pop,ceiling(genome/2),sep="_"),pop,genotype), file=file_out,ncolumns=num_of_mut+2,append=T) 
      
    }
    x<-x+num_of_gen[pop]
  }
  
  if(dist){
    return(mut_dist)
  }else{
    return(list(num_of_pop=num_of_pop,num_of_loci=num_of_mut,num_of_ind=num_of_ind))
  }
}

##############################################################################################################
# function to convert slim output into fasta format
##############################################################################################################
slim2fasta <- function(file_in="slim_output.txt",file_out="slim_fasta.txt",invariant_sites=F,ancestral=T,subsample="all",sample_size=NULL){
  require(seqinr)
  
  lines <- readLines(con=file_in)
  
  pop_line <- which(lines=="Populations:")
  mut_line <- which(lines=="Mutations:")
  gen_line <- which(lines=="Genomes:")
  
  num_of_pop <- mut_line-pop_line-1
  num_of_mut <- gen_line-mut_line-1
  
  num_of_ind <- numeric()
  for (pop in 1:num_of_pop) num_of_ind <- c(num_of_ind,as.numeric(strsplit(lines[pop_line+pop],split=" ")[[1]][2]))
  num_of_gen <- num_of_ind*2
  
  mut_id <- mut_pos <- numeric()
  for (mut in 1:num_of_mut){
    mut_id   <- c(mut_id,as.numeric(strsplit(lines[mut_line+mut],split=" ")[[1]][1]))
    if (invariant_sites) mut_pos  <- c(mut_pos,as.numeric(strsplit(lines[mut_line+mut],split=" ")[[1]][3]))
  }
  
  if (subsample!="all"){
    if (is.null(sample_size)) {
      sample_size <- array(NA,num_of_pop)
      for (pop in 1:num_of_pop){
        if (subsample=="individuals") cat("\n How many individuals do you want to sample in population ",pop,"? ")
        if (subsample=="haplotypes") cat("\n How many haplotypes/gene-copies/haploid-genomes do you want to sample in population ",pop,"? ")
        sample_size[pop] <- as.integer(readLines(n = 1))
      }
    }
    if (subsample=="individuals"){
      sampled_gen<-list()
      for (pop in 1:num_of_pop){
        x <- sample(num_of_ind[pop],sample_size[pop])
        sampled_gen[[pop]] <- sort(c(x*2-1,x*2))
      }  
    }
    if (subsample=="haplotypes"){
      sampled_gen<-list()
      for (pop in 1:num_of_pop) sampled_gen[[pop]] <- sort(sample(num_of_gen[pop],sample_size[pop]))
    }
  }else{
    sampled_gen<-list()
    for (pop in 1:num_of_pop) sampled_gen[[pop]] <- 1:num_of_gen[pop]
  }
  
  
  
  
  if (ancestral){
    write( ">ancestral", file=file_out,ncolumns=1,append=F) 
    if (invariant_sites){
      genotype <- array("A",mut_pos[num_of_mut])
      genotype <- c2s(genotype)
    }else{
      genotype <- array("A",num_of_mut)
      genotype <- c2s(genotype)
    }
    write( genotype, file=file_out,ncolumns=1,append=T) 
  }
  first<-T
  for (pop in 1:num_of_pop){
    x<-0
    for (genome in sampled_gen[[pop]]){

      if (ancestral || !first) append<-T
      write( paste(">","pop_",pop,"_geneCopy_",genome,"_ind_",ceiling(genome/2),sep=""), file=file_out,ncolumns=1,append=append) 
      first<-F
      
      if (invariant_sites){
        genotype <- array("A",mut_pos[num_of_mut])
        genotype[mut_pos[match(as.numeric(strsplit(lines[gen_line+x+genome],split=" ")[[1]][-1]),mut_id)]] <- "T" 
        genotype <- c2s(genotype)
        
      }else{
        genotype <- array("A",num_of_mut)
        genotype[match(as.numeric(strsplit(lines[gen_line+x+genome],split=" ")[[1]][-1]),mut_id)] <- "T" 
        genotype <- c2s(genotype)
      }
      write( genotype, file=file_out,ncolumns=1,append=T) 
      
    }
    x<-x+num_of_gen[pop]
  }
  
  return(list(num_of_pop=num_of_pop,num_of_loci=num_of_mut,num_of_ind=num_of_ind))
}

##############################################################################################################
# function theta.k from pegas modified to evaluate values of theta higher than 100
##############################################################################################################
theta.k <- function (x, n = NULL, k = NULL) {
  if (is.null(n)) {
    if (!is.factor(x)) {
      if (is.numeric(x)) {
        n <- sum(x)
        k <- length(x)
      }
      else x <- factor(x)
    }
    if (is.factor(x)) {
      n <- length(x)
      k <- nlevels(x)
    }
  }
  f <- function(th) th * sum(1/(th + (0:(n - 1)))) - k
  uniroot(f, interval = c(1e-08, 1000))$root
}



######################################################################################
# FUNCTIONS TO CLEAN SLIM OUTPUT DELETING MUTATIONS CO-EXISTING ON THE SAME POSITION #
######################################################################################
# By Arnaud Becheler

############# FUNCTION to identify multiple occurences on the same locus
duplicated2 <- function(x) duplicated(x) | duplicated(x, fromLast=TRUE)

############# FUNCTION to test every mutation possibly occuring on the same X position
TeXevery <- function (neutral){
  
  ### dup_lines indicates the rows where there are global duplication 
  dup_lines <- which(duplicated2(neutral[,"x"])==TRUE)
  fact_x <- as.factor(neutral[dup_lines,"x"])
  
  if (length(levels(fact_x))>0) {     # if there are doublons
    # create a TRUE column to turn later on FALSE which lines are to delete
    keepit <- array(TRUE,nrow(neutral))
    neutral<- cbind(neutral,keepit)
    
    # analyze sub-blocks of co-existing mutations
    for (i in seq(from=1,to=length(levels(fact_x)))) {
      position <- levels(fact_x)[i]
      sub_lines <- which(neutral[,"x"]==position)
      neutral[sub_lines,"keepit"]<- duplicated2(neutral[sub_lines,"derived.x"])
      neutral[sub_lines,"keepit"]<- neutral[sub_lines,"keepit"] | duplicated2(neutral[sub_lines,"derived.y"])
    }
    # keep lines with TRUE in "keepit" column, 
    neutral<-neutral[which(neutral[,"keepit"]==TRUE),]
    neutral<-neutral[,colnames(neutral)!="keepit"]
  }else{ 
    print("there are no doublons")
  }
  return(neutral)
}

############# FUNCTION to clean co-occurences using TeXevery
CleanX <- function (SNPdata) {
  
  m2_line <- which(SNPdata["type"]=="m2")
  
  if (length(m2_line)>0){ #if m2 is present in the table
    kept_lines <- which(SNPdata["x"]==SNPdata[m2_line,"x"])
    kept <- SNPdata[kept_lines,]
    neutral <- SNPdata[-kept_lines,]
    
    neutral <- TeXevery(neutral)
    neutral <- rbind(kept,neutral)
    neutral <- neutral[order(neutral$x,decreasing=FALSE),]
    
  }else{ # if m2 is absent in the table
    neutral <- SNPdata
    neutral <- TeXevery(neutral)
    neutral <- neutral[order(neutral$x,decreasing=FALSE),]
  }
  return(neutral)
}

#################### FUNCTION TO READ SLIM OUTPUTS AND TURN THEM INTO AN SNP_list FORM
Make.SNP.list <- function (file_in_1,file_in_2,file_fixed,populations,sample_size,samp_ind_name) {
  mut_table_1 <- mut_table_2 <- matrix()
  
  ##### NON FIXED MUTATIONS ###########################################
  for (pop in 1:populations) {
    con <- paste("file_in_",pop,sep="")
    lines <- readLines(con=get(con))
    
    pop_line <- which(lines=="Populations:")
    mut_line <- which(lines=="Mutations:")
    gen_line <- which(lines=="Genomes:")
    
    num_of_pop <- mut_line-pop_line-1
    num_of_mut <- gen_line-mut_line-1
    pop_size <- as.numeric(strsplit(lines[pop_line+1],split=" ")[[1]][2])
    assign(x=paste("pop_size_",pop,sep=""),value=pop_size)
    
    
    mut_table <- matrix()
    if (num_of_mut > 0) {
      #read mutations informations
      mut_table <- read.table(get(paste("file_in_",pop,sep="")),skip=mut_line, nrows=num_of_mut,row.names=1)
      colnames(mut_table) <- c("type","x", "s", "h","n_pop")
      count <- matrix(0,num_of_mut,2)
      dimnames(count) <- list(c(), c("derived","ancestral"))
      mut_table <- cbind(mut_table,count)
      
      #sample individuals
      sampled_ind <- sample(x=seq(from=1, to=pop_size),size=sample_size,replace=FALSE,prob=NULL)
      assign(x=paste("sampled_ind_",pop,sep=""),value=sampled_ind)
      
      #read line after line
      for (i in seq(from=1, to=sample_size)) {
        ind_line <- gen_line+(2*sampled_ind[i]-1)
        haplo_1 <- read.table(get(paste("file_in_",pop,sep="")),skip=ind_line-1,nrows=1)
        haplo_2 <- read.table(get(paste("file_in_",pop,sep="")),skip=ind_line,nrows=1)
        
        if (length(haplo_1) > 1) {
          mut_table[as.character(haplo_1[-1]),"derived"]<- mut_table[as.character(haplo_1[-1]),"derived"]+1   
        }
        if (length(haplo_2) > 1) {
          mut_table[as.character(haplo_2[-1]),"derived"]<- mut_table[as.character(haplo_2[-1]),"derived"]+1   
        }
      }
      mut_table[,"ancestral"]<- floor(2*sample_size - mut_table[,"derived"])
      assign(x=paste("mut_table_",pop,sep=""),value=mut_table) 
    }else{
      mut_table <- matrix(NA,1,length(c("type","x", "s", "h","n_pop","derived","ancestral")))
      colnames(mut_table) <- c("type","x", "s", "h","n_pop","derived","ancestral")
      assign(x=paste("mut_table_",pop,sep=""),value=mut_table) 
      
    } 
  }

  
  ############# FIXED MUTATIONS ########################################################
  lines <- readLines(con=file_fixed)
  mut_line <- which(lines=="Mutations:")
  if (length(mut_line) > 0 ) {
    num_of_mut <- length(lines) - mut_line
    fmut_table <- matrix()
    if (num_of_mut > 0) {
      
      #read mutations informations
      fmut_table <- read.table(file_fixed,skip=mut_line, nrows=num_of_mut, row.names=1)
      fmut_table[,"V6"]<-c()
      colnames(fmut_table) <- c("type","x", "s", "h")
      
      count <- matrix(NA,num_of_mut,3)
      dimnames(count) <- list(c(), c("n_pop","derived","ancestral"))
      
      count[,"derived"] <- 2*sample_size
      count[,"ancestral"] <- 0
      fmut_table <- cbind(fmut_table,count)
    }else{
      fmut_table<-matrix(NA,nrow=0,length(c("type","x", "s", "h","n_pop","derived","ancestral")))
      colnames(fmut_table) <- c("type","x", "s", "h","n_pop","derived","ancestral")
      assign(x=paste("fmut_table_",pop,sep=""),value=fmut_table) 
    }
  }else{
    fmut_table<-matrix(NA,nrow=0,length(c("type","x", "s", "h","n_pop","derived","ancestral")))
    colnames(fmut_table) <- c("type","x", "s", "h","n_pop","derived","ancestral")
    assign(x=paste("fmut_table_",pop,sep=""),value=fmut_table)
  }
  
  save(sampled_ind_1,sampled_ind_2, mut_table_1, mut_table_2,fmut_table, file=samp_ind_name)
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
  
  del_lines <- array(TRUE,nrow(SNP_list))
  #clear useless lines and fill NA according to sampling cases
  for (j in seq(from=1,to=nrow(SNP_list))) {
    
    if (is.na(SNP_list[j,"derived.x"])) {
      if (SNP_list[j,"derived.y"]==0) {
        del_lines[j] <- FALSE
      }else{
        SNP_list[j,"derived.x"]<- 0
        SNP_list[j,"ancestral.x"]<- 2*sample_size
      }
      
    }else if (is.na(SNP_list[j,"derived.y"])) {
      if (SNP_list[j,"derived.x"]==0) {
        del_lines[j] <- FALSE
      }else{
        SNP_list[j,"derived.y"]<- 0
        SNP_list[j,"ancestral.y"]<- 2*sample_size
      }
      
    }else{
      if (SNP_list[j,"derived.x"]==2*sample_size && SNP_list[j,"derived.y"]==2*sample_size){
        del_lines[j] <- FALSE
      }else  if (SNP_list[j,"derived.x"]==0 && SNP_list[j,"derived.y"]==0){
        del_lines[j] <- FALSE
      }
    }
  }
  SNP_list <- SNP_list[del_lines,]
  
  
  ######## TRYING TO CLEAN CO EXISTING MUTATIONS ON SAME NUCLEOTIDE ########################
  # Use CleanX
  SNP_list <- CleanX(SNPdata=SNP_list)

  return(SNP_list)
}


##### FUNCTION TO TEST IF M2 HAS THE FREQUENCIES OK FOR THE FSTAT TEST...
m2.gbmaf <- function(loc_list){
  m2_test <-FALSE
  ml <- loc_list[which(loc_list[,"type"]=="m2"),]
  nx <- (ml[1,"derived.x"] + ml[1,"ancestral.x"])
  px <- ml[1,"derived.x"]/nx
  ny <- (ml[1,"derived.y"] + ml[1,"ancestral.y"])
  py <- ml[1,"derived.y"]/ny
  if((px+py)/2 >= 0.01 && (px+py)/2 <= 0.99){
    m2_test <- TRUE
  }else{m2_test <- FALSE}
  
  return(m2_test)
}

############################################# FUNCTIONS FOR FC AND FST ANALISIS ######################################
######## FST locus numerator
Compute.numerator.F_ST <- function(infile) {
  counts <- read.table(infile,skip = 2)
  
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
  FST_num <- (MSP - MSI)
  #FST <- sum(MSP - MSI) / sum(MSP + (nc - 1) * MSI)
  return(FST_num)
}

######## FST locus denominator
Compute.denominator.F_ST <- function(infile) {
  counts <- read.table(infile,skip = 2)
  
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
  FST_denom <- (MSP + (nc - 1) * MSI)
  #FST <- sum(MSP - MSI) / sum(MSP + (nc - 1) * MSI)
  return(FST_denom)
}
######## FC (numerator)
Compute.locus.F_c.num <- function(infile) {
  counts <- read.table(infile,skip = 2)
  
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
Compute.locus.F_c.denom <- function(infile) {
  counts <- read.table(infile,skip = 2)
  
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

Compute.F_ST.N_e <- function (mean_Fst,dT) {
  N_e <- (-dT) / log(1-2*mean_Fst)
  return(N_e)
}

Compute.F_c.N_e <- function(mean_FC,dT,S1,S2) {
  N_e      <- 2*((dT -2)/ (2*(mean_FC - (1/(2*S1)) - (1/(2*S2)))))
  return(N_e)
}

######## Simulations for FC, one locus

# locus<- 1
# Fmat<-Fmat
# freq<-data_freq
# nbsimul<-nbsimul
# dT<-dT
# S1<-S1
# S2<-S2

Drift.simulation.FC <- function(locus, Fmat, freq, nbsimul,  dT, S1, S2){
  NeI     <- as.integer(Fmat[locus,"Ne_FC_multi"])
  px <- rep(freq[locus,1],nbsimul)
  px[which(px==0)]<- 1/(2*S1)
  
  ### DRIFT
  FC_sim <- rep(0,nbsimul)
  py <- rep(0,nbsimul)
  sim <- 1
  while (sim  <= nbsimul) {
    f.drift <- px[sim]
    for (g in 1:dT) {
      f.drift <- (rbinom(1,NeI,f.drift))/(NeI)
    }
    if(((px[sim]+f.drift)/2) >= 0.01 && ((px[sim]+f.drift)/2) <= 0.99) {
      py[sim] <- rbinom(1,2*S2,f.drift)/(2*S2)
      sim <- sim +1
    }
  }
  qx <- 1-px
  qy <- 1-py 
  num1 <- (px-py)*(px-py)
  num2 <- (qx-qy)*(qx-qy)
  denum1 <- ((px+py)/2) - (px*py)
  denum2 <- ((qx+qy)/2) - (qx*qy)
  FC_sim <- (1/2)*((num1/denum1) + (num2/denum2))
  FC_sim[which(FC_sim=='NaN')]<- 0
  
  
  ### PVALUE
  FCobs   <- Fmat[locus,"FC_obs"]
  FC_p_value <- length(which(FC_sim>=FCobs))/nbsimul
  return(FC_p_value)
}

######## Simulations for FST, one locus

Drift.simulation.FST <- function(locus, Fmat, freq, nbsimul,  dT, S1, S2){
  NeI     <- as.integer(Fmat[locus,"Ne_FST_multi"])
  mat <- matrix(NA,nbsimul,4)
  px <- rep(freq[locus,1],nbsimul)
  px[which(px==0)]<- 1/(2*S1)
  mat[,1]<- px
  
  ### DRIFT
  sim <- 1
  while (sim <= nbsimul) {
    f.drift <- px[sim]
    for (g in 1:dT) {
      f.drift <- (rbinom(1,NeI,f.drift))/(NeI)
    }
    if(((px[sim]+f.drift)/2) >= 0.01 && ((px[sim]+f.drift)/2) <= 0.99){
      mat[sim,3] <- rbinom(1,2*S2,f.drift)/(2*S2)
      sim <- sim +1
    }
  }
  mat[,1] <- 2*S1*mat[,1]
  mat[,3] <- 2*S2*mat[,3]
  counts <- apply(X=mat,MARGIN=2,FUN=as.integer)
  counts[,2]<- 2*S1 - counts[,1]
  counts[,4]<- 2*S2 - counts[,3]
  
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
  #FST_num <- (MSP - MSI)
  FST <- (MSP - MSI) / (MSP + (nc - 1) * MSI)
  
  #   FST_sim[which(FST_sim=='NaN')] <- 0
  
  ### PVALUE
  FSTobs   <- Fmat[locus,"FST_obs"]
  FST_p_value<- length(which(FST>=FSTobs))/nbsimul
  return(FST_p_value)
}



####################### Read and analyse data file using F functions #######################################

# data_file <- "SV_1_0_selestim_in"
# data_file <- "false_selestim"
# data_file <- "false_selestim_no_fixation"
# data_file <- "false_selestim_maf_2"

# dT<-20
# nbsimul<-100

FstatFun <- function (data_file,dT,nbsimul) {
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
  
  Fmat[,"Ne_FC_multi"] <- mapply(Compute.F_c.N_e, Fmat[,"FC_multi"], dT,S1,S2)
  Fmat[,"Ne_FST_multi"] <-  mapply(Compute.F_ST.N_e,Fmat[,"FST_multi"], dT)
  
  Fmat[,"FC_p_value"] <- mapply(Drift.simulation.FC, locus=remove, MoreArgs=list(Fmat=Fmat, freq=data_freq, nbsimul=nbsimul, dT=dT, S1=S1, S2=S2))
  Fmat[,"FST_p_value"] <- mapply(Drift.simulation.FST, locus=remove, MoreArgs=list(Fmat=Fmat, freq=data_freq, nbsimul=nbsimul, dT=dT, S1=S1, S2=S2))
  
  #   Fmat[,"FC_q_value"] <- qvalue(Fmat[,"FC_p_value"])$qvalues
  #   Fmat[,"FST_q_value"] <- qvalue(Fmat[,"FST_p_value"])$qvalues
  
  return(Fmat)
}
