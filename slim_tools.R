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

