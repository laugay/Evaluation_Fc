# SIMULATION OF THE PERIOD WITH SELECTION (in between samples)
if (!quiet) cat(paste(Sys.time(),"Selection period for simulation:",simID,"replicate",replic, "\n"))

# Set file names
slim_in_selection   <- paste0(simID,"_",replic, "_selection.txt")
slim_out_selection  <- paste0(simID,"_",replic, "_selection.out")
slim_log_selection  <- paste0(simID,"_",replic, "_selection.log")
slim_init_selection <- paste0(simID,"_",replic, "_selection_init.txt")

# Read slim output from DRIFT period
out_drift_lines      <- readLines(con=slim_out_drift)
out_drift_pop_line   <- which(out_drift_lines=="Populations:")
out_drift_mut_line   <- which(out_drift_lines=="Mutations:")
out_drift_gen_line   <- which(out_drift_lines=="Genomes:")
out_drift_num_of_mut <- out_drift_gen_line-out_drift_mut_line-1

# Get table of mutations from end of drift period
mutation_table <- read.table(file=slim_out_drift, skip=out_drift_mut_line, nrows=out_drift_num_of_mut)

if (selection_mode!="NM") {
  if (selection_mode!="SV"){
    selection_mode<-"SV"
    warning("Adaptation mode undefined, using standing variation") 
  }
  
  # SELECTION ON STANDING VARIATION (changes selection coefficient of on a random locus)
  
  # chose a random locus
  sampled_mut                   <- sample(x=nrow(mutation_table), size=1)
  # change selection coefficient of locus
  levels(mutation_table[,2])    <- c("m1","m2") 
  mutation_table[sampled_mut,2] <- "m2"
  mutation_table[sampled_mut,4] <- sel_coef
  mutation_table[sampled_mut,5] <- dominance_coef
}

# write initialzation file for slim (state of population at starting point of selection period)
write(out_drift_lines[out_drift_pop_line:out_drift_mut_line], slim_init_selection)
write.table(mutation_table, slim_init_selection, append=T, quote=F, col.names=F, row.names=F)
write(out_drift_lines[-(1:(out_drift_gen_line-1))], slim_init_selection, append=T)
#remove(mutation_table)


# write input file for slim
writeMutation          (file=slim_in_selection, number_of_types=2, h=c(0.5,dominance_coef), DFE=c("f","f"), s=c(0,sel_coef), append=F, append_mutation=F)
writeMutationRate      (file=slim_in_selection, u=u)  
writeGenomicElement    (file=slim_in_selection, number_of_types=1, mut_type=list(c("m1","m2")), prop=list(c(1,0)))
writeChromosome        (file=slim_in_selection, element_type="g1", start=1, end=genome_length)
writeRecombinationChrom(file=slim_in_selection, chr_num=chr_num, genome_length=genome_length, r=r, append=T)  
writeGenerations       (file=slim_in_selection, t=selection_period_duration, append=T)
writeDemography        (file=slim_in_selection, type="S", time=1, pop="p1", sigma=sigma)
writeOutput            (file=slim_in_selection, type="A", time=selection_period_duration, filename=slim_out_selection)
writeOutput            (file=slim_in_selection, type="F", time=selection_period_duration, append_output=T)
writeSeed              (file=slim_in_selection, seed=round(runif(1,-2^31,2^31)) )

write("#INITIALIZATION",file=slim_in_selection,ncolumns=1,append=TRUE)
write(slim_init_selection,file=slim_in_selection,ncolumns=1,append=TRUE)

if (!quiet) cat(paste(Sys.time(),"SIMULATION OF THE PERIOD WITH SELECTION STARTS. Simulation:",simID,"replicate",replic, "\n"))
system(paste("./slim",slim_in_selection,">",slim_log_selection))
if (!quiet) cat(paste(Sys.time(),"END OF SIMULATION OF THE PERIOD WITH SELECTION. Simulation:",simID,"replicate",replic, "\n"))
