
# Get a representative fna (nucleotide) file to blast against genomes ---------------

get_consensus_fna <- function(alignment_list, 
                                file_name_suffix = format(Sys.time(), "%Y-%m-%d_%H%M")){
  suppressWarnings({
    # Check if alignment_list is supplied
    if (missing(alignment_list) || length(alignment_list) == 0) {
      warning("No alignment data provided.")
    }
    
    # Set up parallel processing
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    # Paralleled loop to read sequences and compute consensus
    consensus_results <- foreach(f = alignment_list, 
                                 .combine = c, 
                                 .packages = c("Biostrings", "data.table")) %dopar% {
                                   
                                   # Read aligned sequences
                                   nt_set <- readDNAStringSet(f)
                                   
                                   # Find sequences with minimal gaps (count Ns or -)
                                   gap_counts <- vcountPattern("N", nt_set) + vcountPattern("-", nt_set)
                                   most_complete_seq <- nt_set[gap_counts == min(gap_counts)]
                                   
                                   # Count occurrences
                                   seq_counts <- table(most_complete_seq)
                                   
                                   # Find the most common sequence
                                   most_common_seq <- names(which.max(seq_counts))
                                   
                                   # strore as DNA string set
                                   most_common_seq <- DNAStringSet(most_common_seq)
                                   
                                   # get genome name
                                   matches <- which(nt_set == most_common_seq)
                                   
                                   # add name
                                   names(most_common_seq) <- paste(basename(f), names(nt_set[matches[1]]), sep = ";")
                                   
                                   most_common_seq
                                   
                                 }
    
    # Stop cluster after execution
    stopCluster(cl)
  })
  
  # check and replace strings if the start is a gap
  # if string starts with a gap
  bad_strt_idx <- grepl("^-", as.character(consensus_results))
  fix_results <- consensus_results[bad_strt_idx]
  for(i in 1:length(fix_results)){
    algn = tstrsplit(names(fix_results[1]), ";", keep = 1)[[1]]
    nt_set <- readDNAStringSet(alignment_list[grepl(algn, alignment_list)])
    
    bad_strt_idx_sub <- grepl("^-", as.character(nt_set))
    nt_set <- nt_set[!bad_strt_idx_sub]
    
    # Find sequences with minimal gaps (count Ns or -)
    gap_counts <- vcountPattern("N", nt_set) + vcountPattern("-", nt_set)
    most_complete_seq <- nt_set[gap_counts == min(gap_counts)]
    
    # Count occurrences
    seq_counts <- table(most_complete_seq)
    
    # Find the most common sequence
    most_common_seq <- names(which.max(seq_counts))
    
    # strore as DNA string set
    most_common_seq <- DNAStringSet(most_common_seq)
    
    # get genome name
    matches <- which(nt_set == most_common_seq)
    
    # add name
    names(most_common_seq) <- paste(basename(f), names(nt_set[matches[1]]), sep = ";")
    
    most_common_seq
    
  }
  writeXStringSet(consensus_results, filepath = paste0(outdir_dat,"/rep_pangeno_",file_name_suffix,"_2.fna"), format = "fasta")
  
}

# Get genome names --------------------------------------------------------
ecoli_genes_dir <- "C:/Users/carac/Dropbox/Eyre-Walker_work/data/PanX_genes/ecoli/all_genes/"
saureus_genes_dir <- "C:/Users/carac/Dropbox/Eyre-Walker_work/data/PanX_genes/saureus/all_genes/"

# List of input FNA files Saurues
saureus_genes_fnas <- list.files(path = saureus_genes_dir, pattern = "refined_na_aln.fa", full.names = TRUE)

get_consensus_fna(saureus_genes_fnas, file_name_suffix = "saureus")


# List of input FNA files Ecoli
ecoli_genes_fnas <- list.files(path = ecoli_genes_dir, pattern = "refined_na_aln.fa", full.names = TRUE)

get_consensus_fna(ecoli_genes_fnas, file_name_suffix = "ecoli")

