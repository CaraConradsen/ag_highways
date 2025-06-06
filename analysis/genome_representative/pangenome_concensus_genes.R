
pangenome_data_chr <- fread("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/pangenome_data_chr.csv")

# === Run on all files in directory ===
input_dir <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_1695_out/feature_sequences"  # change to your actual folder
fasta_files <- list.files(input_dir, pattern = "nucleotide.fasta", full.names = TRUE)
fasta_files <- as.data.table(fasta_files)

# Subset to AGs
chunk_size <- 100

ag_names <- unique(pangenome_data_chr[pan_grp != "core", gene_family])# remove singletons
ag_locus_tag <- unique(pangenome_data_chr[pan_grp != "core", .(gene_family, locus_tag)])

chunks <- split(ag_names, ceiling(seq_along(ag_names) / chunk_size))

# Apply each chunk's pattern to subset dt
fasta_files_flt <- lapply(chunks, function(chunk) {
  ag_pattern <- paste(chunk, collapse = "|")
  fasta_files[grepl(ag_pattern, fasta_files, perl = TRUE)]
})

# Combine results into one data.table
fasta_files_flt <- rbindlist(fasta_files_flt, use.names = TRUE, fill = TRUE)

#detect and remove duplicates
fasta_files_flt[, n:=.N, by = fasta_files]

fasta_files_flt <- fasta_files_flt[n <2]

fasta_files_flt$n = NULL

fasta_files_flt <- fasta_files_flt[,fasta_files]

# Get a representative fna (nucleotide) file to blast against genomes ---------------
# Will be doing a nblast, so should be able to include 'N' but not "-"

get_consensus_fna <- function(alignment_list, ag_locus_tag = NULL,
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
                                   
                                   gene_fam = gsub(".nucleotide.fasta", "",basename(f))
                                   
                                   # Read aligned sequences
                                   nt_set <- readDNAStringSet(f)
                                   
                                   if(!is.null(ag_locus_tag)){
                                     nt_set <- nt_set[names(nt_set) %in% ag_locus_tag[gene_family == gene_fam, locus_tag]]
                                   }
                                   
                                   # Find sequences with minimal gaps (count Ns or -)
                                   gap_counts <- vcountPattern("-", nt_set)# vcountPattern("N", nt_set) 
                                   
                                   if (!0 %in% gap_counts) {
                                     return(NULL)
                                   }
                                   
                                   most_complete_seq <- nt_set[gap_counts == min(gap_counts)]
                                   
                                   # Count occurrences
                                   seq_counts <- table(most_complete_seq)
                                   
                                   # Find the most common sequence
                                   most_common_seq <- names(which.max(seq_counts))
                                   
                                   # strore as DNA string set
                                   most_common_seq <- DNAStringSet(most_common_seq)
                                   
                                   # add name
                                   names(most_common_seq) <- gene_fam
                                   
                                   most_common_seq
                                   
                                 }
    
    # Stop cluster after execution
    stopCluster(cl)
  })
  
  writeXStringSet(consensus_results, filepath = paste0(outdir_dat,"/rep_pangeno_",file_name_suffix,".fna"), format = "fasta")
  
}



# input data --------------------------------------------------------------

get_consensus_fna(fasta_files_flt, ag_locus_tag, file_name_suffix = "kleb")# 14456 gene families (including paralogs)




