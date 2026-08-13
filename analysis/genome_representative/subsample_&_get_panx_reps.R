
# Read in E. coli/ s. aurues long genomes and then rep.sequences

ew_dir <-"C:/Users/carac/Dropbox/Eyre-Walker_work/adapt_AGs_MGEs/output/data/"

# s.aurues
sa_loci <- readLines(paste0(ew_dir, "rep_pangeno_saureus_2.fna"))
sa_loci <- sa_loci[grepl(">", sa_loci)]
sa_loci <- gsub(">", "", sa_loci)
sa_loci <- do.call(cbind, tstrsplit(sa_loci, ";|\\|", keep = c(1:2)))
sa_loci <- cbind(sa_loci, rep("saureus", nrow(sa_loci)))


# e.coli
ec_loci <- readLines(paste0(ew_dir, "rep_pangeno_ecoli_2.fna"))
ec_loci <- ec_loci[grepl(">", ec_loci)]
ec_loci <- gsub(">", "", ec_loci)
ec_loci <- do.call(cbind, tstrsplit(ec_loci, ";|\\|", keep = c(1:2)))
ec_loci <- cbind(ec_loci, rep("ecoli", nrow(ec_loci)))

annotated_loci = as.data.frame(rbind(sa_loci, ec_loci))

colnames(annotated_loci) = c("alignment","geno_ID", "species")

setDT(annotated_loci)

count_genos <- annotated_loci[,.(n=.N), by=c("species", "geno_ID")]

count_genos <- count_genos[,.(n=.N), by=species]


# Determine new AG genes --------------------------------------------------
pangenome_lng <- unique(fread(paste0(ew_dir, "gene_information_panx_full_lng.csv"), 
                       select = c("alignment", "geno_ID", "species")))

# sample genomes
pangenome_genos =  unique(pangenome_lng[, .(geno_ID, species)])

set.seed(42)

# saureus
sequenced_sa_genos = unique(annotated_loci[species=="saureus", geno_ID])
sampled_sa_genos = sample(setdiff(pangenome_genos[species=="saureus", geno_ID], sequenced_sa_genos),
                          c(485-length(sequenced_sa_genos)))
saureus_genos = c(sequenced_sa_genos, sampled_sa_genos)

# ecoli
sequenced_ec_genos = unique(annotated_loci[species=="ecoli", geno_ID])
sampled_ec_genos = sample(setdiff(pangenome_genos[species=="ecoli", geno_ID], sequenced_ec_genos),
                          c(485-length(sequenced_ec_genos)))
ecoli_genos = c(sequenced_ec_genos, sampled_ec_genos)

panx_accessory_genes = pangenome_lng[geno_ID %chin% c(saureus_genos, ecoli_genos)]

panx_accessory_genes[, number_genomes:= .N, by=c("species", "alignment")]
panx_accessory_genes[, pan_freq := number_genomes/485]

test <- unique(panx_accessory_genes[, .(alignment, number_genomes, species)])
test <- test[,.(n_genes = .N), by = c("species", "number_genomes")]

setorderv(test, cols = c("species", "number_genomes"), c(1,1))
par(mfrow=c(1,2))
with(test[species=="ecoli"],
     plot(number_genomes, n_genes, type = "b", 
          pch = 16, col = "dodgerblue", main = "E. coli, n = 485"))
with(test[species!="ecoli"],
     plot(number_genomes, n_genes, type = "b", 
          pch = 16, col = "firebrick1", main = "S. aureus, n = 485"))

panx_accessory_genes[, pan_type := fcase(pan_freq < 0.99, "ag",
                                         default = "core")]

fwrite(panx_accessory_genes, paste0(outdir_dat, "/subsampled_485_panx_genos_lng.csv"))

# Get AGs sequence -----------------------------------------------
# >g000074_2;representative_genome=SPARK_1282_C1;locus_tag=SPARK_1282_C1_00802;gene_name=sUL1;gene_product=SulP family inorganic anion transporter;number_genomes=466

panx_accessory_genes <- fread(paste0(outdir_dat, "/subsampled_485_panx_genos_lng.csv"))

rep_header <- unique(fread(paste0(ew_dir, "gene_information_panx_full_lng.csv"), 
                              select = c("alignment", "geno_ID", "Gene", "product")))

rep_header <- merge(panx_accessory_genes[pan_type=="ag", .(species, alignment, geno_ID, number_genomes)], 
                    rep_header, all.x = TRUE, by =c("alignment", "geno_ID"))

rep_header[rep_header == ""] <- "NA"

sa_fna <- readDNAStringSet(c(paste0(outdir_dat, "/rep_pangeno_saureus_man_mod.fna")))
ec_fna <- readDNAStringSet(c(paste0(outdir_dat, "/rep_pangeno_ecoli_man_mod.fna")))

to_fix_ec <- readDNAStringSet(c(paste0(outdir_dat, "/rep_pangeno_ecoli_2.fna")))

get_ag_genes = function(dna_string){
  align_info = tstrsplit(names(dna_string), ";", keep=1)[[1]]
  if(align_info %in% unique(rep_header$alignment)){
    locus_tag_info = tstrsplit(names(dna_string), ";|\\|", keep=3)[[1]]
    geno_id_info = tstrsplit(names(dna_string), ";|\\|", keep=2)[[1]]
    
    header_info = rep_header[geno_ID == geno_id_info & alignment == align_info][1,]
    
    fasta_name = paste0(align_info,";representative_genome=",
                        geno_id_info, ";locus_tag=", locus_tag_info,
                        ";gene_name=", header_info$Gene, ";gene_product=",
                        header_info$product, ";number_genomes=", header_info$number_genomes)
    names(dna_string) = fasta_name
    
    return(dna_string)
  }
}


# Get strings and fix gaps ------------------------------------------------
# get fasta and gff directories

sa_fasta_dir <- "C:/Users/carac/Dropbox/Eyre-Walker_work/data/saureus_fastas/"

ec_fasta_dir <- "C:/Users/carac/Dropbox/Eyre-Walker_work/data/ecoli_fastas/"

ew_gffs_dir <- "C:/Users/carac/Dropbox/Eyre-Walker_work/data/ncbi_gffs_downloads/"

# Get AG string sets
sa_ags_string <- foreach(i=1:length(sa_fna),
                         .packages = c("data.table", "Biostrings"),
                         .combine = "c") %do%{
                           get_ag_genes(sa_fna[i])
                         }

# convert to character, then search
bad_sa_idx <- vcountPattern("-", sa_ags_string) > 0 | vcountPattern("N", sa_ags_string) > 0

bad_sa_idx_names <- names(sa_ags_string[bad_sa_idx])
bad_sa_idx_names <- do.call(cbind, tstrsplit(bad_sa_idx_names, ";", keep = c(2:3)))
bad_sa_idx_names <- cbind(bad_sa_idx_names, names(sa_ags_string[bad_sa_idx]))

fixed_sa_genes <- foreach(i = 1:nrow(bad_sa_idx_names),
                          .packages = c("Biostrings", "data.table", "gUtils"), 
                          .combine = c) %do% {
                            
                            print(i)
                            
                            dt_idx <- as.data.table(do.call(rbind,tstrsplit(bad_sa_idx_names[i,], "=", keep = 2)))[,1:2]
                            
                            gr_ginfo <- rtracklayer::import(paste0(ew_gffs_dir, dt_idx[1,1], ".gff"))
                            
                            id <- dt_idx[1,2][[1]]
                            
                            gr_pos <- gr_ginfo[mcols(gr_ginfo)$type == "gene" & mcols(gr_ginfo)$locus_tag == id]
                            
                            if (length(gr_pos) == 0) {
                              gr_pos <- gr_ginfo[
                                (!is.na(mcols(gr_ginfo)$locus_tag) & mcols(gr_ginfo)$locus_tag == id) |
                                  (!is.na(mcols(gr_ginfo)$old_locus_tag) & mcols(gr_ginfo)$old_locus_tag == id)
                              ][1,]
                            }
                            
                            # import FASTA
                            fasta <- readDNAStringSet(paste0(sa_fasta_dir, dt_idx[1,1], ".fasta"))[[1]]
                            
                            # check gr_pos
                            if (length(gr_pos) == 0) {
                              stop(paste("No matching gene found for id:", id, "at row", i))
                            }
                            
                            # get sequence length
                            seq_len <- length(fasta)
                            
                            # get gene coordinates
                            gene_start <- start(gr_pos)
                            gene_end   <- end(gr_pos)
                            
                            # print if out-of-bounds
                            if (gene_start < 1 || gene_end > seq_len) {
                              cat("Out of bounds at i =", i, "id =", id, "\n")
                              cat("FASTA length =", seq_len, " start =", gene_start, " end =", gene_end, "\n")
                            }
                            
                            subseq_gene <- subseq(fasta, start = start(gr_pos), end = end(gr_pos))
                            
                            if(as.character(strand(gr_pos)) == "-"){
                              subseq_gene <- reverseComplement(subseq_gene)
                            }
                            
                            subseq_gene <- DNAStringSet(subseq_gene)
                            
                            names(subseq_gene) <- bad_sa_idx_names[i,3]
                            
                            subseq_gene
                            
                          }

final_sa_ags_string <- c(sa_ags_string[!bad_sa_idx], fixed_sa_genes)

writeXStringSet(final_sa_ags_string, filepath = "./input_data/compare_databases/sa_ag_rep_seqs.fasta", format = "fasta")


# now, ec

ec_ags_string <- foreach(i=1:length(ec_fna),
                         .packages = c("data.table", "Biostrings"),
                         .combine = "c") %do%{
                           get_ag_genes(ec_fna[i])
                         }


# convert to character, then search
bad_ec_idx <- vcountPattern("-", ec_ags_string) > 0 | vcountPattern("N", ec_ags_string) > 0

bad_ec_idx_names <- names(ec_ags_string[bad_ec_idx])
bad_ec_idx_names <- do.call(cbind, tstrsplit(bad_ec_idx_names, ";", keep = c(2:3)))
bad_ec_idx_names <- cbind(bad_ec_idx_names, names(ec_ags_string[bad_ec_idx]))

fixed_ec_genes <- foreach(i = 1:nrow(bad_ec_idx_names),
                          .packages = c("Biostrings", "data.table", "gUtils"),
                          .combine = c) %do% {

                            print(i)

                            dt_idx <- as.data.table(do.call(rbind,tstrsplit(bad_ec_idx_names[i,], "=", keep = 2)))[,1:2]

                            gr_ginfo <- rtracklayer::import(paste0(ew_gffs_dir, dt_idx[1,1], ".gff"))

                            id <- dt_idx[1,2][[1]]

                            gr_pos <- gr_ginfo[mcols(gr_ginfo)$type == "gene" & mcols(gr_ginfo)$locus_tag == id]

                            if (length(gr_pos) == 0) {
                              gr_pos <- gr_ginfo[
                                (!is.na(mcols(gr_ginfo)$locus_tag) & mcols(gr_ginfo)$locus_tag == id) |
                                  (!is.na(mcols(gr_ginfo)$old_locus_tag) & mcols(gr_ginfo)$old_locus_tag == id)
                              ][1,]
                            }

                            # import FASTA
                            fasta <- readDNAStringSet(paste0(ec_fasta_dir, dt_idx[1,1], ".fasta"))[[1]]

                            # check gr_pos
                            if (length(gr_pos) == 0) {
                              stop(paste("No matching gene found for id:", id, "at row", i))
                            }

                            # get sequence length
                            seq_len <- length(fasta)

                            # get gene coordinates
                            gene_start <- start(gr_pos)
                            gene_end   <- end(gr_pos)

                            # print if out-of-bounds
                            if (gene_start < 1 || gene_end > seq_len) {
                              cat("Out of bounds at i =", i, "id =", id, "\n")
                              cat("FASTA length =", seq_len, " start =", gene_start, " end =", gene_end, "\n")
                            }
                            
                            ###HARD STOP
                            if (gene_start < 1 || gene_end > seq_len) {
                              stop(paste("Out of bounds at i =", i,
                                         "id =", id,
                                         "FASTA length =", seq_len,
                                         "start =", gene_start,
                                         "end =", gene_end))
                            }
                            
                            subseq_gene <- subseq(fasta, start = start(gr_pos), end = end(gr_pos))

                            if(as.character(strand(gr_pos)) == "-"){
                              subseq_gene <- reverseComplement(subseq_gene)
                            }

                            subseq_gene <- DNAStringSet(subseq_gene)

                            names(subseq_gene) <- bad_ec_idx_names[i,3]

                            subseq_gene

                          }

final_ec_ags_string <- c(ec_ags_string[!bad_ec_idx], fixed_ec_genes)

# # check:
# 
# dna = ec_fna
# # extract first 3 bases of each sequence
# first_codons <- substr(as.character(dna), 1, 3)
# 
# # find which don't start with ATG or GTG
# bad_idx <- which(!(first_codons %in% c("ATG", "GTG", "TTG")))
# 
# # report indices
# first_codons[bad_idx]
# 
# 



writeXStringSet(final_ec_ags_string, filepath = "./input_data/compare_databases/ec_ag_rep_seqs.fasta", format = "fasta")


# Fix incorrect ec_genes --------------------------------------------------




