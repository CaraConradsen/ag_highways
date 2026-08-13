# compare across databases and species
# what to get all the data in long


# Create base file --------------------------------------------------------
#kpne
kpne_genes <- names(Biostrings::readDNAStringSet("./input_data/compare_databases/kpne_ag_rep_seqs.fasta")) 
kpne_genes <- gsub("representative_genome=|locus_tag=|gene_name=|gene_product=|number_genomes=", "", kpne_genes)
kpne_genes = data.frame(info = kpne_genes)
setDT(kpne_genes)
kpne_genes[, c("gene_family", "geno_ID", "locus_tag", "gene_name", "gene_product", "number_genomes") :=
             tstrsplit(info, ";", fixed = TRUE)]
kpne_genes[,info:=NULL]
kpne_genes$spp = "kpne"

# ec
ec_genes <- names(Biostrings::readDNAStringSet("./input_data/compare_databases/ec_ag_rep_seqs.fasta")) 
ec_genes <- gsub("representative_genome=|locus_tag=|gene_name=|gene_product=|number_genomes=", "", ec_genes)
ec_genes = data.frame(info = ec_genes)
setDT(ec_genes)
ec_genes[, c("gene_family", "geno_ID", "locus_tag", "gene_name", "gene_product", "number_genomes") :=
             tstrsplit(info, ";", fixed = TRUE)]
ec_genes[,info:=NULL]
ec_genes$spp = "ecoli"

# sa
sa_genes <- names(Biostrings::readDNAStringSet("./input_data/compare_databases/sa_ag_rep_seqs.fasta")) 
sa_genes <- gsub("representative_genome=|locus_tag=|gene_name=|gene_product=|number_genomes=", "", sa_genes)
sa_genes = data.frame(info = sa_genes)
setDT(sa_genes)
sa_genes[, c("gene_family", "geno_ID", "locus_tag", "gene_name", "gene_product", "number_genomes") :=
           tstrsplit(info, ";", fixed = TRUE)]
sa_genes[,info:=NULL]
sa_genes$spp = "saureus"

all_genes = rbind(kpne_genes, ec_genes, sa_genes)
rm(kpne_genes); rm(ec_genes); rm(sa_genes)

# Add COG information ------------------------------------------------------
cogs <- fread("C:/Users/carac/Dropbox/Eyre-Walker_work/adapt_AGs_MGEs/output/data/pangenome_bakta.csv",
              select = c("alignment", "COG"))

kpne_cogs <- list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/long_reads_kpne_bakta_anno/",
                        pattern = ".tsv", recursive = TRUE, full.names = TRUE)

# Get bakta information
kpne_cogs <- kpne_cogs[!grepl("hypotheticals|inference|keep_gffs", kpne_cogs)]

kpne_cogs <- lapply(kpne_cogs, function(cg){
  geno_ID <- gsub(".tsv","", basename(cg))
  cg_dat <- fread(cg, select = c("Locus Tag", "DbXrefs"))

  cg_dat[, COG := {
    # extract COG-prefixed parts
    cogs <- str_extract_all(DbXrefs, "COG:[^, ]+")[[1]]
    # remove "COG:" prefix
    cogs_clean <- sub("^COG:", "", cogs)
    # collapse them into one string, e.g. "COG0593, L"
    if (length(cogs_clean) > 0) paste(cogs_clean, collapse = ", ") else NA_character_
  }, by = 1:nrow(cg_dat)]

  cg_dat$geno_ID = geno_ID
  colnames(cg_dat)[1] = "locus_tag"

  cg_dat[grepl("COG", COG),.(geno_ID, locus_tag, COG)]
})

kpne_cogs <- data.table::rbindlist(kpne_cogs)

# Map to pirate gene families

# Add in UHGG -------------------------------------------------------------
read_uhgg_blast_dat <- function(file_name, coverage = 90, percent_identity = 99){
  dat <- fread(file_name)
  colnames(dat) <- c("query","subject","pident","length",
                     "qlen","slen","qcovs","qcovhsp","evalue","bitscore"
  )
  
  # apply filter
  dat <- dat[pident >= percent_identity & qcovs > coverage]
  dat[,c("gene_family", "locus_tag") := tstrsplit(query, ";", fill = TRUE, keep = c(1,2))]
  dat[,locus_tag := gsub("representative_genome=", "", locus_tag)]
  dat[,db:="uhgg"]
  
  return(dat[,.(gene_family,db, subject, qlen, pident)]) # locus_tag
}

uhgg_file_loc = c("./input_data/compare_databases/kpne_uhgg_results.tsv",
                  "./input_data/compare_databases/ec_uhgg_results.tsv",
                  "./input_data/compare_databases/sa_uhgg_results.tsv")

# first get a feel for all data
uhgg_hits = lapply(uhgg_file_loc, read_uhgg_blast_dat, coverage = 0 , percent_identity = 0)
uhgg_hits <- data.table::rbindlist(uhgg_hits)
uhgg_hits[, subject := sub("_.*", "", subject)]

# add loci information

uhgg_metadata = fread("C:/Users/carac/Dropbox/Vos_Lab/ag_highways/input_data/genomes-all_metadata.tsv",
                      select = c("Genome", "Lineage"))
colnames(uhgg_metadata) <- tolower(colnames(uhgg_metadata))

uhgg_metadata[, lineage := gsub("\\b[a-z]__+", "", lineage)]

uhgg_metadata[, c("domain", "phyla", "class", "order", "family", "genus", "species") := tstrsplit(lineage, ";", fill = TRUE, keep = c(1:7))]

uhgg_metadata[, lineage := NULL]

uhgg_metadata[species == "TRUE", species := ""]

colnames(uhgg_metadata)[1] = "subject"

kpneu_uhgg_hits = merge(uhgg_hits,
                        uhgg_metadata, all.x = TRUE, by = "subject")

kpneu_uhgg_hits[species=="", species := NA]



# Add in NCBI and get taxid list for Kronaplots ---------------------------
read_ncbi_blast_dat <- function(file_name, coverage = 90, percent_identity = 99){
  dat <- fread(file_name)
  colnames(dat) <- c("query","subject","pident","length",
                     "qlen","slen","qcovs","qcovhsp","evalue","bitscore",
                     "staxids", "sscinames", "stitle"
  )
  # apply filter
  dat <- dat[pident >= percent_identity & qcovs > coverage]
  dat[, gene_family := query]
  dat[,db:="ncbi"]
  
  return(dat[,.(gene_family, db, subject, qlen,pident, staxids, sscinames)])
}

ncbi_file_loc = c("./input_data/compare_databases/kpne_slurm_ncbi_results.tsv",
                  "./input_data/compare_databases/ec_slurm_ncbi_results.tsv",
                  "./input_data/compare_databases/sa_slurm_ncbi_results.tsv")

# first get a feel for all data
ncbi_hits = lapply(ncbi_file_loc, read_ncbi_blast_dat, coverage = 0 , percent_identity = 0)
ncbi_hits <- data.table::rbindlist(ncbi_hits)
ncbi_hits[, subject := sub("_.*", "", subject)]

ncbi_krona_taxids <- ncbi_hits[,.(gene_family, db, staxids)]

# split and expand into new rows
ncbi_krona_taxids <- ncbi_krona_taxids[, .(staxids = unlist(strsplit(staxids, ";"))), by = c("db","gene_family")]

ncbi_krona_taxids$staxids = as.numeric(ncbi_krona_taxids$staxids)

taxids <- unique(ncbi_krona_taxids$staxids)  # E. coli, K. pneumoniae, S. aureus
# Split by ";", unlist, and convert to numeric

lineages <- taxize::classification(taxids, db="ncbi")

# Convert to a data.frame

lineage_df <- rbindlist(lapply(lineages, function(x) {
  # make a row data.table from the lineage
  dat <- data.frame(t(x$name), stringsAsFactors = FALSE)
  colnames(dat) <- x$rank
  dat$staxids <- if ("species" %in% x$rank) as.numeric(x$id[x$rank == "species"]) else NA
  
  setDT(dat)
  
  # required columns
  col_nam <- c("staxids", "domain","kingdom","phylum","class",
               "order","family","genus","species")
  
  # add missing cols as NA
  missing_cols <- setdiff(col_nam, names(dat))
  for (mc in missing_cols) dat[, (mc) := NA]
  
  # enforce consistent column order
  dat <- dat[, ..col_nam]
  return(dat)
}))

ncbi_hits <- ncbi_hits[,.(gene_family, db, staxids, qlen, pident)]
ncbi_hits <- ncbi_hits[, .(staxids = unlist(strsplit(staxids, ";"))), by = c("gene_family", "db", "qlen", "pident")]

ncbi_hits$staxids = as.numeric(ncbi_hits$staxids)

ncbi_hits = merge(ncbi_hits,
                  unique(lineage_df), all.x = TRUE, by = "staxids")


# now dat for krona_plots

ncbi_krona_taxids <- merge(unique(ncbi_krona_taxids), unique(lineage_df), 
                           all.x = TRUE, by = "staxids")

ncbi_krona_taxids <- ncbi_krona_taxids[!is.na(genus)]

# fwrite(ncbi_krona_taxids[!grepl("na_aln.fa", gene_family)], paste0(outdir_dat, "/kpne_ncbi_taxon.csv"))
colnames(ncbi_krona_taxids)[6] = "phyla"

# GTDB --------------------------------------------------------------------
read_gtdb_blast_dat <- function(file_name, coverage = 90, percent_identity = 99){
  dat <- fread(file_name)
  colnames(dat) <- c("query","subject","pident","length",
                     "qlen","slen","qcovs","qcovhsp","evalue","bitscore"
  )
  # apply filter
  dat <- dat[pident >= percent_identity & qcovs > coverage]
  dat[, gene_family := query]
  dat[, gene_family :=sub(";.*", "", gene_family)]
  dat[,db:="gtdb"]
  dat[,subject:= tstrsplit(subject, "\\|", keep = 1)]
  
  return(dat[,.(gene_family, db, subject, qlen, pident)])
}

gtdb_file_loc = c("./input_data/compare_databases/kpne_gtdb_g_results.tsv",
                  "./input_data/compare_databases/ec_gtdb_g_results.tsv",
                  "./input_data/compare_databases/sa_gtdb_g_results.tsv")


# first get a feel for all data
gtdb_hits = lapply(gtdb_file_loc, read_gtdb_blast_dat, coverage = 0 , percent_identity = 0)
gtdb_hits <- data.table::rbindlist(gtdb_hits)
gtdb_hits[, subject := sub("_.*", "", subject)]

gtdb_taxonomy = fread("C:/Users/carac/Dropbox/Vos_Lab/ag_highways/input_data/bac120_taxonomy_r226.tsv",
                      sep = "\t", header = FALSE)
colnames(gtdb_taxonomy) = c("subject", "lineage")
gtdb_taxonomy[, lineage := gsub("\\b[a-z]__+", "", lineage)]

gtdb_taxonomy[, c("domain", "phyla", "class", "order", "family", "genus", "species") := tstrsplit(lineage, ";", fill = TRUE, keep = c(1:7))]

gtdb_taxonomy[, subject := gsub("GB_|RS_", "", subject)] 

gtdb_taxonomy[, lineage := NULL]

gtdb_hits = merge(gtdb_hits,
                  gtdb_taxonomy, all.x = TRUE, by = "subject")



# Put everything together -------------------------------------------------
keep_cols = c("gene_family", "db", "domain", "phyla", 
              "class", "order", "family", "genus", "species")
all_blast_hits <- rbind(ncbi_krona_taxids[, ..keep_cols],
                        gtdb_hits[, ..keep_cols],
                        kpneu_uhgg_hits[, ..keep_cols])

all_spp_db_dat <- merge(all_genes,all_blast_hits[domain =="Bacteria"],
                        all.x = TRUE,
                        by ="gene_family") 

fwrite(all_spp_db_dat, paste0(outdir_dat, "all_spp_db_dat.csv"))



# Hits barplots -----------------------------------------------------------
every_uhgg_hit <- merge(all_genes,
                        kpneu_uhgg_hits[, .(gene_family, pident, species)],
                        all.x = TRUE, by = "gene_family")
every_uhgg_hit[, rmv := fcase(
  spp == "kpne" & grepl("Klebsiella pneumoni", species), 1,
  spp == "ecoli" & grepl("Escherichia coli", species), 1,
  spp == "saureus" & grepl("Staphylococcus aureus", species), 1,
  default = 0
)]
every_uhgg_hit[rmv == 1, pident:= 0]

every_gtdb_hit <- merge(all_genes,
                        gtdb_hits[, .(gene_family, pident, species)],
                        all.x = TRUE, by = "gene_family")
every_gtdb_hit[, rmv := fcase(
  spp == "kpne" & grepl("Klebsiella pneumoni", species), 1,
  spp == "ecoli" & grepl("Escherichia coli", species), 1,
  spp == "saureus" & grepl("Staphylococcus aureus", species), 1,
  default = 0
)]
every_gtdb_hit[rmv == 1, pident:= 0]


every_ncbi_hit <- merge(all_genes,
                        ncbi_hits[, .(gene_family, pident, species)],
                        all.x = TRUE, by = "gene_family")
every_ncbi_hit[, rmv := fcase(
  spp == "kpne" & grepl("Klebsiella pneumoni", species), 1,
  spp == "ecoli" & grepl("Escherichia coli", species), 1,
  spp == "saureus" & grepl("Staphylococcus aureus", species), 1,
  default = 0
)]
every_ncbi_hit[rmv == 1, pident:= 0]


get_db_bins <- function(hit_dt,db_name, breaks = c(0,1,80, seq(82, 100, by = 2))){
  hit_dt[is.na(pident), pident := 0]
  hit_dt <- hit_dt[, .(max_pident = max(pident)), by = c("spp", "gene_family")]
  
  # Create a binned column
  hit_dt[, bin := cut(max_pident, breaks = breaks, include.lowest = TRUE, right = TRUE)]
  
  # Count how many in each bin
  hit_bins <- hit_dt[, .N, by = c("spp", "bin")][order(bin)]
  hit_bins[, perc := 100 * N / sum(N), by = spp]
  
  hit_bins$db = db_name
  
  return(hit_bins)
}

gtdb_bins <- get_db_bins(every_gtdb_hit, "gtdb")
ncbi_bins <- get_db_bins(every_ncbi_hit, "ncbi")
uhgg_bins <- get_db_bins(every_uhgg_hit, "uhgg")

all_bins = rbind(gtdb_bins, ncbi_bins, uhgg_bins)

cols12 <- c(
  "#00152C", "#012E52", "#01446B", "#015B83", "#01719B","#0286AC",
   "#049AB5", "#24ADB0", "#4BBF9C", "#78CF80", "#A6DB6D", "#D2E65E"
)

create_stacked_barchart <- function(db_bins_dt, plot_ord = c("uhgg", "ncbi", "gtdb")){
  
  hit_mat <- sapply(plot_ord, function(ord){
    db_bins_dt[db == ord][order(bin), perc]
  })

  # Create stacked bar chart and capture bar midpoints
  bp <- barplot(hit_mat,
                horiz = TRUE,
                beside = FALSE,         # stacked bars
                col = cols12,
                border = "white",
                legend.text = TRUE,
                main = "",
                yaxt = "n")             # suppress default y-axis
  
  # Now add custom y-axis with correct alignment
  axis(side = 2, at = bp, labels = plot_ord, las = 2, tick = FALSE)
}



mat <- matrix(c(1,2,5,
                1,3,5,
                1,4,5,
                6,6,6), nrow = 4, ncol = 3, byrow = TRUE)
layout(mat, widths = c(1,6,1.25),
       heights = c(4,4,4,1), respect = FALSE)

plot.new()
par(mar = c(3,3,0.5,1))
create_stacked_barchart(all_bins[all_bins$spp=="kpne"])
axis(expression(italic(K.~pneumoniae)), side = 2, at = 1.85, 
     xpd = TRUE, line = 2.75, cex.axis = 1.25,
     tick = FALSE, las = 2)
text(98,3.15, "# AGs is 21,455", pos = 2)
create_stacked_barchart(all_bins[all_bins$spp=="ecoli"])
axis(expression(italic(E.~coli)), side = 2, at = 1.85, 
     xpd = TRUE, line = 2.75, cex.axis = 1.25,
     tick = FALSE, las = 2)
text(98,3.15, "# AGs is 26,160", pos = 2)
create_stacked_barchart(all_bins[all_bins$spp=="saureus"])
axis("Percent of accessory genes", side  = 1, at = 50,
     xpd = TRUE, line = 2, tick = FALSE, cex.axis = 1.25)
axis(expression(italic(S.~aureus)), side = 2, at = 1.85, 
     xpd = TRUE, line = 2.75, cex.axis = 1.25,
     tick = FALSE, las = 2)
text(98,3.15, "# AGs is 4,515", pos = 2)

# legend
par(mar = c(13,6,13,3))
labs <- as.character(unique(all_bins[,bin]))
labs <- gsub("\\[|\\]|\\(|\\)", "", labs)
labs <- gsub(",", "-", labs)
labs <- paste0(labs, "%")

plot(c(0, 1), c(0,12), type = "n", xlab = "", ylab = "",
     main = "", xaxt ="n", yaxt ="n", 
     yaxs="i", xaxs="i", xpd = TRUE)
i <- 1*(0:11)
rect(0, 0+i, 1, 1+i, col = cols12, border = "white")
axis(2, at = 0.5 + 0:11, labels = labs, las = 2, cex.axis = 0.8)
box(col="white", lwd = 2)
text(0.5,13, 
     "Maximum %\nof identical\nnucleotide identity\n per AG", 
     xpd = TRUE, cex = 1.1)
