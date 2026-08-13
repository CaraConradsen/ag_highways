# Kronaplots


# dependencies --------------------------------------------------------------

# load the functions for plotting Krona charts
source('C:/Users/carac/Dropbox/Vos_Lab/ag_highways/analysis/embed_krona-main/embed_krona_wsl.R')

# UHGG
kpne_uhgg_blast_hits <- fread("./input_data/compare_databases/kpne_uhgg_results.tsv")

colnames(kpne_uhgg_blast_hits) <- c("query","genome","percent_identity","length",
                                    "qlen","slen","qcovs","qcovhsp","evalue","bitscore"
)

kpne_uhgg_blast_hits <- kpne_uhgg_blast_hits[percent_identity >=99 & qcovs >90]

kpne_uhgg_blast_hits[,genome := tstrsplit(genome, "_", fill = TRUE, keep = 1)]

uhgg_metadata = fread("C:/Users/carac/Dropbox/Vos_Lab/ag_highways/input_data/genomes-all_metadata.tsv",
                      select = c("Genome", "Species_rep", "Lineage"))
colnames(uhgg_metadata) <- tolower(colnames(uhgg_metadata))

uhgg_metadata[, lineage := gsub("\\b[a-z]__+", "", lineage)]

uhgg_metadata[, c("domain", "phyla", "class", "order", "family", "genus", "species") := tstrsplit(lineage, ";", fill = TRUE, keep = c(1:7))]

uhgg_metadata[, lineage := NULL]

uhgg_metadata[species == "TRUE", species := ""]

kpneu_uhgg_hits = merge(kpne_uhgg_blast_hits[,.(query, percent_identity, genome)],
                        uhgg_metadata, all.x = TRUE, by = "genome")

kpneu_uhgg_hits[species=="", species := NA]

kpneu_uhgg_hits[, counts:=.N, by = c("domain","phyla","class",
                                     "order","family","genus","species")]

kpne_uhgg_krona_dat <- unique(kpneu_uhgg_hits[,.(domain, phyla, class, order, family, genus, species, counts)])

# remove pneumonia
kpne_uhgg_krona_dat <- kpne_uhgg_krona_dat[species!="Klebsiella pneumoniae"]

plot_krona(kpne_uhgg_krona_dat[,1:7], kpne_uhgg_krona_dat[,.(counts)],
           output ="./krona_files/kpne_uhgg.html")

# NCBI Krona plot
kpne_NCBI_blast_hits <- fread("./input_data/compare_databases/kpne_slurm_ncbi_results.tsv")

colnames(kpne_NCBI_blast_hits) <- c("qseqid","sseqid", "percent_identity",
                                    "length", "mismatch", "gapopen", "qstart",
                                    "qend", "sstart", "send", "evalue", "bitscore", "stitle")

kpne_NCBI_blast_hits <- kpne_NCBI_blast_hits[percent_identity >=99 & qcovs >90]

kpne_NCBI_blast_hits[,genome := tstrsplit(genome, "_", fill = TRUE, keep = 1)]


# GTDB --------------------------------------------------------------------
kpne_GTDB_blast_hits <- fread("./input_data/compare_databases/kpne_gtdb_g_results.tsv")

colnames(kpne_GTDB_blast_hits) <- c("query","subject","percent_identity","length",
                                    "qlen","slen","qcovs","qcovhsp","evalue","bitscore"
)

kpne_GTDB_blast_hits <- kpne_GTDB_blast_hits[percent_identity >=99 & qcovs >90]

# tidy names
kpne_GTDB_blast_hits[,c("assmbly") := tstrsplit(subject , "\\|", fill = TRUE, keep = 1)]

# get taxonomy data
gtdb_taxonomy = fread("C:/Users/carac/Dropbox/Vos_Lab/ag_highways/input_data/bac120_taxonomy_r226.tsv",
                      sep = "\t", header = FALSE)
colnames(gtdb_taxonomy) = c("assmbly", "lineage")
gtdb_taxonomy[, lineage := gsub("\\b[a-z]__+", "", lineage)]

gtdb_taxonomy[, c("domain", "phyla", "class", "order", "family", "genus", "species") := tstrsplit(lineage, ";", fill = TRUE, keep = c(1:7))]

gtdb_taxonomy[, assmbly := gsub("GB_|RS_", "", assmbly)] 

gtdb_taxonomy[, lineage := NULL]

kpne_GTDB_blast_hits = merge(kpne_GTDB_blast_hits[,.(query, assmbly)],
                             gtdb_taxonomy, all.x = TRUE, by = "assmbly")

kpne_GTDB_blast_hits <- kpne_GTDB_blast_hits[domain=="Bacteria"]

kpne_GTDB_blast_hits[, counts:=.N, by = c("domain","phyla","class",
                                     "order","family","genus","species")]

kpne_GTDB_blast_hits <- unique(kpne_GTDB_blast_hits[,.(domain, phyla, class, order, family, genus, species, counts)])

# remove pneumonia
kpne_GTDB_blast_hits <- kpne_GTDB_blast_hits[species!="Klebsiella pneumoniae"]

plot_krona(kpne_GTDB_blast_hits[,1:7], kpne_GTDB_blast_hits[,.(counts)],
           output ="./krona_files/kpne_gtdb_g.html")
# NCBI --------------------------------------------------------------------
kpne_ncbi_hits <- fread(, paste0(outdir_dat, "/kpne_ncbi_taxon.csv"))

kpne_ncbi_hits <- kpne_ncbi_hits[domain=="Bacteria"]

kpne_ncbi_hits[, counts:=.N, by = c("domain","phylum","class",
                                          "order","family","genus","species")]

kpne_ncbi_hits <- unique(kpne_ncbi_hits[,.(domain, phylum, class, order, family, genus, species, counts)])
# remove pneumonia
kpne_ncbi_hits <- kpne_ncbi_hits[species!="Klebsiella pneumoniae"]

plot_krona(kpne_ncbi_hits[,1:7], kpne_ncbi_hits[,.(counts)],
           output ="./krona_files/kpne_ncbi.html")
