# How much of the genome is covered?

rep_kpne_genes = names(readDNAStringSet("./input_data/compare_databases/kpne_ag_rep_seqs.fasta"))

rep_kpne_genes <- data.frame(gene_family = rep_kpne_genes)

setDT(rep_kpne_genes)

rep_kpne_genes[, c("gene_family", "gene_product", "number_genomes") := tstrsplit(gene_family, ";", fill = TRUE, keep = c(1,5, 6))]

rep_kpne_genes[, gene_product := gsub("gene_product=", "", gene_product)]

rep_kpne_genes[, number_genomes := as.integer(gsub("number_genomes=", "", number_genomes))]

rep_kpne_genes[, ag_type := fcase(number_genomes >= 456, "core",
                                 default = "ag")]

kpne_blast <- fread("./input_data/compare_databases/kpne_slurm_ncbi_results.tsv")
colnames(kpne_blast) <- c("query","subject","pident","length",
                          "qlen","slen","qcovs","qcovhsp","evalue","bitscore",
                          "staxids", "sscinames", "stitle"
)

# no self hits
kpne_blast <- kpne_blast[!grepl("Klebsiella pneumoniae", sscinames),]

kpne_blast <- kpne_blast[qcovs >= 70]

kpne_best <- kpne_blast[, .(max_pident = max(pident),
                            max_qcovs  = max(qcovs)),
                        by = query]

colnames(kpne_best)[1] = "gene_family"

kpne_best <- merge(rep_kpne_genes, kpne_best,
                   all.x = TRUE, by = "gene_family")

kpne_best[is.na(max_pident), max_pident := 0]


# Basic: bins across the full 0–100 range
kpne_best[, pident_bin := cut(max_pident,
                              breaks = seq(0, 100, by = 1),
                              include.lowest = TRUE,
                              right = TRUE)]



# plot --------------------------------------------------------------------

png(paste0(outdir_fig, "/kpne_perc_hits.png"),
    width  = 25.7,
    height = 17.0,
    units  = "cm",
    res    = 300)    


make_pident_pie <- function(dt, title_prefix) {
  par(mar = c(0, 2, 4, 2))
  # Classify into three categories
  cat <- ifelse(dt$pident_bin == "(99,100]", "High (99-100]",
                ifelse(dt$pident_bin == "[0,1]",   "No hit",
                       "Low < 98"))
  # Preserve a sensible plotting order
  cat <- factor(cat, levels = c("High (99-100]", "Low < 98", "No hit"))
  
  counts <- table(cat)
  n_total <- sum(counts)
  pct <- round(100 * counts / n_total, 1)
  
  labels <- paste0(names(counts), "\n", counts, " (", pct, "%)")
  
  pie(counts,
      labels = labels,
      main   = paste0(title_prefix, " (N = ",format(n_total, big.mark = ","), ")"),
      col    = c("skyblue2","skyblue4","grey20"), 
      border = c("skyblue2","skyblue4","grey20"),
      xpd = TRUE)
}

mat <- matrix(1:6, 
              byrow = TRUE, ncol = 3)
layout(mat, heights = c(2,5))

# 1. All ag_types
make_pident_pie(kpne_best,
                "All genes")

# 2. Core only
make_pident_pie(kpne_best[ag_type == "core"],
                "Core genes")

# 3. Accessory only
make_pident_pie(kpne_best[ag_type == "ag"],
                "Accessory genes")


make_pident_dist <- function(dt, title_prefix) {
  par(mar = c(4, 4, 0, 0))
  
  counts <- table(droplevels(dt$pident_bin))
  n_total <- sum(counts)
  pct <- round(100 * counts / n_total, 1)
  bar_labels <- paste0(counts, "\n(", pct, "%)")
  ylim_max <- max(counts) * 1.15
  
  bp <- barplot(counts,
                main = "",
                xlab = "Max percent identity bin",
                ylab = "Number of gene families",
                col  = "skyblue4",
                border = NA,
                las  = 2,
                cex.names = 0.7,
                ylim = c(0, ylim_max))
  
  if(any(names(counts) %in%  "(99,100]" == TRUE)){
  # Overlay the (99,100] bar in a different colour
  i <- which(names(counts) == "(99,100]")
  rect(xleft   = bp[i] - 0.5,
       xright  = bp[i] + 0.5,
       ybottom = 0,
       ytop    = counts[i],
       col     = "skyblue2",
       border  = NA)
  }
  
  if(any(names(counts) %in%  "[0,1]" == TRUE)){
    # Overlay the (99,100] bar in a different colour
    i <- which(names(counts) == "[0,1]")
    rect(xleft   = bp[i] - 0.5,
         xright  = bp[i] + 0.5,
         ybottom = 0,
         ytop    = counts[i],
         col     = "grey20",
         border  = NA)
  }

}

# 1. All ag_types
make_pident_dist(kpne_best,
                 "all genes")

# 2. Core only
make_pident_dist(kpne_best[ag_type == "core"],
                 "core genes")

# 3. Accessory only
make_pident_dist(kpne_best[ag_type == "ag"],
                 "accessory genes")

dev.off()