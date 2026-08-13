# look at UHGG Hits


blast_dir <- list.files("C:/Users/carac/Dropbox/Vos_Lab/ag_highways/input_data/blast_output",
                        full.names = TRUE)

blast_columns <- c("Query", "Subject", "Percent_Identity", "Alignment_Length", 
                   "Mismatches", "Gaps", "Gap_Openings", "Subject_Start", 
                   "Subject_End", "Query_Length", "E_value", "Bit_Score")

# register parallel
cl <- makeCluster(num_cores)
registerDoParallel(cl)

blast_protein_results <- foreach(i = 1:length(blast_dir),
                                 .combine = rbind,
                                 .packages = c("stringr", "data.table")) %dopar% {
                                   
                                   data_name = str_split(basename(blast_dir[i]), pattern="_")[[1]][2]
                                   
                                   temp_blast_dt <- fread(blast_dir[i])
                                   
                                   colnames(temp_blast_dt) = blast_columns
                                   
                                   temp_blast_dt[, db_grp := data_name]
                                   
                                   # bit score selection
                                   temp_blast_dt <- temp_blast_dt[, .SD[Bit_Score == max(Bit_Score)], by = Query]
                                   
                                   # lowest evalue
                                   temp_blast_dt <- temp_blast_dt[, .SD[E_value == min(E_value)], by = Query]
                                   
                                   # percent identity
                                   temp_blast_dt <- temp_blast_dt[, .SD[Percent_Identity == max(Percent_Identity)], by = Query]
                                   
                                   # alignment length
                                   temp_blast_dt <- temp_blast_dt[, .SD[Alignment_Length == max(Alignment_Length)], by = Query]
                                   
                                   # alignment length
                                   temp_blast_dt <- temp_blast_dt[, .SD[Query_Length == max(Query_Length)], by = Query]
                                   
                                   # count number of subjects
                                   temp_blast_dt <- temp_blast_dt[, row_N:=.N, by=c("Query")]
                                   temp_blast_dt <- temp_blast_dt[, N:=.N, by=c("Query", "Subject")]
                                   
                                   # get the most frequent
                                   temp_blast_dt_freq <- temp_blast_dt[row_N>1, .SD[which.max(N)], by = Query]
                                   
                                   temp_blast_dt <- rbind(temp_blast_dt_freq, temp_blast_dt[row_N<2,])
                                   
                                   temp_blast_dt[,c("N", "row_N"):= NULL]
                                   
                                   temp_blast_dt
                                   
                                   
                                 }

# Stop the parallel backend after the loop is done
stopCluster(cl)

blast_protein_results[, Genome := tstrsplit(Subject, "_", fill = TRUE, keep = 1)]
blast_protein_results[, Subject := NULL]

fwrite(blast_protein_results, paste0(outdir_dat, "/blast_protein_results.csv"))



# add species information -------------------------------------------------
uhgg_metadata = fread("C:/Users/carac/Dropbox/Vos_Lab/ag_highways/input_data/genomes-all_metadata.tsv",
                      select = c("Genome", "Species_rep", "Lineage"))

uhgg_metadata[, Lineage := gsub("\\b[a-z]__+", "", Lineage)]

uhgg_metadata[, c("domain", "phyla", "class", "order", "family", "genus", "species") := tstrsplit(Lineage, ";", fill = TRUE, keep = c(1:7))]

uhgg_metadata[, Lineage := NULL]

uhgg_metadata[species == "TRUE", species := ""]

kpneu_uhgg_hits = merge(blast_protein_results[,.(Query, Percent_Identity, Genome)],
                        uhgg_metadata, all.x = TRUE, by = "Genome")


# plots -------------------------------------------------------------------

with(kpneu_uhgg_hits[species != 'Klebsiella pneumoniae'],
     hist(Percent_Identity, breaks = 1000))




# Diversity plot ----------------------------------------------------------

diversity_plot <- kpneu_uhgg_hits[species != 'Klebsiella pneumoniae' & Percent_Identity > 95]

# determine genus by family
genus <- diversity_plot[, .(genus_n = .N), by = c("family", "genus")]
genus[, genus_frac := genus_n/ sum(genus_n)]# Genus pie (outer ring), nested within families

setorderv(genus, c("family", "genus_n"), c(1,-1))

# Prepare family totals
fam <- genus[, .(family_n = .N), by = family]
fam[, family_col := rainbow(nrow(fam), alpha = 0.6)]

# merge together
diversity_plot <- merge(genus,fam,  by = "family", all.y = TRUE)
diversity_plot[, family_frac := sum(genus_frac), by = family]

# Add colours for the genus

# Helper function to darken a colour (alpha preserved)
darken_colours <- function(base_col, n) {
  rgb_vals <- col2rgb(base_col, alpha = TRUE) / 255
  # Darkening factors from 1 (original) to ~0.5 (darker)
  dark_factors <- seq(1, 0.5, length.out = n)
  apply(matrix(dark_factors), 1, function(f) {
    rgb(rgb_vals[1] * f, rgb_vals[2] * f, rgb_vals[3] * f, alpha = rgb_vals[4])
  })
}

# Assign darker genus colours within each family
diversity_plot[, genus_col := darken_colours(family_col[1], .N), by = family]

# Start plotting

layout(matrix(c(1,2),nrow=1),
       width = c(3,1))

par(mar = c(2, 1, 1, 2), xpd=TRUE)
plot.new()
plot.window(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))

start <- 0
outer_start <- 0
for(f in unique(diversity_plot$family)){
  # family
  radius_inner = 0.0; radius_outer = 0.5
  end <- start + diversity_plot[family == f, family_frac[1]] * 2 * pi
  theta <- seq(start, end, length.out = 100)
  x <- c(radius_inner * cos(rev(theta)), radius_outer * cos(theta))
  y <- c(radius_inner * sin(rev(theta)), radius_outer * sin(theta))
  polygon(x, y, col = diversity_plot[family == f, family_col[1]], border = "black")
  
  if(diversity_plot[family == f, family_frac[1]] > 0.1){
    text(mean(x), max(y)/2, cex= 0.75,
         diversity_plot[family == f, family[1]])
  }
  
  # Outer genus loop
  radius_inner = 0.5; radius_outer = 1
  
  for(g in diversity_plot[family == f, genus]){
    outer_end <- outer_start + diversity_plot[family == f & genus == g, genus_frac] * 2 * pi
    theta <- seq(outer_start, outer_end, length.out = 100)
    x <- c(radius_inner * cos(rev(theta)), radius_outer * cos(theta))
    y <- c(radius_inner * sin(rev(theta)), radius_outer * sin(theta))
    polygon(x, y, col = diversity_plot[family == f & genus == g, 
                                       genus_col], 
            border = "black")
    if(diversity_plot[family == f & genus == g, genus_frac] > 0.005){
      theta_mid <- (outer_start + outer_end) / 2
      r_label <- radius_outer * 1.15  # push label just outside the outer edge
      x_label <- r_label * cos(theta_mid)
      y_label <- r_label * sin(theta_mid)
      text(x_label, y_label, labels = g, cex = 0.75,
           adj = ifelse(cos(theta_mid) > 0, 0, 1))
      lines(x = c(x_label-0.01, cos(theta_mid)), 
            y = c(y_label-0.01, sin(theta_mid)), lwd = 0.8)
    }
  
    outer_start <- outer_end
  }
  

  start <- end
}

text(0.2,1.2, cex = 1.5,
     "non-Klebsiella pneumoniae; Percent Identity > 95")

# family legend
plot.new()
with(unique(diversity_plot[family!="Enterobacteriaceae", .(family, family_col)]),
     legend(0,0.85, legend = family,
            xpd = TRUE, title = "Non-Enterobacteriaceae families",
            fill = family_col, cex = 0.6, bty = "n")
)

# Species ANI distributions


# Ridge plot --------------------------------------------------------------
ridge_plot <- kpneu_uhgg_hits[species != 'Klebsiella pneumoniae' & 
                                Percent_Identity > 95 & family == "Enterobacteriaceae"]

kleb_ridge_plot <- ridge_plot[grep("Klebsiella", genus)]
kleb_ridge_plot[, n:= .N, by = species]
kleb_ridge_plot <- kleb_ridge_plot[n>3]
kleb_ridge_plot[, n:= .N, by = species]
kleb_ridge_plot[, kleb_plot_species := paste0(species," (n = ", n, ")")]

# library
library(ggridges)
library(ggplot2)

# Diamonds dataset is provided by R natively
#head(diamonds)

# basic example
ggplot(kleb_ridge_plot, aes(x = Percent_Identity, y = kleb_plot_species, fill = kleb_plot_species)) +
  geom_density_ridges(alpha = 0.5,  # 1. Transparency
                      scale = 1.5  # 3. Increase offset between ridges
                      ) +
  scale_x_reverse(limits = c(100, 95)) +  # 2. Reverse x-axis: from 10000 to 0
  theme_ridges() + 
  theme(legend.position = "none") +
  labs(y = "Species", x = "Percent identity", title = "Klebsiella")



# Raoultella
raul_ridge_plot <- ridge_plot[grep("Raoultella", genus)]
raul_ridge_plot[, n:= .N, by = species]
# raul_ridge_plot <- raul_ridge_plot[n>10]
raul_ridge_plot[, n:= .N, by = species]
raul_ridge_plot[, kleb_plot_species := paste0(species," (n = ", n, ")")]

# Diamonds dataset is provided by R natively
#head(diamonds)

# basic example
ggplot(raul_ridge_plot, aes(x = Percent_Identity, y = kleb_plot_species, fill = kleb_plot_species)) +
  geom_density_ridges(alpha = 0.5,  # 1. Transparency
                      scale = 1.5  # 3. Increase offset between ridges
  ) +
  scale_x_reverse(limits = c(100, 95)) +  # 2. Reverse x-axis: from 10000 to 0
  theme_ridges() + 
  theme(legend.position = "none") +
  labs(y = "Species", x = "Percent identity", title = "Raoultella")

# Enterobacter
enterobacter_ridge_plot <- ridge_plot[genus=='Enterobacter']
enterobacter_ridge_plot[, n:= .N, by = species]
# enterobacter_ridge_plot <- enterobacter_ridge_plot[n>10]
enterobacter_ridge_plot[, n:= .N, by = species]
enterobacter_ridge_plot[, kleb_plot_species := paste0(species," (n = ", n, ")")]

# Diamonds dataset is provided by R natively
#head(diamonds)

# basic example
ggplot(enterobacter_ridge_plot, aes(x = Percent_Identity, y = kleb_plot_species, fill = kleb_plot_species)) +
  geom_density_ridges(alpha = 0.5,  # 1. Transparency
                      scale = 1.5  # 3. Increase offset between ridges
  ) +
  scale_x_reverse(limits = c(100, 95)) +  # 2. Reverse x-axis: from 10000 to 0
  theme_ridges() + 
  theme(legend.position = "none") +
  labs(y = "Species", x = "Percent identity", title = "Enterobacter")

