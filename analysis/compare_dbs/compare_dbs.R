# compare across databases and species

all_spp_db_dat <- fread(paste0(outdir_dat, "/all_spp_db_dat.csv"))
all_spp_db_dat[, p_a :=fcase(is.na(genus), 0,
                             default = 1)]

#hits
hit_stats <- unique(all_spp_db_dat[,.(spp,db,gene_family,p_a)])

hit_stats[, (n=.N), by = c("spp", "db", "p_a")]



# How many blast hits? ----------------------------------------------------

# Sample data
data <- matrix(c(3, 2, 4,
                 4, 1, 3,
                 2, 5, 1),
               nrow = 3, byrow = TRUE)

# Add row and column names (optional)
rownames(data) <- c("A", "B", "C")
colnames(data) <- c("X", "Y", "Z")

# Create stacked bar chart
barplot(data, horiz = TRUE,
        beside = FALSE,         # FALSE = stacked bars
        col = c("skyblue", "orange", "forestgreen"),
        legend.text = TRUE,
        main = "Stacked Bar Chart (Base R)")
