list.of.packages <- c("data.table", "tidyverse","cowplot","class","plotly","ggpointdensity")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)
library(tidyverse)
library(cowplot)
library(class)
library(plotly)
library(ggpointdensity)

setwd("~/Johns Hopkins/mccoy_lab - Documents/projects/distinguishing_meiotic_and_mitotic_origins_of_trisomies_v2/LASER/")

ref <- fread("laser.RefPC.coord")

igsr <- fread("igsr_samples.txt")[, c("Sample name", "Superpopulation code")] %>%
  setnames(., c("indivID", "superpop"))

ref <- merge(ref, igsr, by = "indivID") %>%
  .[, popID := NULL]


# simulate first-generation admixed individuals

admix <- function(pop1, pop2, ref_pca, index) {
  admixed_sample <- rbind(ref_pca[superpop == pop1][sample(.N, 1)],
        ref_pca[superpop == pop2][sample(.N, 1)]) %>%
    .[, 2:5] %>%
    colMeans()
  return(data.table(indivID = paste("sample", index, sep = "_"),
                    PC1 = admixed_sample[1],
                    PC2 = admixed_sample[2],
                    PC3 = admixed_sample[3],
                    PC4 = admixed_sample[4],
                    superpop = paste(pop1, pop2, "admixed", sep = "_")))
}

set.seed(1)
ref <- rbind(ref,
             rbindlist(lapply(1:200, function(x) admix("EUR", "EAS", ref, x))),
             rbindlist(lapply(201:400, function(x) admix("EUR", "AFR", ref, x))),
             rbindlist(lapply(401:600, function(x) admix("EUR", "SAS", ref, x))),
             rbindlist(lapply(601:800, function(x) admix("EAS", "SAS", ref, x))))
      
# load real data

seq <- fread("laser.SeqPC.coord")
seq <- seq[complete.cases(seq)]

seq[, inferred_superpop := knn(ref[, grepl("PC", colnames(ref)), with = FALSE],
                               test = seq[, grepl("PC", colnames(seq)), with = FALSE], 
                               cl = ref$superpop, k = 10, l = 0, prob = TRUE)]

table(seq$inferred_superpop)

a <- ggplot(data = ref, aes(x = PC1, y = PC2, color = superpop)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_discrete(name = "1KGP Superpop.") +
  xlim(-300, 250) +
  ylim(-300, 250) +
  NULL

b <- ggplot() +
  geom_pointdensity(data = seq, aes(x = PC1, y = PC2), size = 0.1) +
  #geom_label_repel(data = seq, aes(x = PC1, y = PC2, label = indivID), size = 2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  scale_color_viridis_c(name = "n embryos", limits = range(0, 300)) +
  xlim(-300, 250) +
  ylim(-300, 250) +
  NULL

c <- ggplot(data = ref, aes(x = PC2, y = PC3, color = superpop)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_discrete(name = "1KGP Superpop.") +
  xlim(-300, 250) +
  ylim(-300, 250) +
  NULL

d <- ggplot() +
  geom_pointdensity(data = seq, aes(x = PC2, y = PC3), size = 0.1) +
  #geom_label_repel(data = seq, aes(x = PC1, y = PC2, label = indivID), size = 2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  scale_color_viridis_c(name = "n embryos", limits = range(0, 300)) +
  xlim(-300, 250) +
  ylim(-300, 250) +
  NULL

#plot_grid(a, b, c, d, labels = c("A.", "B.", "C.", "D."), rel_widths = c(1, 0.9), nrow = 2)

e <- ggplot() +
  geom_point(data = seq, aes(x = PC1, y = PC2, color = inferred_superpop), size = 0.1) +
  #geom_label_repel(data = seq, aes(x = PC1, y = PC2, label = indivID), size = 2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  xlim(-300, 250) +
  ylim(-300, 250) +
  NULL

f <- ggplot() +
  geom_point(data = seq, aes(x = PC2, y = PC3, color = inferred_superpop), size = 0.1) +
  #geom_label_repel(data = seq, aes(x = PC1, y = PC2, label = indivID), size = 2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  xlim(-300, 250) +
  ylim(-300, 250) +
  NULL

plot_grid(a, e, c, f, labels = c("A.", "B.", "C.", "D."), rel_widths = c(1, 1), nrow = 2)

table(seq$inferred_superpop) %>% as.data.table()

# plot_ly(x = seq$PC1, y = seq$PC2, z = seq$PC3, type = "scatter3d", mode = "markers", color = seq$inferred_superpop, size = 2)
# plot_ly(x = ref$PC1, y = ref$PC2, z = ref$PC3, type = "scatter3d", mode = "markers", color = ref$superpop, size = 2)

fwrite(seq[, -1], file = "zouves_inferred_ancestry.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
