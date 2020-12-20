library(data.table)
library(tidyverse)
library(plotrix)

dt <- fread('~/Downloads/roc.csv', header = FALSE) %>%
  setnames(., c("z", paste("bucket", 1:10, sep = "_"), "depth", "read_length", "window_size", "min_reads", "max_reads", "minimal_score", "max_HF", "ref_ancestry")) %>%
  pivot_longer(!c("z", "depth", "read_length", "window_size", "min_reads", "max_reads", "minimal_score", "max_HF", "ref_ancestry"), names_to = "bucket", values_to = "xy") %>%
  as.data.table() %>%
  .[, xy := gsub("\\(|\\)", "", xy)] %>%
  .[, c("x", "y") := tstrsplit(xy, ",", fixed = TRUE)] %>%
  .[, x := as.numeric(x)] %>%
  .[, y := as.numeric(y)] %>%
  .[, xy := NULL]

# ggplot(data = dt, aes(x = x, y = y, color = factor(depth))) +
#   geom_point() +
#   theme_bw() +
#   geom_abline(slope = 1, lty = "dashed") +
#   xlab("BFPR") +
#   ylab("BTPR")

dt_summary <- dt %>%
  group_by(z, depth, read_length, ref_ancestry) %>%
  #group_by(z, depth, read_length, window_size, min_reads, max_reads, minimal_score, max_HF, ref_ancestry) %>%
  summarize(., mean_x = mean(x), mean_y = mean(y), se = std.error(y))
  
ggplot(data = dt_summary, aes(x = mean_x, y = mean_y, ymin = mean_y - se, ymax = mean_y + se, fill = factor(depth), color = factor(depth))) +
  geom_ribbon(alpha = 0.25, color = NA) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_abline(slope = 1, lty = "dashed") +
  xlab("BFPR") +
  ylab("BTPR") +
  scale_color_brewer(name = "Depth of \nCoverage", palette = "Set2") +
  scale_fill_brewer(name = "Depth of \nCoverage", palette = "Set2") +
  facet_grid(ref_ancestry ~ read_length, labeller = label_both)


dt_summary$ref_ancestry <- factor(dt_summary$ref_ancestry, levels = c("EUR", "EAS"))
ggplot(data = dt_summary, aes(x = mean_x, y = mean_y, ymin = mean_y - se, ymax = mean_y + se, fill = factor(depth), color = factor(depth), lty = ref_ancestry)) +
  geom_ribbon(alpha = 0.25, color = NA) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_abline(slope = 1, lty = "dashed") +
  xlab("BFPR") +
  ylab("BTPR") +
  scale_color_brewer(name = "Depth of \nCoverage", palette = "Set2") +
  scale_fill_brewer(name = "Depth of \nCoverage", palette = "Set2") +
  facet_grid(. ~ read_length, labeller = label_both)
  
  
