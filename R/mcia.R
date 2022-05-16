library("omicade4")
library("gliomaData")
library("ggplot2")

A <- readRDS("data.RDS")
A <- lapply(A, t)
out <- mcia(A)
saveRDS(out, "mcia.RDS")
meta <- read.delim("data/from_zip/metadata.txt", row.names = 1)
m <- meta[rownames(out$mcoa$SynVar), ]
data_plot <- cbind(out$mcoa$SynVar, m)
theme_set(theme_bw())
theme_update(strip.background = element_blank(), 
             panel.grid.minor = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())
data_plot %>% 
  ggplot() +
  geom_point(aes(SynVar1, SynVar2, col = Location, shape = Location)) +
  labs(title = "MCIA in the Morgan dataset")
ggsave("Figures/mcia.png")
ggsave("~/Documents/projects/thesis/images/morgan-mcia.png", width = 170,
       units = "mm", dpi = 300, bg = "white")

p1 <- data_plot %>% 
  ggplot() +
  geom_roc(aes(d = Location, m = SynVar1), n.cuts = 0) +
  style_roc()
p1 +
  annotate("text", x = .5, y = 0.5, label = paste("AUC =", round(calc_auc(p1)[, "AUC"], 3))) +
  labs(title = "AUC for MCIA", subtitle = "Classification by location")
ggsave("Figures/mcia_AUC.png")
