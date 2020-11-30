library("omicade4")
library("gliomaData")
library("ggplot2")

A <- readRDS("data.RDS")
A <- lapply(A, t)
out <- mcia(A)
saveRDS(out, "mcia.RDS")
data_plot <- cbind(out$mcoa$SynVar, m)
data_plot %>% 
  ggplot() +
  geom_point(aes(SynVar1, SynVar2, col = Location, shape = Location)) +
  labs(title = "MCIA on the pouchitis dataset")
ggsave("mcia.png")

p1 <- data_plot %>% 
  ggplot() +
  geom_roc(aes(d = Location, m = SynVar1), n.cuts = 0) +
  style_roc()
p1 +
  annotate("text", x = .5, y = 0.5, label = paste("AUC =", round(calc_auc(p1)[, "AUC"], 3))) +
  labs(title = "AUC for MCIA", subtitle = "Classification by location")
ggsave("mcia_AUC.png")
