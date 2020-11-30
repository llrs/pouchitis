library("tidyr")
library("dplyr")
library("ggplot2")
library("forcats")
library("UpSetR")
library("grid")
library("org.Hs.eg.db")
library("clusterProfiler")
library("integration")
library("stringr")
library("patchwork")
library("plotROC")

m0 <- readRDS("model0_edge.RDS")
m0i <- readRDS("imodel0_edge.RDS")
m1 <- readRDS("model1_edge.RDS")
m1.1 <- readRDS("model1.1_edge.RDS")
m1.1i <- readRDS("imodel1.1_edge.RDS")
m1.2 <- readRDS("model1.2_edge.RDS")
m1.2i <- readRDS("imodel1.2_edge.RDS")
m2 <- readRDS("model2_edge.RDS")
m2i <- readRDS("imodel2_edge.RDS")
m2.1 <- readRDS("model2.1_edge.RDS")
m2.1i <- readRDS("imodel2.1_edge.RDS")
m2.2 <- readRDS("model2.2_edge.RDS")
m2.2i <- readRDS("imodel2.2_edge.RDS")
m2.3 <- readRDS("model2.2_edge.RDS")

m0i$AVE$AVE_inner - m0$AVE$AVE_inner
m1.1i$AVE$AVE_inner - m1.1$AVE$AVE_inner
m1.2i$AVE$AVE_inner - m1.2$AVE$AVE_inner
m2.1i$AVE$AVE_inner - m2.1$AVE$AVE_inner

# Samples ####

m0GE <- tidyer(m0$Y[[1]], "0", "GE")
m0iGE <- tidyer(m0i$Y[[1]], "0 i", "GE")
m1GE <- tidyer(m1$Y[[1]], "1", "GE")
m1.1GE <- tidyer(m1.1$Y[[1]], "1.1", "GE")
m1.1iGE <- tidyer(m1.1i$Y[[1]], "1.1 i", "GE")
m1.2GE <- tidyer(m1.2$Y[[1]], "1.2", "GE")
m1.2iGE <- tidyer(m1.2i$Y[[1]], "1.2 i", "GE")
m2GE <- tidyer(m2$Y[[1]], "2", "GE")
m2.1GE <- tidyer(m2.1$Y[[1]], "2.1", "GE")
m2.1iGE <- tidyer(m2.1i$Y[[1]], "2.1 i", "GE")
m2.2GE <- tidyer(m2.2$Y[[1]], "2.2", "GE")
m2.2iGE <- tidyer(m2.2i$Y[[1]], "2.2 i", "GE")
m2.3GE <- tidyer(m2.3$Y[[1]], "2.3", "GE")

m0M <- tidyer(m0$Y[[2]], "0", "M")
m0iM <- tidyer(m0i$Y[[2]], "0 i", "M")
m1M <- tidyer(m1$Y[[2]], "1", "M")
m1.1M <- tidyer(m1.1$Y[[2]], "1.1", "M")
m1.1iM <- tidyer(m1.1i$Y[[2]], "1.1 i", "M")
m1.2M <- tidyer(m1.2$Y[[2]], "1.2", "M")
m1.2iM <- tidyer(m1.2i$Y[[2]], "1.2 i", "M")
m2M <- tidyer(m2$Y[[2]], "2", "M")
m2.1M <- tidyer(m2.1$Y[[2]], "2.1", "M")
m2.1iM <- tidyer(m2.1i$Y[[2]], "2.1 i", "M")
m2.2M <- tidyer(m2.2$Y[[2]], "2.2", "M")
m2.2iM <- tidyer(m2.2i$Y[[2]], "2.2 i", "M")
m2.3M <- tidyer(m2.3$Y[[2]], "2.3", "M")

meta <- read.delim("data/from_zip/metadata.txt", row.names = 1)

# This is specific to having M before GE (the columnames wouldn't match oterwise)
merger <- function(x, y) {
  d <- cbind(y, x[, 3, drop = FALSE])
  merge(d, meta, by.x = "Rownames", by.y = "row.names")
}


df <- rbind(
  merger(m0M, m0GE),
  merger(m0iM, m0iGE),
  merger(m1M, m1GE),
  merger(m1.1M, m1.1GE),
  merger(m1.1iM, m1.1iGE),
  merger(m1.2M, m1.2GE),
  merger(m1.2iM, m1.2iGE),
  merger(m2M, m2GE),
  merger(m2.1M, m2.1GE),
  merger(m2.1iM, m2.1iGE),
  merger(m2.2M, m2.2GE),
  merger(m2.2iM, m2.2iGE),
  merger(m2.3M, m2.3GE))
saveRDS(df, "summary_models.RDS")
df <- readRDS("summary_models.RDS")

# Set theme without background on the labels
theme_set(theme_bw())
theme_update(strip.background = element_blank(), 
             panel.grid.minor = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank())

df <- tibble::as_tibble(df)
sex <- df %>% 
  filter(!grepl(" i$", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Gender)) +
  facet_wrap(~Model, scales = "free", nrow = 2) + 
  # guides(col = FALSE) +
  labs(x = "Transcriptome", y = "Microbiome", col = "Sex") 
loc <- df %>% 
  filter(!grepl(" i$", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Location)) +
  facet_wrap(~Model, scales = "free", nrow = 2) + 
  # guides(col = FALSE) +
  labs(x = "Transcriptome", y = "Microbiome") 

p <- sex/loc  + plot_annotation(tag_levels = "A")
ggsave("Figures/Figure_5.png", units = "mm", width = 170, 
       height = 170, dpi = 300)


# ROC ####
p1 <- df %>% 
  filter(!endsWith(Model, " i"),
         Component == "comp1") %>% 
  ggplot() +
  geom_roc(aes(d = Location, m = GE, col = Model), n.cuts =  0) +
  style_roc()
p1 +
  annotate(geom = "text", x = 1, y = seq(0.65, by = -0.05, length.out = 8),
           label = paste("AUC = ", round(calc_auc(p1)[, "AUC"], 3))) +
  labs(title = "AUC for RGCCA models", subtitle = "Classifying by location")
ggsave("Figures/AUC_models.png")

# Weights ####
m0GE <- tidyer(m0$a[[1]], "0", "GE")
m0iGE <- tidyer(m0i$a[[1]], "0 i", "GE")
m1.1GE <- tidyer(m1.1$a[[1]], "1.1", "GE")
m1.1iGE <- tidyer(m1.1i$a[[1]], "1.1 i", "GE")
m1.2GE <- tidyer(m1.2$a[[1]], "1.2", "GE")
m1.2iGE <- tidyer(m1.2i$a[[1]], "1.2 i", "GE")
m2GE <- tidyer(m2$a[[1]], "2", "GE")
m2.1GE <- tidyer(m2.1$a[[1]], "2.1", "GE")
m2.1iGE <- tidyer(m2.1i$a[[1]], "2.1 i", "GE")
m2.2GE <- tidyer(m2.2$a[[1]], "2.2", "GE")
m2.2iGE <- tidyer(m2.2i$a[[1]], "2.2 i", "GE")
m2.3GE <- tidyer(m2.3$a[[1]], "2.3", "GE")

m0M <- tidyer(m0$a[[2]], "0", "M")
m0iM <- tidyer(m0i$a[[2]], "0 i", "M")
m1.1M <- tidyer(m1.1$a[[2]], "1.1", "M")
m1.1iM <- tidyer(m1.1i$a[[2]], "1.1 i", "M")
m1.2M <- tidyer(m1.2$a[[2]], "1.2", "M")
m1.2iM <- tidyer(m1.2i$a[[2]], "1.2 i", "M")
m2M <- tidyer(m2$a[[2]], "2", "M")
m2.1M <- tidyer(m2.1$a[[2]], "2.1", "M")
m2.1iM <- tidyer(m2.1i$a[[2]], "2.1 i", "M")
m2.2M <- tidyer(m2.2$a[[2]], "2.2", "M")
m2.2iM <- tidyer(m2.2i$a[[2]], "2.2 i", "M")
m2.3M <- tidyer(m2.3$a[[2]], "2.3", "M")



dfGE <- rbind(m0GE, m0iGE, m1.1GE, m1.1iGE, m1.2GE, m1.2iGE, m2GE, m2.1GE, m2.1iGE, m2.2GE, m2.2iGE, m2.3GE)
dfM <- rbind(m0M, m0iM, m1.1M, m1.1iM, m1.2M, m1.2iM, m2M, m2.1M, m2.1iM, m2.2M, m2.2iM, m2.3M)

keepGE <- dfGE %>% 
  # filter(!grepl(" i", Model)) %>%
  filter(Component == "V1" & GE != 0) %>%
  mutate(Presence = if_else(GE != 0, 1, 0)) %>% 
  dplyr::select(-Component, -GE, Rownames) %>% 
  group_by(Rownames) %>% 
  spread(Model, Presence) %>% 
  as.data.frame()
rownames(keepGE) <- keepGE$Rownames
keepGE <- keepGE[, -grep("Rownames", colnames(keepGE))]
keepGE[is.na(keepGE)] <- 0

keepM <- dfM %>% 
  # filter(!grepl(" i", Model)) %>%
  filter(Component == "comp1" & M != 0) %>% 
  mutate(Presence = if_else(M != 0, 1, 0)) %>% 
  dplyr::select(-Component, -M, Rownames) %>% 
  group_by(Rownames) %>% 
  spread(Model, Presence) %>% 
  ungroup() %>% 
  as.data.frame()
rownames(keepM) <- keepM$Rownames
keepM <- keepM[, -grep("Rownames", colnames(keepM))]
keepM[is.na(keepM)] <- 0



text_sizes <- c(1.3, 1.3, 1, 1, 1.5, 1.5)
upset(keepGE, order.by = "freq", nsets = 6, 
      sets = rev(colnames(keepGE)), keep.order = TRUE,
      line.size = NA, text.scale = text_sizes, scale.sets = "identity")
