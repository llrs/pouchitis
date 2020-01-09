library("integration")
library("experDesign")
library("RGCCA")
library("metagenomeSeq")

metadata <- read.delim("data/input/gene expression metadata gold.txt")

index <- use_index(metadata$Batch)
pheno <- metadata[, c("ID_1", "Gender", "Location", "ISCORE", "Antibiotics", "Outcome")]
batch_effect <- evaluate_index(index, pheno)
orig_effect <- evaluate_orig(pheno)

otus <- read.delim("data/from_zip/otu_table.txt", skip = 1)
expr <- read.delim("data/from_zip/gene_expression2.tab", row.names = 1)
metadata <- read.delim("data/from_zip/metadata.txt", row.names = 1)

in_otus <- metadata$MID %in% colnames(otus)
in_expr <- rownames(metadata) %in% colnames(expr)

m <- metadata[in_otus & in_expr, ]
m <- droplevels(m)

expr <- expr[, colnames(expr) %in% rownames(m)]
otus <- otus[, colnames(otus) %in% m$MID]

all0 <- apply(otus, 1, function(x){all(x == 0)})
otus <- otus[!all0, ]

otus <- filter_RNAseq(norm_RNAseq(otus))
# MR_i <- newMRexperiment(otus)
# MR_i <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
# otus <- MRcounts(MR_i, norm = TRUE, log = TRUE)

Demographics <- model_RGCCA(m, c("Gender", "ID_1", "Antibiotics"))
Location <- model_RGCCA(m, c("Location", "ISCORE"))
A <- list(RNAseq =  t(expr), "16S" = t(otus))
shrinkage <- sapply(A, tau.estimate)
C <- matrix(0, ncol = 2, nrow = 2, dimnames = list(names(A), names(A)))
model0 <- subSymm(C, "RNAseq", "16S", 1)
model0 <- subSymm(model0, "RNAseq", "RNAseq", 1)
m0 <- sgcca(A, model0, c1 = shrinkage, ncomp = rep(2, length(A)))
m0 <- improve.sgcca(m0, names(A))
saveRDS(m0, "imodel0_edge.RDS")

meta <- model_RGCCA(m, c("Gender", "ID_1", "Antibiotics", "Location", "ISCORE"))
A1.2 <- A
A1.2$meta <- meta
C <- matrix(0, ncol = length(A1.2), nrow = length(A1.2), dimnames = list(names(A1.2), names(A1.2)))
model1 <- subSymm(C, "RNAseq", "RNAseq", 1)
model1 <- subSymm(model1, "RNAseq", "meta", 1)
model1 <- subSymm(model1, "meta", "16S", 1)
m1 <- sgcca(A1.2, model1, c1 = c(shrinkage, 1), ncomp = rep(2, length(A1.2)))
m1 <- improve.sgcca(m1, names(A1.2))
plot(m1$Y$RNAseq[, 1], m1$Y$`16S`[, 1], pch = 16, col = m$Location)
saveRDS(m1, "imodel1_edge.RDS")

model1.1 <- subSymm(model1, "RNAseq", "16S", 0)
m1.1 <- sgcca(A1.2, model1.1, c1 = c(shrinkage, 1), ncomp = rep(2, length(A1.2)))
m1.1 <- improve.sgcca(m1.1, names(A1.2))
saveRDS(m1.1, "imodel1.1_edge.RDS")
plot(m1.1$Y$RNAseq[, 1], m1.1$Y$`16S`[, 1], pch = 16, col = m$ISCORE)


designs1.2 <- weight_design(weights = 11, size = 3)
k <- vapply(designs1.2, correct, logical(1L))
designs1.2 <- designs1.2[k]

Ab <- lapply(A1.2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

testing <- function(x, ...) {
  try({
    result.sgcca <- RGCCA::sgcca(C = x, 
                                 verbose = FALSE, 
                                 scale = FALSE,
                                 ...)
    analyze(result.sgcca)}, 
    silent = TRUE)
}

out <- sapply(designs1.2, testing, A = Ab, c1 = c(shrinkage, 1), USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "imodel1.2_testing_edge.RDS")


columns <- grep("var", colnames(out))
model1.2 <- symm(model1.1, out[which.max(out$AVE_inner), columns])
m1.2 <- sgcca(A1.2, model1.2, c1 = c(shrinkage, 1), ncomp = rep(2, 3))
m1.2 <- improve.sgcca(m1.2, names(A1.2))
saveRDS(m1.2, "imodel1.2_edge.RDS")

# Models 2 ####

A2.2 <- c(A1.2[1:2], list(Location = Location, Demographics = Demographics))
shrinkage2 <- c(shrinkage, 1, 1)
names(shrinkage2) <- names(A2.2)
C <- matrix(
  0, ncol = length(A2.2), nrow = length(A2.2),
  dimnames = list(names(A2.2), names(A2.2))
)
model2 <- subSymm(C, "RNAseq", "RNAseq", 1)
model2 <- subSymm(model2, "RNAseq", "Demographics", 1)
model2 <- subSymm(model2, "RNAseq", "Location", 1)
model2 <- subSymm(model2, "16S", "Demographics", 1)
model2 <- subSymm(model2, "16S", "Location", 1)
m2 <- sgcca(A2.2, model2, c1 = shrinkage2, ncomp = rep(2, length(A2.2)))
m2 <- improve.sgcca(m2, names(A2.2))
saveRDS(m2, "imodel2_edge.RDS")

Ab <- lapply(A2.2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
designs2.1 <- weight_design(weights = 3, size = 4)
k <- vapply(designs2.1, correct, logical(1L))
designs2.1 <- designs2.1[k]
out <- sapply(designs2.1, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "imodel2.1_testing_edge.RDS")

columns <- grep("var", colnames(out))
model2.1 <- symm(model2, out[which.max(out$AVE_inner), columns])
m2.1 <- sgcca(A2.2, model2.1, c1 = shrinkage2, ncomp = rep(2, 4))
m2.1 <- improve.sgcca(m2.1, names(A2.2))
saveRDS(m2.1, "imodel2.1_edge.RDS")

designs2.2 <- weight_design(weights = 11, size = 4, 
                            diff0 = which(lower.tri(model2.1) & model2.1 != 0))
out <- sapply(designs2.2, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "imodel2.2_testing_edge.RDS")

model2.2 <- symm(C, out[which.max(out$AVE_inner), grep("var", colnames(out))])
m2.2 <- sgcca(A2.2, model2.2, c1 = shrinkage2, ncomp = rep(2, 4))
m2.2 <- improve.sgcca(m2.2, names(A2.2))
saveRDS(m2.2, "imodel2.2_edge.RDS")

if (model2.2["RNAseq", "16S"] == 0)  {
  out <- readRDS("model2.1_testing_edge.RDS")
  columns <- grep("var", colnames(out))
  sub_out <- out[out$var12 != 0, ]
  
  model2.3 <- symm(model2, sub_out[which.max(sub_out$AVE_inner), columns])
  designs2.3 <- weight_design(weights = 11, size = 4, 
                              diff0 = which(lower.tri(model2.3) & model2.3 != 0))
  out <- sapply(designs2.3, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
  out <- simplify2array(out[lengths(out) != 1])
  out <- as.data.frame(t(out))
  saveRDS(out, "imodel2.3_testing_edge.RDS")
  
  
  model2.3 <- symm(C, out[which.max(out$AVE_inner), grep("var", colnames(out))])
  m2.3 <- sgcca(A2.2, model2.3, c1 = shrinkage2, ncomp = rep(2, 4))
  m2.3 <- improve.sgcca(m2.3, names(A2.2))
  saveRDS(m2.3, "imodel2.3_edge.RDS")
}

l <- list.files(pattern = "model.*[0-9]_edge.RDS")
# sapply(l, function(x){readRDS(x)$C})
ls <- sapply(l, function(x){readRDS(x)$AVE$AVE_inner})
lsa <- gsub("_edge\\.RDS", "", l)
pos <- order(as.numeric(gsub("model", "", lsa)))
ls[, pos]

# l <- list.files(pattern = "model.*[0-9].RDS")
# sapply(l, function(x){readRDS(x)$C})
# ls <- sapply(l, function(x){readRDS(x)$AVE$AVE_inner})
# lsa <- gsub("\\.RDS", "", l)
models <- as.numeric(gsub("model", "", lsa))
pos <- order(models)
models[pos]
s <- lapply(seq_along(l), function(x) {
  y <- readRDS(l[x])
  df <- data.frame("RNAseq" = y$Y[[1]][, 1], "16S" = y$Y[[2]][, 1], 
                   "ID" = 1:nrow(y$Y[[1]]), Model = models[x], check.names = FALSE)
  cbind(df, m[, c("Gender", "Location", "Outcome", "ISCORE", "Antibiotics")])
})

s3 <- do.call(rbind, s)

rn <- tidyr::spread(s3[, -2], key = "Model", value = "RNAseq") 
sn <- tidyr::spread(s3[, -1], key = "Model", value = "16S") 

write.csv(rn, "RNASeq_models.csv", row.names = FALSE)
write.csv(sn, "16S_models.csv", row.names = FALSE)

theme_set(theme_bw())
theme_update(strip.background = element_blank())

ggplot(s3) +
  geom_point(aes(RNAseq, `16S`, col = Outcome)) +
  facet_wrap(~Model, nrow = 2, scales = "free")
ggplot(s3) +
  geom_point(aes(RNAseq, `16S`, col = ISCORE)) +
  facet_wrap(~Model, nrow = 2, scales = "free")
ggplot(s3) +
  geom_point(aes(RNAseq, `16S`, col = Gender)) +
  facet_wrap(~Model, nrow = 2, scales = "free")
ggplot(s3) +
  geom_point(aes(RNAseq, `16S`, col = Location)) +
  facet_wrap(~Model, nrow = 2, scales = "free")

## Create some index and permutations to test the variability on the models.



# Boots ####
seed <- 487129
set.seed(seed)
index <- vector("list", length = 1000)
for (i in seq_len(1000)) {
  index[[i]] <- sample(nrow(A2.2[[1]]), replace = TRUE)
}
saveRDS(index, "index_locale.RDS")
index <- readRDS("index_locale.RDS")

base_boot <- function(index, A, C) {
  STAB <- vector("list", length = length(A))
  AVE <- vector("numeric", length = 2)
  names(AVE) <- c("inner", "outer")
  names(STAB) <- names(A)
  
  A <- subsetData(A, index)
  
  try({
    res <- sgcca(A, C, scheme = "centroid", scale = TRUE)
    AVE["inner"] <- res$AVE$AVE_inner
    AVE["outer"] <- res$AVE$AVE_outer
    for (j in seq_along(A)) {
      STAB[[j]] <- res$a[[j]][, 1]
    }
    
  }, silent = TRUE)
  list(AVE = AVE, STAB = STAB)
}
model0 <- readRDS("model0.RDS")
boot0 <- lapply(index, base_boot, A = A[1:2], C = model0$C)
saveRDS(boot0, "boot0.RDS")
model1.2 <- readRDS("model1.2.RDS")
boot1.2 <- lapply(index, base_boot, A = A1.2, C = model1.2$C)
saveRDS(boot1.2, "boot1.2.RDS")
model2.2 <- readRDS("model2.2.RDS")
boot2.2 <- lapply(index, base_boot, A = A2.2, C = model2.2$C)
saveRDS(boot2.2, "boot2.2.RDS")


ave0 <- cbind.data.frame(t(sapply(boot0, function(x)x$AVE)), Model = "0")
ave1.2 <- cbind.data.frame(t(sapply(boot1.2, function(x)x$AVE)), Model = "1.2")
ave2.2 <- cbind.data.frame(t(sapply(boot2.2, function(x)x$AVE)), Model = "2.2")
aves <- rbind(ave0, ave1.2, ave2.2)
aves <- filter(aves, inner != 0 & outer != 0)
m0 <- data.frame(inner = model0$AVE$AVE_inner[1], 
                 outer = model0$AVE$AVE_outer[1],
                 Model = "0")
m1.2 <- data.frame(inner = model1.2$AVE$AVE_inner[1], 
                   outer = model1.2$AVE$AVE_outer[1],
                   Model = "1.2")
m2.2 <- data.frame(inner = model2.2$AVE$AVE_inner[1], 
                   outer = model2.2$AVE$AVE_outer[1],
                   Model = "2.2")
ggplot(aves, aes(inner, outer, col = Model)) +
  geom_point(alpha = 0.25) +
  stat_ellipse() +
  geom_point(data = m0) +
  geom_point(data = m1.2) +
  geom_point(data = m2.2)
