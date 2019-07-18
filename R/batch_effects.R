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

expr <- expr[, colnames(expr) %in% rownames(m)]
otus <- otus[, colnames(otus) %in% m$MID]

expr <- expr
# otus <- filter_RNAseq(norm_RNAseq(otus))

MR_i <- newMRexperiment(otus)
MR_i <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
otus <- MRcounts(MR_i, norm = TRUE, log = TRUE)

expr.data <- expr-apply(expr, 1, median)
pca.expr <- prcomp(t(expr.data), center = FALSE, scale. = FALSE)
s <- summary(pca.expr)
plot(pca.expr$x, pch = 16, col = m$Location, 
     xlab = paste("PC1", round(s$importance[2, 1]*100, 4), "%"),
     ylab = paste("PC2", round(s$importance[2, 2]*100, 4), "%"))
plot(pca.expr$x, pch = 16, col = m$Outcome, 
     xlab = paste("PC1", round(s$importance[2, 1]*100, 4), "%"),
     ylab = paste("PC2", round(s$importance[2, 2]*100, 4), "%"))
plot(pca.expr$x, pch = 16, col = m$Batch, 
     xlab = paste("PC1", round(s$importance[2, 1]*100, 4), "%"),
     ylab = paste("PC2", round(s$importance[2, 2]*100, 4), "%"))
plot(pca.expr$x, pch = 16, col = m$Gender, 
     xlab = paste("PC1", round(s$importance[2, 1]*100, 4), "%"),
     ylab = paste("PC2", round(s$importance[2, 2]*100, 4), "%"))

otus.data <- otus-apply(otus, 1, median)
pca <- prcomp(t(otus.data), center = TRUE, scale. = TRUE)
s <- summary(pca)
plot(pca$x, pch = 16, col = m$Location, 
     xlab = paste("PC1", round(s$importance[2, 1]*100, 4), "%"),
     ylab = paste("PC2", round(s$importance[2, 2]*100, 4), "%"))
plot(pca$x, pch = 16, col = m$Outcome, 
     xlab = paste("PC1", round(s$importance[2, 1]*100, 4), "%"),
     ylab = paste("PC2", round(s$importance[2, 2]*100, 4), "%"))
plot(pca$x, pch = 16, col = m$Gender, 
     xlab = paste("PC1", round(s$importance[2, 1]*100, 4), "%"),
     ylab = paste("PC2", round(s$importance[2, 2]*100, 4), "%"))
plot(pca$x, pch = 16, col = m$Batch, 
     xlab = paste("PC1", round(s$importance[2, 1]*100, 4), "%"),
     ylab = paste("PC2", round(s$importance[2, 2]*100, 4), "%"))

Demographics <- model_RGCCA(m, c("Gender", "ID_1", "Antibiotics"))
Location <- model_RGCCA(m, c("Location", "ISCORE"))
A <- list(RNAseq =  t(expr), "16S" = t(otus))
shrinkage <- sapply(A, tau.estimate)
C <- matrix(0, ncol = 2, nrow = 2, dimnames = list(names(A), names(A)))
model0 <- subSymm(C, "RNAseq", "16S", 1)
m0 <- sgcca(A, model0, c1 = shrinkage, ncomp = rep(2, length(A)))
m0 <- improve.sgcca(m0, names(A))
saveRDS(m0, "model0.RDS")

meta <- model_RGCCA(m, c("Gender", "ID_1", "Antibiotics", "Location", "ISCORE"))
A$meta <- meta
C <- matrix(0, ncol = length(A), nrow = length(A), dimnames = list(names(A), names(A)))
model1 <- subSymm(C, "RNAseq", "16S", 1)
model1 <- subSymm(model1, "RNAseq", "meta", 1)
model1 <- subSymm(model1, "meta", "16S", 1)
m1 <- sgcca(A, model1, c1 = c(shrinkage, 1), ncomp = rep(2, length(A)))
m1 <- improve.sgcca(m1, names(A))
plot(m1$Y$RNAseq[, 1], m1$Y$`16S`[, 1], pch = 16, col = m$Location)
saveRDS(m1, "model1.RDS")

model1.1 <- subSymm(model1, "RNAseq", "16S", 0)
m1.1 <- sgcca(A, model1.1, c1 = c(shrinkage, 1), ncomp = rep(2, length(A)))
m1.1 <- improve.sgcca(m1.1, names(A))
saveRDS(m1.1, "model1.1.RDS")
plot(m1.1$Y$RNAseq[, 1], m1.1$Y$`16S`[, 1], pch = 16, col = m$ISCORE)


designs1.2 <- weight_design(weights = 11, size = 3)
k <- vapply(designs1.2, correct, logical(1L))
designs1.2 <- designs1.2[k]

Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

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
saveRDS(out, "model1.2_testing.RDS")


columns <- grep("var", colnames(out))
model1.2 <- symm(model1.1, out[which.max(out$AVE_inner), columns])
m1.2 <- sgcca(A, model1.2, c1 = c(shrinkage, 1), ncomp = rep(2, 3))
m1.2 <- improve.sgcca(m1.2, names(A))
saveRDS(m1.2, "model1.2.RDS")

# Models 2 ####

A <- c(A[1:2], list(Location = Location, Demographics = Demographics))
shrinkage2 <- c(shrinkage, 1, 1)
names(shrinkage2) <- names(A)
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model2 <- subSymm(C, "16S", "RNAseq", 1)
model2 <- subSymm(model2, "RNAseq", "Demographics", 1)
model2 <- subSymm(model2, "RNAseq", "Location", 1)
model2 <- subSymm(model2, "16S", "Demographics", 1)
model2 <- subSymm(model2, "16S", "Location", 1)
m2 <- sgcca(A, model2, c1 = shrinkage2, ncomp = rep(2, length(A)))
m2 <- improve.sgcca(m2, names(A))
saveRDS(m2, "model2.RDS")

Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
designs2.1 <- weight_design(weights = 3, size = 4)
k <- vapply(designs2.1, correct, logical(1L))
designs2.1 <- designs2.1[k]
out <- sapply(designs2.1, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "model2.1_testing.RDS")

columns <- grep("var", colnames(out))
model2.1 <- symm(model2, out[which.max(out$AVE_inner), columns])
m2.1 <- sgcca(A, model2.1, c1 = shrinkage2, ncomp = rep(2, 4))
m2.1 <- improve.sgcca(m2.1, names(A))
saveRDS(m2.1, "model2.1.RDS")

designs2.2 <- weight_design(weights = 11, size = 4, 
                            diff0 = which(lower.tri(model2.1) & model2.1 != 0))
out <- sapply(designs2.2, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "model2.2_testing.RDS")

model2.2 <- symm(C, out[which.max(out$AVE_inner), grep("var", colnames(out))])
m2.2 <- sgcca(A, model2.2, c1 = shrinkage2, ncomp = rep(2, 4))
m2.2 <- improve.sgcca(m2.2, names(A))
saveRDS(m2.2, "model2.2.RDS")

if (model2.2["RNAseq", "16S"] == 0)  {
  out <- readRDS("model2.1_testing.RDS")
  columns <- grep("var", colnames(out))
  sub_out <- out[out$var12 != 0, ]
  
  model2.3 <- symm(model2, sub_out[which.max(sub_out$AVE_inner), columns])
  designs2.3 <- weight_design(weights = 11, size = 4, 
                              diff0 = which(lower.tri(model2.3) & model2.3 != 0))
  out <- sapply(designs2.3, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
  out <- simplify2array(out[lengths(out) != 1])
  out <- as.data.frame(t(out))
  saveRDS(out, "model2.3_testing.RDS")
  
  
  model2.3 <- symm(C, out[which.max(out$AVE_inner), grep("var", colnames(out))])
  m2.3 <- sgcca(A, model2.3, c1 = shrinkage2, ncomp = rep(2, 4))
  m2.3 <- improve.sgcca(m2.3, names(A))
  saveRDS(m2.3, "model2.3.RDS")
}

l <- list.files(pattern = "model.*[0-9]_edge.RDS")
sapply(l, function(x){readRDS(x)$C})
ls <- sapply(l, function(x){readRDS(x)$AVE$AVE_inner})
lsa <- gsub("_edge\\.RDS", "", l)
pos <- order(as.numeric(gsub("model", "", lsa)))
ls[, pos]

l <- list.files(pattern = "model.*[0-9].RDS")
sapply(l, function(x){readRDS(x)$C})
ls <- sapply(l, function(x){readRDS(x)$AVE$AVE_inner})
lsa <- gsub("\\.RDS", "", l)
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

ggplot(s3) +
  geom_point(aes(RNAseq, `16S`, col = Outcome)) +
  facet_wrap(~Model, nrow = 2, scales = "free") +
  theme_bw()
ggplot(s3) +
  geom_point(aes(RNAseq, `16S`, col = ISCORE)) +
  facet_wrap(~Model, nrow = 2, scales = "free") +
  theme_bw()
ggplot(s3) +
  geom_point(aes(RNAseq, `16S`, col = Gender)) +
  facet_wrap(~Model, nrow = 2, scales = "free") +
  theme_bw()
ggplot(s3) +
  geom_point(aes(RNAseq, `16S`, col = Location)) +
  facet_wrap(~Model, nrow = 2, scales = "free") +
  theme_bw()
