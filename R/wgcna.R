library("WGCNA")
wgcna <- cbind(A[[1]], A[[2]])
type <- c(rep("RNA", ncol(A[[1]])), rep("DNA", ncol(A[[2]])))
sft <- pickSoftThreshold(wgcna)
enableWGCNAThreads()
net <- blockwiseModules(wgcna,
  power = 6,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "wgcna",
  verbose = 3
)


tab <- table(net$colors, type)

keep <- apply(tab[-1, ], 1, function(x){all(x != 0)})
