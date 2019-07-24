library("pROC")
X <- read.table("16S_models.csv",
  sep = ",", dec = ".", header = TRUE,
  row.names = 1, check.names = FALSE
)
####  https://www.r-bloggers.com/simple-roc-plots-with-ggplot2-part-1/
#' Produces x and y co-ordinates for ROC curve plot
#'
#' Calculates the confusion matrix values
#' @param grp Factor of labels classifying subject status
#' @param pred Numerical values of each observation
#' @return List with 2 components:
#'         roc = data.frame with x and y co-ordinates of plot
#'         stats = data.frame containing: area under ROC curve, p value, upper
#'         and lower 95% confidence interval
rocdata <- function(grp, pred) {
  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  }

  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }

  cut <- unique(pred)
  tp <- sapply(cut, function(x) sum(pred > x & grp == levels(grp)[2], na.rm = TRUE))
  fn <- sapply(cut, function(x) sum(pred < x & grp == levels(grp)[2], na.rm = TRUE))
  fp <- sapply(cut, function(x) sum(pred > x & grp == levels(grp)[1], na.rm = TRUE))
  tn <- sapply(cut, function(x) sum(pred < x & grp == levels(grp)[1], na.rm = TRUE))
  predicted_positive <- tp + fp
  predicted_negative <- fn + tn
  condition_positive <- tp + fn
  condition_negative <- fp + tn
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  ppv <- tp / (tp + fp)
  roc <- data.frame(FPR = fpr, TPR = tpr, PPV = ppv)
  roc <- roc[order(roc$FPR, roc$TPR), ]

  i <- 2:nrow(roc)
  auc <- (roc$FPR[i] - roc$FPR[i - 1]) %*% (roc$TPR[i] + roc$TPR[i - 1]) / 2

  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc^2) / (1 + auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) - 1) * (q1 - auc^2)) + (
    (length(neg) - 1) * (q2 - auc^2))) / (length(pos) * length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)

  se.auc.null <- sqrt((1 + length(pos) + length(neg)) / (
    12 * length(pos) * length(neg)))
  z <- (auc - 0.5) / se.auc.null
  p <- 2 * pnorm(-abs(z))

  stats <- data.frame(
    auc = auc,
    p.value = p,
    ci.upper = ci.upper,
    ci.lower = ci.lower
  )

  return(list(roc = roc, stats = stats))
}

####################################### 16S ##################################
X <- read.table("16S_models.csv",
  sep = ",", dec = ".", header = TRUE,
  row.names = 1, check.names = FALSE
)

flag <- rep(0, nrow(X))
flag[which(X$Location == "Pouch")] <- 1
var <- X[, 6:13]

pdf("easyROC_16S.pdf")

for (i in 1:dim(var)[2]) {
  roc <- rocdata(flag, var[, i])

  cat(colnames(var)[i], " roc asociated p-value ", roc$stats[, 2], "\n")
  print(auc(flag, var[, i]))

  plot.roc(flag, var[, i], 
           main = paste(colnames(var)[i], "\npvalue", roc$stats[, 2]),
           print.auc = TRUE, ci = TRUE)
}

dev.off()
########################################################################
####################################### rnaseq  ##################################
X <- read.table("RNASeq_models.csv",
  sep = ",", dec = ".", header = TRUE,
  row.names = 1, check.names = FALSE
)

flag <- matrix(0, length(X$Location))
flag[which(X$Location == "Pouch")] <- 1
var <- X[, 6:13]

pdf("easyROC_RNAseq.pdf")

for (i in 1:dim(var)[2]) {
  roc <- rocdata(flag, var[, i])

  auc(flag, var[, i])

  cat(colnames(var)[i], " roc asociated p-value ", roc$stats[, 2], "\n")
  print(auc(flag, var[, i]))

  plot.roc(flag, var[, i], main = colnames(var)[i], print.auc = TRUE, ci = T)
}
dev.off()
