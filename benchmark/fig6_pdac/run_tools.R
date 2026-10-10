# Fig-6 benchmark: run Scissor / scAB / DEGAS on replicate r of the PDAC Xenium subsample.
# Usage: Rscript run_tools.R <method> <r>
# Each tool uses its own published defaults and its own high-risk selection rule.
# Output: fig6_pdac/results/<method>_<r>.csv with columns cell, high (0/1), score.
suppressMessages({library(Seurat); library(Matrix)})
options(Seurat.object.assay.version = "v3")   # Scissor/scAB read assay@data (v3 slots)
args <- commandArgs(TRUE); method <- args[1]; r <- as.integer(args[2])
set.seed(r)
IN <- "fig6_pdac/inputs"; OUT <- "fig6_pdac/results"; dir.create(OUT, showWarnings = FALSE)
counts <- as.matrix(read.csv(sprintf("%s/sc_counts_%d.csv", IN, r), row.names = 1, check.names = FALSE))
meta   <- read.csv(sprintf("%s/sc_meta_%d.csv", IN, r), row.names = 1)
bulk   <- as.matrix(read.csv(sprintf("%s/bulk_expr.csv", IN), row.names = 1, check.names = FALSE))
pheno  <- read.csv(sprintf("%s/pheno.csv", IN), row.names = 1)
pheno$time <- pmax(pheno$time, 1)             # Cox/guanrank need positive times
t0 <- Sys.time()

if (method == "scissor") {
  library(Scissor)
  # Scissor::Seurat_preprocessing defaults, without the t-SNE/UMAP it does not use
  sc <- CreateSeuratObject(counts = counts, min.cells = 400, min.features = 0)
  sc <- NormalizeData(sc, verbose = FALSE); sc <- FindVariableFeatures(sc, verbose = FALSE)
  sc <- ScaleData(sc, verbose = FALSE); sc <- RunPCA(sc, features = VariableFeatures(sc), verbose = FALSE)
  sc <- FindNeighbors(sc, dims = 1:10, verbose = FALSE)
  ph <- as.matrix(pheno[colnames(bulk), c("time", "status")])
  res <- Scissor(bulk, sc, ph, family = "cox", alpha = NULL, cutoff = 0.2,
                 Save_file = sprintf("%s/scissor_tmp_%d.RData", OUT, r))
  score <- setNames(rep(0, ncol(counts)), colnames(counts))
  score[res$Scissor_pos] <- 1; score[res$Scissor_neg] <- -1
  high <- as.integer(score == 1)                # Scissor+ = associated with worse survival
  file.remove(sprintf("%s/scissor_tmp_%d.RData", OUT, r))
  cat("alpha", res$para$alpha, "\n")

} else if (method == "scab") {
  library(scAB)
  sc <- CreateSeuratObject(counts = counts)
  sc <- NormalizeData(sc, verbose = FALSE); sc <- FindVariableFeatures(sc, nfeatures = 3000, verbose = FALSE)
  sc <- ScaleData(sc, verbose = FALSE); sc <- RunPCA(sc, features = VariableFeatures(sc), verbose = FALSE)
  sc <- FindNeighbors(sc, dims = 1:40, verbose = FALSE)            # scAB::run_seurat defaults
  ph <- as.matrix(pheno[colnames(bulk), c("time", "status")])
  obj <- create_scAB(sc, bulk, ph, method = "survival")
  K <- select_K(obj); cat("K", K, "\n")
  fit <- scAB(Object = obj, K = K)
  sc <- findSubset(sc, scAB_Object = fit, tred = 2)
  high <- as.integer(sc$scAB_select == "scAB+ cells")
  # continuous score: max standardised loading over the modules scAB selected cells from
  H <- fit$H; Z <- (H - rowMeans(H)) / apply(H, 1, sd)
  score <- apply(Z, 2, max)

} else if (method == "degas") {
  library(DEGAS)
  initDEGAS(); set_seed_term(r)
  DEGAS.pyloc <<- Sys.getenv("DEGAS_PYTHON", "python")   # python with TensorFlow for DEGAS
  tmp <- sprintf("%s/degas_tmp_%d/", OUT, r); dir.create(tmp, showWarnings = FALSE)
  scX  <- preprocessCounts(counts)                       # cells x genes
  patX <- preprocessCounts(2^bulk - 1)                   # bulk is log2; DEGAS expects counts scale
  patLab <- as.matrix(pheno[colnames(bulk), c("time", "status")])
  m <- runCCMTLBag(scX, NULL, patX, patLab, tmp, "BlankCox", "DenseNet", 3, 5)
  haz <- as.numeric(predClassBag(m, scX, "pat"))
  # DEGAS PDAC spatial example: spatial kNN smoothing (k=6), min-max scale, > 0.5
  haz <- knnSmooth(haz, as.matrix(meta[colnames(counts), c("x_centroid", "y_centroid")]), k = 6)
  score <- scaleFunc(haz); high <- as.integer(score > 0.5)
  unlink(tmp, recursive = TRUE)
}

# every score vector above is already in colnames(counts) order
out <- data.frame(cell = colnames(counts), high = high, score = as.numeric(score))
write.csv(out, sprintf("%s/%s_%d.csv", OUT, method, r), row.names = FALSE)
cat(method, r, "high", sum(high), "of", length(high), "time",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
