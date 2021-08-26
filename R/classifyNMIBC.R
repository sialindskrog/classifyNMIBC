#' Transcriptomic UROMOL2021 classes of non-muscle-invasive bladder cancer (NMIBC) samples
#'
#' This package performs Pearson nearest-centroid classification according to the transcriptomic classes of NMIBC based on gene expression profiles.
#'
#' @param x dataframe with unique genes in rows and samples to be classified in columns (or single named vector of gene expression values).
#' RNA-seq data needs to be log-transformed and micro-array data should be normalized. Gene names may be supplied as Ensembl gene IDs or HUGO gene symbols.
#' @param minCor a numeric value specifying a minimal threshold for best Pearson's correlation between a sample's gene expression profile and centroids profiles.
#' A sample showing no correlation above this threshold will remain unclassifed and prediction results will be set to NA. Default minCor value is 0.2.
#' @param gene_id a character value specifying the type of gene identifiers used for the names/rownames of 'x',
#' ensembl_gene_ID for Ensembl gene IDs or hgnc_symbol for HUGO gene symbols.
#'
#' @return A dataframe with classification results. The NMIBC_Class column returns the predicted class labels of the samples.
#' The cor_pval column returns the p-value associated with the Pearson's correlation between the sample and the nearest centroid.
#' The separationLevel gives a measure (ranging from 0 to 1) of how a sample is representative of its consensus class,
#' with 0 meaning the sample is too close to the other classes to be confidently assigned to one class label, and 1 meaning the sample is very representative of its class.
#' This separationLevel is measured as follows: (correlation to nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid correlation.
#' The remaing four columns return the Pearson's correlation values for each sample and each centroid.
#' NMIBC_Class predictions are set to NA if the minCor condition is not verified.
#'
#' @export
#'
#' @examples
#' data(test_data)
#' classifyNMIBC(test_data)
#'
#' @note This classifier was built similarly to the published tool for the consensus classes of MIBC:
#' Kamoun, A et. al. A Consensus Molecular Classification of Muscle-invasive Bladder Cancer. Eur Uro (2019), doi: https://doi.org/10.1016/j.eururo.2019.09.006
#'
#' Alternatively, you can classifiy samples through our online web application: http://nmibc-class.dk



classifyNMIBC <- function(x, minCor = .2, gene_id = c("ensembl_gene_ID", "hgnc_symbol")[1]){

  data(centroids)

  classes <- c("Class_1", "Class_2a", "Class_2b", "Class_3")

  gkeep <- intersect(centroids[, gene_id], rownames(x))
  if (length(gkeep) == 0) stop("Empty intersection between profiled genes and the genes used for classification.\n Make sure that gene names correspond to the type of identifiers specified by the gene_id argument")
  if (length(gkeep) < 0.6 * nrow(centroids)) warning("Input gene expression profiles include less than 60% of the genes used for classification. Results may not be relevant")

  cor.dat <- as.data.frame(cor(x[gkeep, ], centroids[match(gkeep, centroids[, gene_id]), classes], use = "complete.obs"), row.names = colnames(x))

  # Best correlated centroid
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(y){classes[which.max(y)]})
  cor.dat$corToNearest <- apply(cor.dat[, classes], 1, max)
  cor.dat$cor_pval <- sapply(colnames(x), function(smp){
    cor.test(x[gkeep, smp], centroids[match(gkeep, centroids[, gene_id]), cor.dat[smp, "nearestCentroid"]])$p.value
  })

  # Separation level metrics
  cor.dat$deltaSecondNearest <- apply(cor.dat$corToNearest - cor.dat[, classes], 1, function(x){sort(x)[2]})
  cor.dat$deltaMed <- apply(cor.dat$corToNearest - cor.dat[, classes], 1, median)
  cor.dat$separationLevel <- cor.dat$deltaSecondNearest/cor.dat$deltaMed

  cor.dat$NMIBC_class <- cor.dat$nearestCentroid

  # Set to NA if best correlation < minCor
  try(cor.dat[which(cor.dat$corToNearest < minCor), "NMIBC_Class"] <-  NA)
  try(cor.dat[which(cor.dat$corToNearest < minCor), "separationLevel"] <-  NA)

  cor.dat <- cor.dat[, c("NMIBC_class" , "cor_pval", "separationLevel", classes)]
  return(cor.dat)
}

