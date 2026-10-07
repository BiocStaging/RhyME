#' PQL Cosinor Estimation and Rhythmicity Testing
#'
#' Fits beta-regression cosinor models with optional subject random
#' intercepts using penalized quasi-likelihood (PQL), an optional empirical
#' Bayes (EB) prior on log precision, and a null-model score test for
#' rhythmicity.
#'
#' @details Input and output are `SummarizedExperiment` objects. The fixed
#' effects are an intercept, optional adjustment covariates, and sine/cosine
#' terms. Omitting `subject_col` fits a fixed-effects-only model for
#' independent samples.
#'
#' @seealso [fitCosinor()], [plotCosinorFit()]
#'
#' @keywords internal
#' @useDynLib RhyME, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom SummarizedExperiment assay assayNames colData rowData rowData<-
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom ggplot2 .data
"_PACKAGE"
