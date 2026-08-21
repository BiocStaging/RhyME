#' @keywords internal
"_PACKAGE"

#' @importFrom stats sd pchisq p.adjust plogis rbeta rnorm runif var qlogis nlminb optimHess
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom SummarizedExperiment assay assayNames colData rowData rowData<-
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom TMB MakeADFun
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_wrap labs theme_bw .data
#' @useDynLib RhyME, .registration = TRUE
NULL
