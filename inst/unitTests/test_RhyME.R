test_RhyME <- function() {
    example <- readRDS(system.file("extdata", "example_data.rds",
        package = "RhyME"))
    fit <- fitCosinor(example$se[1:5, ], subject_col = "subject",
        prior = "none", BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    result <- rowData(fit)

    checkEquals(dim(result), c(5L, 13L))
    checkTrue(all(is.finite(result$b0)))
    checkTrue(all(is.finite(result$pvalue)))
    checkTrue(all(result$pvalue >= 0 & result$pvalue <= 1))
    checkTrue(all(result$converged))
    checkTrue(inherits(plotCosinorFit(fit), "ggplot"))
}
