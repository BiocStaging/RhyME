library(RhyME)

example <- readRDS(system.file("extdata", "example_data.rds", package = "RhyME"))
# A quick smoke run; EB estimation requires at least 1,000 CpGs.
fit <- fitCosinor(example$se[seq_len(20), ], subject_col = "subject",
                  prior = "none")
print(head(SummarizedExperiment::rowData(fit), 5))
print(S4Vectors::metadata(fit)$RhyME$fitSummary)
print(plotCosinorFit(fit))

# To estimate the EB prior on the complete example:
# fit_eb <- fitCosinor(example$se, subject_col = "subject", prior = "auto")
