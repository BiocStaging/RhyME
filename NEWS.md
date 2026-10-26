# RhyME 0.99.2

- Respect the registered BiocParallel backend and prefix result columns with
  RhyME_ to preserve existing row annotations.
- Add unit tests for automatic EB estimation, p-value separation, BH
  adjustment, and phase recovery.
- Update the vignette with rhythmic and non-rhythmic examples, evaluated
  automatic EB fits, MetaCycle comparisons, and methylation array input.
- Declare MetaCycle and minfi in Suggests, and bundle the parameters and
  sampling design needed to regenerate the example data.
- Use Bioconductor installation instructions in the README.

# RhyME 0.99.1

- Replace TMB-based fitting with penalized quasi-likelihood (PQL), using a
  null-model score test for rhythmicity and a separate full-model fit for
  rhythm parameters.

# RhyME 0.99.0

- Initial submission with TMB-based beta cosinor models, optional subject
  random intercepts, and empirical Bayes precision shrinkage.
