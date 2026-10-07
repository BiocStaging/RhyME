.fitPQLBatch <- function(beta, design, subject = NULL, prior_phi = NULL) {
    n <- nrow(beta)
    first <- .fitOneCpG(beta[1L, ], design, subject, prior_phi)
    result <- lapply(first, rep, times = n)
    for (i in seq_len(n)[-1L]) {
        fit <- .fitOneCpG(beta[i, ], design, subject, prior_phi)
        for (j in seq_along(result)) result[[j]][i] <- fit[[j]]
    }
    as.data.frame(result, row.names = rownames(beta), stringsAsFactors = FALSE)
}

.fitPQLPilotBatch <- function(beta, design, subject = NULL) {
    vapply(seq_len(nrow(beta)), function(i) {
        tryCatch({
            input <- .preparePQLData(beta[i, ], design, subject)
            fit <- .fitPQLModel(input, null_model = TRUE)
            if (isTRUE(fit$converged)) fit$log_phi else NA_real_
        }, error = function(e) NA_real_)
    }, numeric(1))
}

.fitOneCpG <- function(y, X, subject = NULL, prior_phi = NULL) {
    tryCatch({
        input <- .preparePQLData(y, X, subject)
        null <- .fitPQLModel(input, TRUE, prior_phi)
        full <- .fitPQLModel(input, FALSE, prior_phi)
        number <- function(fit, name) {
            if (is.null(fit[[name]])) NA_real_ else as.numeric(fit[[name]])
        }
        parameters <- vapply(c("b0", "b_sin", "b_cos", "log_phi", "tau"),
                             function(name) number(full, name), numeric(1))
        p <- number(null, "p_value")
        usable_p <- isTRUE(null$converged) && is.finite(p) && p >= 0 && p <= 1
        usable_parameters <- isTRUE(full$converged) && all(is.finite(parameters))
        if (!usable_parameters) parameters[] <- NA_real_
        c(as.list(parameters), list(
            pvalue = if (usable_p) p else NA_real_,
            score_chisq = number(null, "score_chisq"),
            converged = isTRUE(null$converged) && isTRUE(full$converged)
        ))
    }, error = function(e) list(
        b0 = NA_real_, b_sin = NA_real_, b_cos = NA_real_,
        log_phi = NA_real_, tau = NA_real_, pvalue = NA_real_,
        score_chisq = NA_real_, converged = FALSE
    ))
}

.preparePQLData <- function(y, X, subject = NULL) {
    keep <- is.finite(y)
    if (!all(keep)) {
        y <- y[keep]
        X <- X[keep, , drop = FALSE]
        if (!is.null(subject)) subject <- as.integer(factor(subject[keep])) - 1L
        .checkPQLDesign(X, subject)
    }
    list(y = as.numeric(y), design = X,
         subject = if (is.null(subject)) integer(length(y)) else subject,
         fixed_only = is.null(subject))
}

# Catch failures separately for the null and full fits.
.fitPQLModel <- function(input, null_model, prior = NULL) {
    tryCatch({
        pql_working_lmm_cpp(
            input$y, input$design, input$subject,
            null_model = null_model, use_prior = !is.null(prior),
            prior_mean = if (is.null(prior)) 0 else prior$mean,
            prior_sd = if (is.null(prior)) 1 else prior$sd,
            fixed_only = input$fixed_only
        )
    }, error = function(e) list(
        converged = FALSE, status = "error", error = conditionMessage(e)
    ))
}
