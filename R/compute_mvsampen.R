# =============================================================================
# Internal helpers
# =============================================================================

.mv_embed_matrix <- function(x, m) {
    n <- nrow(x)
    p <- ncol(x)
    n_embed <- n - m + 1L
    if (n_embed <= 1L) {
        return(NULL)
    }

    out <- matrix(NA_real_, nrow = n_embed, ncol = p * m)
    col_idx <- 1L

    for (lag in 0:(m - 1L)) {
        block <- x[(1L + lag):(n_embed + lag), , drop = FALSE]
        out[, col_idx:(col_idx + p - 1L)] <- block
        col_idx <- col_idx + p
    }

    out
}

.mv_match_count <- function(embed, tol) {
    if (is.null(embed)) {
        return(0L)
    }
    n <- nrow(embed)
    if (n <= 1L) {
        return(0L)
    }

    count <- 0L
    for (i in 1:(n - 1L)) {
        diffs <- abs(sweep(
            embed[(i + 1L):n, , drop = FALSE],
            2,
            embed[i, ],
            "-"
        ))
        dmax <- apply(diffs, 1, max)
        count <- count + sum(dmax <= tol)
    }
    count
}

.mv_sampen_single <- function(x, m, r) {
    x <- x[stats::complete.cases(x), , drop = FALSE]
    n <- nrow(x)

    if (n <= (m + 1L)) {
        return(list(
            entropy = NA_real_,
            pairs_m = 0L,
            pairs_m1 = 0L,
            matches_m = 0L,
            matches_m1 = 0L,
            n_used = n,
            status = "too_short"
        ))
    }

    x_scaled <- scale(x)
    x_scaled[is.na(x_scaled)] <- 0

    emb_m <- .mv_embed_matrix(x_scaled, m)
    emb_m1 <- .mv_embed_matrix(x_scaled, m + 1L)

    n_m <- nrow(emb_m)
    n_m1 <- nrow(emb_m1)
    pairs_m <- choose(n_m, 2)
    pairs_m1 <- choose(n_m1, 2)

    if (pairs_m <= 0L || pairs_m1 <= 0L) {
        return(list(
            entropy = NA_real_,
            pairs_m = pairs_m,
            pairs_m1 = pairs_m1,
            matches_m = 0L,
            matches_m1 = 0L,
            n_used = n,
            status = "too_short"
        ))
    }

    matches_m <- .mv_match_count(emb_m, tol = r)
    matches_m1 <- .mv_match_count(emb_m1, tol = r)

    b <- matches_m / pairs_m
    a <- matches_m1 / pairs_m1

    entropy <- if (b <= 0) {
        NA_real_
    } else if (a <= 0) {
        Inf
    } else {
        -log(a / b)
    }

    list(
        entropy = entropy,
        pairs_m = pairs_m,
        pairs_m1 = pairs_m1,
        matches_m = matches_m,
        matches_m1 = matches_m1,
        n_used = n,
        status = "ok"
    )
}


# =============================================================================
# Exported functions
# =============================================================================

#' Compute Multivariate Sample Entropy per Participant
#'
#' @description
#' Computes multivariate sample entropy (mvSampEn) from a pre-built
#' state-space embedding, separately for each participant.
#'
#' **Input must be an `experdyn_statespace` object** from [make_statespace()].
#' The function does not perform any additional embedding — each row of the
#' state-space matrix is already a trajectory point that encodes the optimal
#' delay `τ` and dimension `d` chosen by [estimate_embedding_parameters()]
#' and applied by [make_statespace()]. SampEn is then computed directly on
#' this sequence of trajectory points:
#'
#' \itemize{
#'   \item **B** = proportion of all point-pairs \eqn{(i \ne j)} whose
#'         Chebyshev distance is \eqn{\le r} (same-neighbourhood).
#'   \item **A** = proportion of consecutive point-pairs where *both* the
#'         current and next state-space point are simultaneously within
#'         \eqn{r} (trajectory-neighbourhood).
#'   \item \eqn{\text{mvSampEn} = -\log(A / B)}.
#' }
#'
#' For each participant, entropy is first computed within each subseries and
#' then aggregated (weighted by number of template pairs) to yield a single
#' participant-level estimate.
#'
#' @param data An `experdyn_statespace` object returned by [make_statespace()].
#' @param r `[numeric(1): 0.2]`\cr
#'   Matching tolerance. Data are z-scored within each subseries before
#'   matching, so `r` is on a standardised scale.
#' @param participant_col `[character(1): "ID"]`\cr
#'   Participant identifier column.
#' @param series_col `[character(1): "Series_ID"]`\cr
#'   Subseries identifier column.
#' @param time_col `[character(1): "time"]`\cr
#'   Optional time column used for ordering within subseries when present.
#'
#' @return An object of class `experdyn_mvsampen` with:
#'
#'   | Element | Description |
#'   |---------|-------------|
#'   | `participant_entropy` | One row per participant with aggregated entropy. |
#'   | `series_entropy` | One row per participant-subseries entropy estimate. |
#'   | `dimensions` | Statespace embedding dimensions used as input. |
#'   | `r` | Matching tolerance used. |
#'   | `variables` | Variables included in entropy estimation. |
#'
#' @seealso
#' * [make_statespace()] for building the state-space passed to this function.
#' * [estimate_embedding_parameters()] to choose `delay` and `dimensions`.
#' * [compute_rqa()] for recurrence-based complexity measures.
#'
#' @export
compute_mvsampen <- function(
    data,
    r = 0.2,
    participant_col = "ID",
    series_col = "Series_ID",
    time_col = "time"
) {
    if (!inherits(data, "experdyn_statespace")) {
        stop(
            paste0(
                "`data` must be an `experdyn_statespace` object returned by ",
                "`make_statespace()`. Raw data frames are not accepted: build ",
                "the state-space first so that the optimal delay \u03c4 and ",
                "embedding dimension are encoded in the input."
            ),
            call. = FALSE
        )
    }
    df <- data$embedded
    if (!participant_col %in% colnames(df)) {
        df[[participant_col]] <- "All"
    }
    if (!series_col %in% colnames(df)) {
        df[[series_col]] <- paste(df[[participant_col]], 1L, sep = "_")
    }

    lag_cols <- grep("_lag\\d+$", colnames(df), value = TRUE)
    if (length(lag_cols) == 0L) {
        stop(
            "No lagged embedding columns found in the statespace object.",
            call. = FALSE
        )
    }

    lag_values <- as.integer(sub(".*_lag", "", lag_cols))
    dimensions <- length(unique(lag_values))

    lag_info <- data.frame(
        col = lag_cols,
        variable = sub("_lag\\d+$", "", lag_cols),
        lag = lag_values,
        stringsAsFactors = FALSE
    )
    variables <- unique(lag_info$variable)

    # Use all lag columns: each row is a full state-space point encoding delay τ.
    # m_internal = 1 means each state-space point is its own template;
    # extending to m+1 = 2 checks whether consecutive points are also nearby,
    # capturing trajectory-level (rather than point-level) predictability.
    x_cols <- lag_cols
    m_internal <- 1L

    if (!is.numeric(r) || length(r) != 1L || r <= 0) {
        stop("`r` must be a single positive numeric value.", call. = FALSE)
    }

    if (time_col %in% colnames(df)) {
        df[[time_col]] <- as.POSIXct(df[[time_col]])
        ord <- order(df[[participant_col]], df[[series_col]], df[[time_col]])
        df <- df[ord, , drop = FALSE]
    }

    keys <- unique(df[, c(participant_col, series_col), drop = FALSE])
    series_rows <- vector("list", nrow(keys))

    for (i in seq_len(nrow(keys))) {
        pid <- keys[[participant_col]][i]
        sid <- keys[[series_col]][i]

        idx <- df[[participant_col]] == pid & df[[series_col]] == sid
        x <- as.matrix(df[idx, x_cols, drop = FALSE])

        stats <- .mv_sampen_single(x, m = m_internal, r = r)

        series_rows[[i]] <- data.frame(
            ID = pid,
            Series_ID = sid,
            n_used = stats$n_used,
            entropy = stats$entropy,
            pairs_m = stats$pairs_m,
            pairs_m1 = stats$pairs_m1,
            matches_m = stats$matches_m,
            matches_m1 = stats$matches_m1,
            status = stats$status,
            stringsAsFactors = FALSE
        )
    }

    series_entropy <- dplyr::bind_rows(series_rows)

    # Compute SD from the per-series entropy values first, independently of the
    # weighted-mean aggregate, to avoid any name-masking within summarise().
    entropy_sd <- series_entropy |>
        dplyr::group_by(ID) |>
        dplyr::summarise(
            entropy_sd_series = stats::sd(
                entropy[is.finite(entropy)],
                na.rm = TRUE
            ),
            .groups = "drop"
        )

    participant_entropy <- series_entropy |>
        dplyr::group_by(ID) |>
        dplyr::summarise(
            n_series = dplyr::n(),
            n_series_valid = sum(is.finite(entropy), na.rm = TRUE),
            entropy = ifelse(
                sum(is.finite(entropy) & pairs_m > 0) > 0,
                stats::weighted.mean(
                    entropy[is.finite(entropy) & pairs_m > 0],
                    w = pairs_m[is.finite(entropy) & pairs_m > 0],
                    na.rm = TRUE
                ),
                NA_real_
            ),
            .groups = "drop"
        ) |>
        dplyr::left_join(entropy_sd, by = "ID")

    structure(
        list(
            participant_entropy = participant_entropy,
            series_entropy = series_entropy,
            dimensions = as.integer(dimensions),
            r = r,
            variables = variables
        ),
        class = "experdyn_mvsampen"
    )
}

# =============================================================================
# Methods for experdyn_mvsampen
# =============================================================================

#' Print method for `experdyn_mvsampen` objects
#'
#' @param x An `experdyn_mvsampen` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.experdyn_mvsampen <- function(x, ...) {
    cat("===== Multivariate Sample Entropy =====\n")
    cat(sprintf("  Statespace dimensions  : %d\n", x$dimensions))
    cat(sprintf("  Tolerance (r)         : %.3f\n", x$r))
    cat(sprintf("  Variables used        : %d\n", length(x$variables)))
    cat(sprintf("  Participants          : %d\n", nrow(x$participant_entropy)))
    cat(sprintf("  Series                : %d\n", nrow(x$series_entropy)))
    cat("\n")
    print(x$participant_entropy, row.names = FALSE)
    cat("\n")
    cat("  Detailed per-series results: use `$series_entropy`\n")
    cat("========================================\n")
    invisible(x)
}
