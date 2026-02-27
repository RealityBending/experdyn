# =============================================================================
# Internal helpers
# =============================================================================

#' Count run lengths of TRUE values meeting a minimum length
#'
#' @param x Logical vector.
#' @param min_len Minimum run length to retain.
#' @return Integer vector of qualifying run lengths.
#' @keywords internal
.rqa_count_runs <- function(x, min_len) {
    if (!any(x, na.rm = TRUE)) {
        return(integer(0L))
    }
    r <- rle(as.logical(x))
    r$lengths[r$values & !is.na(r$values) & r$lengths >= min_len]
}

#' Core RQA computation for a single embedded matrix
#'
#' @param x Numeric matrix where each row is a state-space point.
#' @param r Euclidean distance threshold (applied to z-scored data).
#' @param min_diag Minimum diagonal line length for DET/L/ENTR.
#' @param min_vert Minimum vertical line length for LAM/TT.
#'
#' @return Named list of RQA measures plus the sparse recurrence data frame.
#' @keywords internal
.rqa_single <- function(x, r, min_diag, min_vert) {
    x <- x[stats::complete.cases(x), , drop = FALSE]
    n <- nrow(x)

    if (n < (min_diag + 2L)) {
        return(list(
            n_obs = n,
            RR = NA_real_,
            DET = NA_real_,
            L = NA_real_,
            Lmax = NA_integer_,
            ENTR = NA_real_,
            LAM = NA_real_,
            TT = NA_real_,
            Vmax = NA_integer_,
            recurrence_df = NULL,
            status = "too_short"
        ))
    }

    x_scaled <- scale(x)
    x_scaled[is.na(x_scaled)] <- 0

    # n x n Euclidean distance matrix → binary recurrence matrix
    d <- as.matrix(stats::dist(x_scaled, method = "euclidean"))
    R <- d <= r
    diag(R) <- FALSE

    n_recur <- sum(R)
    RR <- n_recur / (n * (n - 1L))

    # --- Diagonal lines (DET, L, Lmax, ENTR) ---
    diag_lengths <- integer(0L)
    for (k in seq_len(n - 1L)) {
        idx_a <- seq_len(n - k)
        idx_b <- seq.int(k + 1L, n)
        diag_lengths <- c(
            diag_lengths,
            .rqa_count_runs(R[cbind(idx_a, idx_b)], min_diag), # upper
            .rqa_count_runs(R[cbind(idx_b, idx_a)], min_diag) # lower
        )
    }
    n_det_pts <- sum(diag_lengths)
    DET <- if (n_recur > 0L && n_det_pts > 0L) n_det_pts / n_recur else 0
    L <- if (length(diag_lengths) > 0L) mean(diag_lengths) else NA_real_
    Lmax <- if (length(diag_lengths) > 0L) max(diag_lengths) else NA_integer_
    ENTR <- if (length(diag_lengths) > 1L) {
        p_k <- tabulate(diag_lengths) / length(diag_lengths)
        p_k <- p_k[p_k > 0]
        -sum(p_k * log(p_k))
    } else {
        NA_real_
    }

    # --- Vertical lines (LAM, TT, Vmax) ---
    vert_lengths <- integer(0L)
    for (j in seq_len(n)) {
        vert_lengths <- c(vert_lengths, .rqa_count_runs(R[, j], min_vert))
    }
    n_lam_pts <- sum(vert_lengths)
    LAM <- if (n_recur > 0L && n_lam_pts > 0L) n_lam_pts / n_recur else 0
    TT <- if (length(vert_lengths) > 0L) mean(vert_lengths) else NA_real_
    Vmax <- if (length(vert_lengths) > 0L) max(vert_lengths) else NA_integer_

    # Sparse recurrence representation for plotting (only TRUE positions)
    rec_idx <- which(R, arr.ind = TRUE)
    recurrence_df <- if (nrow(rec_idx) > 0L) {
        data.frame(i = rec_idx[, 1L], j = rec_idx[, 2L])
    } else {
        data.frame(i = integer(0L), j = integer(0L))
    }

    list(
        n_obs = n,
        RR = RR,
        DET = DET,
        L = L,
        Lmax = as.integer(Lmax),
        ENTR = ENTR,
        LAM = LAM,
        TT = TT,
        Vmax = as.integer(Vmax),
        recurrence_df = recurrence_df,
        status = "ok"
    )
}


# =============================================================================
# Exported functions
# =============================================================================

#' Compute Recurrence Quantification Analysis per Participant
#'
#' @description
#' Computes standard Recurrence Quantification Analysis (RQA) measures from a
#' multivariate state-space embedding, separately for each participant.
#'
#' RQA characterises the temporal organisation of a dynamical system by
#' examining how often and in what patterns the system revisits previous states.
#' The function accepts output from [make_statespace()] directly — each row of
#' the embedded matrix is treated as a point in the reconstructed state space.
#'
#' When multiple subseries exist per participant (e.g., from overnight breaks
#' detected by [resample_timeseries()]), RQA is first computed within each
#' subseries and then averaged to yield participant-level estimates.
#'
#' @param data An `experdyn_statespace` object returned by [make_statespace()],
#'   or a data frame containing lagged columns named `*_lag0`, `*_lag1`, etc.
#' @param r `[numeric(1): 1.0]`\cr
#'   Euclidean distance threshold for the recurrence matrix. Data are z-scored
#'   within each subseries before distances are computed, so `r` is on a
#'   standardised scale.
#' @param min_diag `[integer(1): 2]`\cr
#'   Minimum diagonal line length included in DET, L, and ENTR calculations.
#' @param min_vert `[integer(1): 2]`\cr
#'   Minimum vertical line length included in LAM and TT calculations.
#' @param participant_col `[character(1): "ID"]`\cr
#'   Participant identifier column.
#' @param series_col `[character(1): "Series_ID"]`\cr
#'   Subseries identifier column.
#' @param time_col `[character(1): "time"]`\cr
#'   Optional time column used for ordering within subseries.
#' @param verbose `[logical(1): TRUE]`\cr
#'   If `TRUE`, prints a short completion summary.
#'
#' @return An object of class `experdyn_rqa` — a named list with:
#'
#'   | Element | Description |
#'   |---------|-------------|
#'   | `participant_rqa` | One row per participant with cross-series mean RQA measures. |
#'   | `series_rqa` | One row per participant-subseries with per-series RQA measures. |
#'   | `recurrence_data` | Sparse data frame of recurrence point coordinates for plotting. |
#'   | `r` | Matching tolerance used. |
#'   | `min_diag` | Minimum diagonal line length used. |
#'   | `min_vert` | Minimum vertical line length used. |
#'
#'   **RQA measures** (columns in `participant_rqa` and `series_rqa`):
#'
#'   | Measure | Description |
#'   |---------|-------------|
#'   | `RR` | Recurrence Rate — proportion of recurrent points. |
#'   | `DET` | Determinism — proportion of recurrences forming diagonal lines ≥ `min_diag`. |
#'   | `L` | Mean diagonal line length. |
#'   | `Lmax` | Longest diagonal line. |
#'   | `ENTR` | Shannon entropy of diagonal line length distribution. |
#'   | `LAM` | Laminarity — proportion of recurrences forming vertical lines ≥ `min_vert`. |
#'   | `TT` | Trapping Time — mean vertical line length. |
#'   | `Vmax` | Longest vertical line. |
#'
#' @seealso
#' * [make_statespace()] for building the embedding passed to this function.
#' * [plot_rqa()] for visualising the recurrence structure.
#' * [compute_mvsampen()] for sample entropy, a complementary complexity measure.
#'
#' @examples
#' set.seed(1)
#' n <- 30
#' t_obs <- as.POSIXct("2025-01-01 09:00:00") +
#'   cumsum(sample(60 * c(60, 90), n, replace = TRUE))
#' df <- data.frame(
#'   time        = t_obs,
#'   Component_1 = cumsum(rnorm(n)),
#'   Component_2 = cumsum(rnorm(n))
#' )
#' rs  <- resample_timeseries(df, time_col = "time", verbose = FALSE, plot = FALSE)
#' ss  <- make_statespace(rs)
#' rqa <- compute_rqa(ss, r = 1.2, verbose = FALSE)
#' rqa
#'
#' @export
compute_rqa <- function(
    data,
    r = 1.0,
    min_diag = 2L,
    min_vert = 2L,
    participant_col = "ID",
    series_col = "Series_ID",
    time_col = "time",
    verbose = TRUE
) {
    from_statespace <- inherits(data, "experdyn_statespace")
    df <- if (from_statespace) data$embedded else data

    if (!is.data.frame(df)) {
        stop(
            "`data` must be an `experdyn_statespace` object or an embedded data frame.",
            call. = FALSE
        )
    }
    if (!is.numeric(r) || length(r) != 1L || r <= 0) {
        stop("`r` must be a single positive numeric value.", call. = FALSE)
    }

    if (!participant_col %in% colnames(df)) {
        df[[participant_col]] <- "All"
    }
    if (!series_col %in% colnames(df)) {
        df[[series_col]] <- paste(df[[participant_col]], 1L, sep = "_")
    }

    lag_cols <- grep("_lag\\d+$", colnames(df), value = TRUE)
    if (length(lag_cols) == 0L) {
        stop(
            "No lagged embedding columns found. Provide output from `make_statespace()`.",
            call. = FALSE
        )
    }

    # Order by time when present
    if (time_col %in% colnames(df)) {
        df[[time_col]] <- as.POSIXct(df[[time_col]])
        ord <- order(df[[participant_col]], df[[series_col]], df[[time_col]])
        df <- df[ord, , drop = FALSE]
    }

    min_diag <- as.integer(min_diag)
    min_vert <- as.integer(min_vert)

    keys <- unique(df[, c(participant_col, series_col), drop = FALSE])
    series_rows <- vector("list", nrow(keys))
    recurrence_chunks <- vector("list", nrow(keys))

    for (i in seq_len(nrow(keys))) {
        pid <- keys[[participant_col]][i]
        sid <- keys[[series_col]][i]

        idx <- df[[participant_col]] == pid & df[[series_col]] == sid
        x <- as.matrix(df[idx, lag_cols, drop = FALSE])

        res <- .rqa_single(x, r = r, min_diag = min_diag, min_vert = min_vert)

        series_rows[[i]] <- data.frame(
            ID = pid,
            Series_ID = sid,
            n_obs = res$n_obs,
            RR = res$RR,
            DET = res$DET,
            L = res$L,
            Lmax = res$Lmax,
            ENTR = res$ENTR,
            LAM = res$LAM,
            TT = res$TT,
            Vmax = res$Vmax,
            status = res$status,
            stringsAsFactors = FALSE
        )

        if (!is.null(res$recurrence_df) && nrow(res$recurrence_df) > 0L) {
            recurrence_chunks[[i]] <- dplyr::mutate(
                res$recurrence_df,
                ID = pid,
                Series_ID = sid
            )
        }
    }

    series_rqa <- dplyr::bind_rows(series_rows)

    participant_rqa <- series_rqa |>
        dplyr::filter(status == "ok") |>
        dplyr::group_by(ID) |>
        dplyr::summarise(
            n_series = dplyr::n(),
            n_series_valid = sum(!is.na(RR)),
            # Weight averaged measures by subseries length so that longer
            # (more data-rich) subseries contribute proportionally more.
            RR = {
                ok <- !is.na(RR)
                if (!any(ok)) {
                    NA_real_
                } else {
                    stats::weighted.mean(RR[ok], w = n_obs[ok])
                }
            },
            DET = {
                ok <- !is.na(DET)
                if (!any(ok)) {
                    NA_real_
                } else {
                    stats::weighted.mean(DET[ok], w = n_obs[ok])
                }
            },
            L = {
                ok <- !is.na(L)
                if (!any(ok)) {
                    NA_real_
                } else {
                    stats::weighted.mean(L[ok], w = n_obs[ok])
                }
            },
            Lmax = {
                x <- Lmax[!is.na(Lmax)]
                if (length(x) == 0L) NA_integer_ else as.integer(max(x))
            },
            ENTR = {
                ok <- !is.na(ENTR)
                if (!any(ok)) {
                    NA_real_
                } else {
                    stats::weighted.mean(ENTR[ok], w = n_obs[ok])
                }
            },
            LAM = {
                ok <- !is.na(LAM)
                if (!any(ok)) {
                    NA_real_
                } else {
                    stats::weighted.mean(LAM[ok], w = n_obs[ok])
                }
            },
            TT = {
                ok <- !is.na(TT)
                if (!any(ok)) {
                    NA_real_
                } else {
                    stats::weighted.mean(TT[ok], w = n_obs[ok])
                }
            },
            Vmax = {
                x <- Vmax[!is.na(Vmax)]
                if (length(x) == 0L) NA_integer_ else as.integer(max(x))
            },
            .groups = "drop"
        )

    non_null <- !vapply(recurrence_chunks, is.null, logical(1L))
    recurrence_data <- dplyr::bind_rows(recurrence_chunks[non_null])

    if (verbose) {
        cat("--- RQA Complete ---\n")
        cat(sprintf("  Participants  : %d\n", nrow(participant_rqa)))
        cat(sprintf("  Series        : %d\n", nrow(series_rqa)))
        cat(sprintf("  Tolerance (r) : %.3f\n", r))
        cat(sprintf("  min_diag / min_vert : %d / %d\n", min_diag, min_vert))
        cat("--------------------\n")
    }

    structure(
        list(
            participant_rqa = participant_rqa,
            series_rqa = series_rqa,
            recurrence_data = recurrence_data,
            r = r,
            min_diag = min_diag,
            min_vert = min_vert
        ),
        class = "experdyn_rqa"
    )
}


#' Plot Recurrence Plots Faceted by Participant or Series
#'
#' @description
#' Renders a classical recurrence plot — a square binary image where a filled
#' point at position *(i, j)* means the system revisited state *j* at time *i*.
#' Diagonal lines indicate deterministic structure; vertical lines indicate
#' laminar (intermittent) behaviour.
#'
#' The function accepts the `experdyn_rqa` object produced by [compute_rqa()],
#' which stores the sparse recurrence coordinates computed during analysis.
#'
#' @param data An `experdyn_rqa` object returned by [compute_rqa()].
#' @param facet_by `[character(1): "participant"]`\cr
#'   Whether to create one panel per `"participant"` (default) or per
#'   `"series"` (subseries-level, useful for exploratory inspection).
#' @param max_n `[integer(1): 300]`\cr
#'   Maximum time index retained per panel for plotting. Points with *i* or *j*
#'   greater than `max_n` are trimmed. Set to `Inf` (or `0`) to plot all points.
#' @param point_colour `[character(1): "black"]`\cr
#'   Fill colour of recurrent points.
#' @param background_colour `[character(1): "white"]`\cr
#'   Background (non-recurrent) colour.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @seealso [compute_rqa()] for computing the RQA object passed here.
#'
#' @export
plot_rqa <- function(
    data,
    facet_by = c("participant", "series"),
    max_n = 300L,
    point_colour = "black",
    background_colour = "white"
) {
    if (!inherits(data, "experdyn_rqa")) {
        stop(
            "`data` must be an `experdyn_rqa` object returned by `compute_rqa()`.",
            call. = FALSE
        )
    }

    facet_by <- match.arg(facet_by)

    rec_df <- data$recurrence_data
    if (is.null(rec_df) || nrow(rec_df) == 0L) {
        stop(
            "No recurrence data available. Re-run `compute_rqa()` to regenerate.",
            call. = FALSE
        )
    }

    # Trim to max_n per panel
    if (!is.null(max_n) && is.finite(max_n) && max_n > 0) {
        rec_df <- rec_df[rec_df$i <= max_n & rec_df$j <= max_n, , drop = FALSE]
    }

    if (nrow(rec_df) == 0L) {
        stop(
            sprintf(
                "No recurrence points remain after trimming to max_n = %d.",
                max_n
            ),
            call. = FALSE
        )
    }

    facet_col <- if (facet_by == "participant") "ID" else "Series_ID"

    # Ensure every participant/series gets a panel, even those with zero
    # recurrences (otherwise facet_wrap silently drops them).
    all_levels <- if (facet_by == "participant") {
        as.character(data$participant_rqa$ID)
    } else {
        as.character(data$series_rqa$Series_ID)
    }
    if (
        length(all_levels) > 0L &&
            facet_col %in% colnames(rec_df) &&
            nrow(rec_df) > 0L
    ) {
        rec_df[[facet_col]] <- factor(rec_df[[facet_col]], levels = all_levels)
    }

    p <- ggplot2::ggplot(
        rec_df,
        ggplot2::aes(x = j, y = i)
    ) +
        ggplot2::geom_tile(fill = point_colour, width = 1, height = 1) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
            expand = c(0, 0)
        ) +
        ggplot2::labs(
            x = "Time index j",
            y = "Time index i"
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(
                fill = background_colour,
                colour = NA
            ),
            panel.grid = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_line(linewidth = 0.3),
            strip.text = ggplot2::element_text(face = "bold"),
            aspect.ratio = 1
        )

    if (facet_col %in% colnames(rec_df)) {
        p <- p +
            ggplot2::facet_wrap(
                stats::as.formula(paste0("~", facet_col)),
                scales = "free",
                drop = FALSE # show panels even for participants with no recurrences
            )
    }

    p
}


# =============================================================================
# Methods for experdyn_rqa
# =============================================================================

#' Print method for `experdyn_rqa` objects
#'
#' Displays a concise summary of RQA results produced by [compute_rqa()].
#'
#' @param x An `experdyn_rqa` object.
#' @param digits `[integer(1): 3]`\cr Number of decimal places shown.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.experdyn_rqa <- function(x, digits = 3L, ...) {
    cat("===== Recurrence Quantification Analysis =====\n")
    cat(sprintf("  Tolerance (r)          : %.3f\n", x$r))
    cat(sprintf("  Min diagonal (min_diag): %d\n", x$min_diag))
    cat(sprintf("  Min vertical (min_vert): %d\n", x$min_vert))
    cat(sprintf("  Participants           : %d\n", nrow(x$participant_rqa)))
    cat(sprintf("  Series                 : %d\n", nrow(x$series_rqa)))
    cat("\n")

    tbl <- x$participant_rqa
    fmt_cols <- c("RR", "DET", "L", "Lmax", "ENTR", "LAM", "TT", "Vmax")
    fmt_cols <- fmt_cols[fmt_cols %in% colnames(tbl)]
    tbl[fmt_cols] <- lapply(tbl[fmt_cols], round, digits = digits)
    print(as.data.frame(tbl), row.names = FALSE)

    cat("\n  Per-series results : use `$series_rqa`\n")
    cat("  Recurrence data    : use `$recurrence_data` (for custom plots)\n")
    cat("===============================================\n")
    invisible(x)
}
