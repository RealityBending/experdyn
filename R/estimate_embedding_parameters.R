# =============================================================================
# Exported functions
# =============================================================================

#' Estimate Optimal State-Space Embedding Parameters
#'
#' @description
#' Estimates delay and embedding dimension (`m`) for each
#' variable time series using functions from `nonlinearTseries`.
#'
#' The function is compatible with outputs from [resample_timeseries()] and
#' automatically accommodates participant-level subseries (`Series_ID`) when
#' available.
#'
#' @param data A data frame, or an object returned by [resample_timeseries()].
#'   All numeric columns other than `ID`, `Series_ID`, and `time` are treated
#'   as variables to analyse.
#' @param max_delay `[integer(1): 10]`\cr
#'   Maximum delay considered for AMI-based delay estimation.
#' @param max_dim `[integer(1): 10]`\cr
#'   Maximum embedding dimension considered for false-nearest-neighbours
#'   estimation.
#' @param ami_selection_method `[character(1): "first.e.decay"]`\cr
#'   Selection rule passed to [nonlinearTseries::timeLag()].
#' @param fnn_threshold `[numeric(1): 0.05]`\cr
#'   False-neighbour threshold passed to
#'   [nonlinearTseries::estimateEmbeddingDim()].
#' @param plot `[logical(1): TRUE]`\cr
#'   If `TRUE`, stores a ggplot summary of per-series delay and `m` estimates
#'   in `$plot`.
#'
#' @return An object of class `experdyn_embedding` — a named list with:
#'
#'   | Element | Description |
#'   |---------|-------------|
#'   | `recommended_delay` | Overall recommended delay (median across estimates). |
#'   | `recommended_dim` | Overall recommended embedding dimension (max across variable medians). |
#'   | `variable_summary` | Per-variable summary table (central tendency + variability). |
#'   | `series_estimates` | Per-series/per-variable raw estimates. |
#'   | `overall_variability` | Overall variability metrics for delay and `m`. |
#'   | `n_variables` | Number of analysed variables. |
#'   | `n_series` | Number of analysed series/subseries. |
#'   | `plot` | A `ggplot` object (or `NULL` if `plot = FALSE`). |
#'
#' @export
estimate_embedding_parameters <- function(
    data,
    max_delay = 10,
    max_dim = 10,
    ami_selection_method = "first.e.decay",
    fnn_threshold = 0.05,
    plot = TRUE
) {
    # --- 1. Resolve input ---
    if (inherits(data, "experdyn_resampled")) {
        df <- data$resampled
    } else {
        df <- data
    }

    if (!is.data.frame(df)) {
        stop(
            "`data` must be a data frame or an `experdyn_resampled` object.",
            call. = FALSE
        )
    }

    # Variables: all numeric columns that are not metadata identifiers.
    metadata_cols <- c(
        "ID",
        "Series_ID",
        "time",
        "time_diff",
        "subseries_num",
        "is_gap",
        ".row_index",
        "row_index"
    )
    variable_cols <- setdiff(
        colnames(df)[vapply(df, is.numeric, logical(1L))],
        metadata_cols
    )
    if (length(variable_cols) == 0) {
        stop("No analysable numeric variable columns found.", call. = FALSE)
    }

    # Ensure participant and subseries identifiers exist
    if (!"ID" %in% colnames(df)) {
        df <- dplyr::mutate(df, ID = "All")
    }
    if (!"Series_ID" %in% colnames(df)) {
        df <- df |>
            dplyr::group_by(ID) |>
            dplyr::mutate(Series_ID = paste(ID, 1L, sep = "_")) |>
            dplyr::ungroup()
    }

    # Optional ordering by time if present
    if ("time" %in% colnames(df)) {
        df <- df |>
            dplyr::mutate(time = as.POSIXct(time)) |>
            dplyr::arrange(ID, Series_ID, time)
    } else {
        df <- df |>
            dplyr::group_by(ID, Series_ID) |>
            dplyr::mutate(.row_index = dplyr::row_number()) |>
            dplyr::ungroup()
    }

    # --- 2. Estimate parameters per variable x subseries ---
    series_keys <- df |>
        dplyr::distinct(ID, Series_ID)

    estimates <- vector(
        "list",
        length = nrow(series_keys) * length(variable_cols)
    )
    idx <- 1L

    for (series_row in seq_len(nrow(series_keys))) {
        current_id <- series_keys$ID[[series_row]]
        current_series <- series_keys$Series_ID[[series_row]]

        series_df <- df |>
            dplyr::filter(ID == current_id, Series_ID == current_series)

        for (var in variable_cols) {
            ts <- series_df[[var]]
            ts <- ts[!is.na(ts)]

            # delay via AMI
            delay_est <- suppressWarnings(tryCatch(
                nonlinearTseries::timeLag(
                    ts,
                    technique = "ami",
                    selection.method = ami_selection_method,
                    lag.max = max_delay,
                    do.plot = FALSE
                ),
                error = function(e) NA_real_
            ))
            if (is.na(delay_est) || delay_est <= 0) {
                delay_est <- 1
            }

            # m via FNN
            dim_est <- suppressWarnings(tryCatch(
                nonlinearTseries::estimateEmbeddingDim(
                    ts,
                    time.lag = as.integer(delay_est),
                    max.embedding.dim = max_dim,
                    threshold = fnn_threshold,
                    do.plot = FALSE
                ),
                error = function(e) NA_real_
            ))
            if (is.na(dim_est) || dim_est <= 0) {
                dim_est <- 2
            }

            estimates[[idx]] <- data.frame(
                ID = current_id,
                Series_ID = current_series,
                Variable = var,
                n_points = length(ts),
                delay = as.numeric(delay_est),
                dim = as.numeric(dim_est),
                stringsAsFactors = FALSE
            )
            idx <- idx + 1L
        }
    }

    series_estimates <- dplyr::bind_rows(estimates)

    # --- 3. Summarise per variable ---
    variable_summary <- series_estimates |>
        dplyr::group_by(Variable) |>
        dplyr::summarise(
            n_series = dplyr::n(),
            delay_median = stats::median(delay, na.rm = TRUE),
            delay_mean = mean(delay, na.rm = TRUE),
            delay_sd = stats::sd(delay, na.rm = TRUE),
            delay_iqr = stats::IQR(delay, na.rm = TRUE),
            dim_median = stats::median(dim, na.rm = TRUE),
            dim_mean = mean(dim, na.rm = TRUE),
            dim_sd = stats::sd(dim, na.rm = TRUE),
            dim_iqr = stats::IQR(dim, na.rm = TRUE),
            .groups = "drop"
        )

    recommended_delay <- as.integer(round(stats::median(
        series_estimates$delay,
        na.rm = TRUE
    )))
    recommended_dim <- as.integer(ceiling(max(
        variable_summary$dim_median,
        na.rm = TRUE
    )))

    overall_variability <- list(
        delay_sd = stats::sd(series_estimates$delay, na.rm = TRUE),
        delay_iqr = stats::IQR(series_estimates$delay, na.rm = TRUE),
        dim_sd = stats::sd(series_estimates$dim, na.rm = TRUE),
        dim_iqr = stats::IQR(series_estimates$dim, na.rm = TRUE)
    )

    # --- 4. Optional plot ---
    p <- NULL
    if (plot) {
        plot_data <- series_estimates |>
            tidyr::pivot_longer(
                cols = c(delay, dim),
                names_to = "Parameter",
                values_to = "Estimate"
            ) |>
            dplyr::mutate(
                Parameter = dplyr::case_match(
                    Parameter,
                    "delay" ~ "Delay",
                    "dim" ~ "Dimension (m)",
                    .default = Parameter
                )
            )

        p <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(x = Variable, y = Estimate)
        ) +
            ggplot2::geom_boxplot(
                fill = "grey90",
                colour = "grey35",
                outlier.shape = NA,
                linewidth = 0.3
            ) +
            ggplot2::geom_jitter(
                ggplot2::aes(colour = as.factor(ID)),
                width = 0.12,
                height = 0.1,
                alpha = 0.7,
                shape = 4,
                stroke = 0.9,
                size = 2.0
            ) +
            ggplot2::facet_wrap(
                ~Parameter,
                scales = "free_y",
                ncol = 1
            ) +
            ggplot2::scale_y_continuous(
                breaks = function(x) {
                    seq(
                        floor(min(x, na.rm = TRUE)),
                        ceiling(max(x, na.rm = TRUE)),
                        by = 1
                    )
                }
            ) +
            ggplot2::labs(
                x = "Variable",
                y = "Estimate",
                colour = "Participant"
            ) +
            ggplot2::theme_minimal()
    }

    structure(
        list(
            recommended_delay = recommended_delay,
            recommended_dim = recommended_dim,
            variable_summary = variable_summary,
            series_estimates = series_estimates,
            overall_variability = overall_variability,
            n_variables = length(variable_cols),
            n_series = nrow(series_keys),
            plot = p
        ),
        class = "experdyn_embedding"
    )
}

# =============================================================================
# Methods for experdyn_embedding
# =============================================================================

#' Print method for `experdyn_embedding` objects
#'
#' Displays recommended embedding parameters and variability diagnostics from
#' [estimate_embedding_parameters()].
#'
#' @param x An `experdyn_embedding` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.experdyn_embedding <- function(x, ...) {
    cat("===== Embedding Parameter Estimation =====\n")
    cat(sprintf("  Recommended delay       : %d\n", x$recommended_delay))
    cat(sprintf("  Recommended dim (m)     : %d\n", x$recommended_dim))
    cat(sprintf("  Variables analysed      : %d\n", x$n_variables))
    cat(sprintf("  Series/subseries used   : %d\n", x$n_series))
    cat("\n")

    cat("  Variability across all variable-series estimates:\n")
    cat(sprintf(
        "    Delay SD / IQR        : %.2f / %.2f\n",
        x$overall_variability$delay_sd,
        x$overall_variability$delay_iqr
    ))
    cat(sprintf(
        "    Dim SD / IQR          : %.2f / %.2f\n",
        x$overall_variability$dim_sd,
        x$overall_variability$dim_iqr
    ))
    cat("\n")

    cat("  Per-variable summary:\n")
    print(
        x$variable_summary |>
            dplyr::select(
                Variable,
                n_series,
                delay_median,
                delay_iqr,
                dim_median,
                dim_iqr
            ),
        row.names = FALSE
    )

    cat("\n")
    cat("  Plot                  : use `$plot` to inspect estimates\n")
    cat("==========================================\n")
    invisible(x)
}
