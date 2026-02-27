# =============================================================================
# Exported functions
# =============================================================================

#' Build Multivariate State-Space Embedding
#'
#' @description
#' Constructs a delayed-coordinate multivariate state-space embedding for each
#' participant subseries.
#'
#' The function is compatible with objects returned by [resample_timeseries()].
#' When multiple participants/subseries are present, each subseries is embedded
#' independently and then combined.
#'
#' @param data A data frame, or an object returned by [resample_timeseries()].
#' @param delay `[integer(1): 1]`\cr
#'   Delay (in rows) between lagged coordinates.
#' @param dimensions `[integer(1): 2]`\cr
#'   Number of delayed coordinates per variable.
#'   Variable columns are detected automatically as all columns except common
#'   metadata identifiers (`ID`, `Series_ID`, `time`, and internal helper
#'   columns). Non-numeric detected variable columns trigger an error.
#'
#' @return An object of class `experdyn_statespace` — a named list with:
#'
#'   | Element | Description |
#'   |---------|-------------|
#'   | `embedded` | Embedded data frame with lagged variables and series identifiers. |
#'   | `delay` | Delay used for embedding. |
#'   | `dimensions` | Number of embedding dimensions used. |
#'   | `variables` | Names of embedded variables. |
#'   | `n_series_input` | Number of participant subseries processed. |
#'   | `n_series_embedded` | Number of subseries with sufficient points to embed. |
#'
#' @export
make_statespace <- function(
    data,
    delay = 1,
    dimensions = 2
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

    if (
        !is.numeric(delay) ||
            length(delay) != 1 ||
            delay < 1 ||
            delay != as.integer(delay)
    ) {
        stop("`delay` must be a single positive integer.", call. = FALSE)
    }
    if (
        !is.numeric(dimensions) ||
            length(dimensions) != 1 ||
            dimensions < 1 ||
            dimensions != as.integer(dimensions)
    ) {
        stop("`dimensions` must be a single positive integer.", call. = FALSE)
    }

    # Ensure identifiers
    if (!"ID" %in% colnames(df)) {
        df <- dplyr::mutate(df, ID = "All")
    }
    if (!"Series_ID" %in% colnames(df)) {
        df <- df |>
            dplyr::group_by(ID) |>
            dplyr::mutate(Series_ID = paste(ID, 1L, sep = "_")) |>
            dplyr::ungroup()
    }

    # Resolve variable columns automatically
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
    variable_cols <- setdiff(colnames(df), metadata_cols)

    if (length(variable_cols) == 0) {
        stop(
            "No variable columns selected for state-space embedding.",
            call. = FALSE
        )
    }

    is_numeric_var <- vapply(df[variable_cols], is.numeric, logical(1))
    if (!all(is_numeric_var)) {
        bad_cols <- variable_cols[!is_numeric_var]
        stop(
            sprintf(
                "Detected non-numeric variable columns: %s",
                paste(bad_cols, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    # Order if time exists
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

    # --- 2. Embed each subseries independently ---
    series_keys <- df |>
        dplyr::distinct(ID, Series_ID)

    embedded_chunks <- vector("list", nrow(series_keys))
    embedded_count <- 0L

    for (series_i in seq_len(nrow(series_keys))) {
        current_id <- series_keys$ID[[series_i]]
        current_series <- series_keys$Series_ID[[series_i]]

        series_df <- df |>
            dplyr::filter(ID == current_id, Series_ID == current_series)

        ts_matrix <- as.matrix(series_df[, variable_cols, drop = FALSE])
        n_obs <- nrow(ts_matrix)
        n_var <- ncol(ts_matrix)

        n_embed <- n_obs - (dimensions - 1) * delay
        if (n_embed <= 0) {
            next
        }

        embed_mat <- matrix(
            NA_real_,
            nrow = n_embed,
            ncol = dimensions * n_var
        )

        col_names <- character(dimensions * n_var)
        col_counter <- 1L

        for (var_i in seq_len(n_var)) {
            for (d in seq_len(dimensions)) {
                lag <- (d - 1L) * delay
                embed_mat[, col_counter] <- ts_matrix[
                    (1 + lag):(n_embed + lag),
                    var_i
                ]
                col_names[col_counter] <- paste0(
                    variable_cols[[var_i]],
                    "_lag",
                    lag
                )
                col_counter <- col_counter + 1L
            }
        }

        embed_df <- as.data.frame(embed_mat, stringsAsFactors = FALSE)
        colnames(embed_df) <- col_names

        # Reference row/time corresponds to latest lag in each embedded vector
        ref_idx <- seq_len(n_embed) + (dimensions - 1L) * delay

        embed_df <- dplyr::mutate(
            embed_df,
            ID = current_id,
            Series_ID = current_series,
            row_index = ref_idx,
            .before = 1
        )

        if ("time" %in% colnames(series_df)) {
            embed_df <- dplyr::mutate(
                embed_df,
                time = series_df$time[ref_idx],
                .after = Series_ID
            )
        }

        embedded_count <- embedded_count + 1L
        embedded_chunks[[embedded_count]] <- embed_df
    }

    if (embedded_count == 0L) {
        embedded <- data.frame()
    } else {
        embedded <- dplyr::bind_rows(embedded_chunks[seq_len(embedded_count)])
    }

    structure(
        list(
            embedded = embedded,
            delay = as.integer(delay),
            dimensions = as.integer(dimensions),
            variables = variable_cols,
            n_series_input = nrow(series_keys),
            n_series_embedded = embedded_count
        ),
        class = "experdyn_statespace"
    )
}

# =============================================================================
# Methods for experdyn_statespace
# =============================================================================

#' Print method for `experdyn_statespace` objects
#'
#' Displays a concise summary of the multivariate embedded state space.
#'
#' @param x An `experdyn_statespace` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.experdyn_statespace <- function(x, ...) {
    cat("========= Multivariate State Space =========\n")
    cat(sprintf("  Delay (rows)         : %d\n", x$delay))
    cat(sprintf("  Dimensions           : %d\n", x$dimensions))
    cat(sprintf("  Variables embedded   : %d\n", length(x$variables)))
    cat(sprintf("  Series processed     : %d\n", x$n_series_input))
    cat(sprintf("  Series embedded      : %d\n", x$n_series_embedded))
    cat(sprintf("  Embedded rows        : %d\n", nrow(x$embedded)))
    cat("  Data frame           : use `$embedded` to extract\n")
    cat("============================================\n")
    invisible(x)
}
