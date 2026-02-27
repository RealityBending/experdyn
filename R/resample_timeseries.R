# =============================================================================
# Internal helpers
# =============================================================================

#' Safe Natural Cubic Spline Interpolation
#'
#' A thin wrapper around [stats::spline()] that gracefully handles `NA` values
#' and series that are too short to be interpolated. Output is clamped to the
#' range of valid observations — no extrapolation beyond the data.
#'
#' @param x Numeric vector of known *x* positions (e.g., time in seconds).
#' @param y Numeric vector of known *y* values (possibly containing `NA`s).
#' @param xout Numeric vector of output positions at which to evaluate the spline.
#' @param method Character string passed to [stats::spline()]. Default `"natural"`.
#' @param min_pts Minimum number of non-`NA` observations required; if fewer
#'   are available, a vector of `NA_real_` of length `length(xout)` is returned.
#'
#' @return Numeric vector of interpolated values, the same length as `xout`.
#'   Positions outside `[min(x_valid), max(x_valid)]` are returned as `NA`.
#'
#' @keywords internal
.safe_spline <- function(x, y, xout, method = "natural", min_pts = 4) {
    valid_idx <- !is.na(y)
    x_valid <- x[valid_idx]
    y_valid <- y[valid_idx]

    if (length(y_valid) < min_pts) {
        return(rep(NA_real_, length(xout)))
    }

    out <- rep(NA_real_, length(xout))
    in_range <- xout >= min(x_valid) & xout <= max(x_valid)
    if (any(in_range)) {
        out[in_range] <- stats::spline(
            x = x_valid,
            y = y_valid,
            xout = xout[in_range],
            method = method
        )$y
    }
    out
}

#' Safe Linear Interpolation
#'
#' A thin wrapper around [stats::approx()] that gracefully handles `NA` values
#' and series that are too short to be interpolated. Output is clamped to the
#' range of valid observations — no extrapolation beyond the data.
#'
#' @param x Numeric vector of known *x* positions.
#' @param y Numeric vector of known *y* values (possibly containing `NA`s).
#' @param xout Numeric vector of output positions.
#' @param min_pts Minimum number of non-`NA` observations required.
#'
#' @return Numeric vector of interpolated values, the same length as `xout`.
#'   Positions outside `[min(x_valid), max(x_valid)]` are returned as `NA`.
#'
#' @keywords internal
.safe_approx <- function(x, y, xout, min_pts = 2) {
    valid_idx <- !is.na(y)
    x_valid <- x[valid_idx]
    y_valid <- y[valid_idx]

    if (length(y_valid) < min_pts) {
        return(rep(NA_real_, length(xout)))
    }

    # approx() with rule=1 returns NA outside the range by default
    stats::approx(x = x_valid, y = y_valid, xout = xout, rule = 1)$y
}


# =============================================================================
# Exported functions
# =============================================================================

#' Resample and Interpolate Irregularly Sampled ESM Time Series
#'
#' @description
#' Converts irregularly spaced time series data — typical of Experience Sampling
#' Method (ESM) / Ecological Momentary Assessment (EMA) studies — to a **regular
#' time grid** using natural cubic spline interpolation.
#'
#' The function automatically detects large temporal gaps in the signal and
#' splits the data into contiguous **subseries** before interpolating, so that
#' artefacts are not introduced across overnight breaks or missing assessment
#' periods.
#'
#' @param data A data frame containing at minimum:
#'   * a timestamp column (identified by `time_col`);
#'   * one or more **numeric variable columns** — all numeric columns other
#'     than `time_col` and `participant_col` are treated as variables to
#'     interpolate (no specific naming convention is required);
#'   * optionally, a participant identifier column (identified by `participant_col`).
#' @param time_col `[character(1): "time"]`\cr
#'   Name of the column in `data` that contains the timestamps. Values must be
#'   coercible to `POSIXct`.
#' @param participant_col `[character(1) or NULL: NULL]`\cr
#'   Name of the column in `data` that identifies participants. If `NULL`
#'   (default), all observations are treated as coming from a single participant
#'   labelled `"All"`.
#' @param max_gap `[numeric(1): 6]`\cr
#'   Maximum allowable gap **in hours** between two consecutive observations
#'   before a new subseries is started. Gaps larger than this value (e.g.,
#'   overnight breaks between daily bursts of notifications) mark a boundary
#'   between subseries.
#' @param target_frequency `[numeric(1): 60]`\cr
#'   Desired spacing **in minutes** between consecutive observations on the
#'   regular output grid. For example, `target_frequency = 60` produces one
#'   observation per hour; `target_frequency = 30` produces one every 30 minutes.
#' @param interpolation_method `[character(1): "spline"]`\cr
#'   Interpolation method to use. Either `"spline"` (natural cubic spline via
#'   [stats::spline()]) or `"linear"` (piecewise linear via [stats::approx()]).
#'   Both methods operate strictly within the range of valid (non-`NA`)
#'   observations — no extrapolation is performed beyond the data.
#' @param min_points `[integer(1): 4]`\cr
#'   Minimum number of valid (non-`NA`) observations required within a
#'   subseries to attempt interpolation. Subseries with fewer points are
#'   silently dropped from the output.
#' @param verbose `[logical(1): TRUE]`\cr
#'   If `TRUE`, prints a short processing summary to the console.
#' @param plot `[logical(1): TRUE]`\cr
#'   If `TRUE`, a [ggplot2::ggplot()] object is built showing the interpolated
#'   time series (lines) overlaid with the original raw observations (crosses).
#'   Participants are distinguished by colour; multiple variables
#'   are shown in separate facets. The plot is stored as `$plot` in the returned
#'   object.
#' @param plot_offset `[logical(1): FALSE]`\cr
#'   If `TRUE`, each subseries is shifted so that its x-axis starts at 0 (hours
#'   elapsed since the beginning of that subseries). This makes subseries and
#'   participants directly comparable regardless of their absolute timestamps.
#'
#' @return An object of class `experdyn_resampled` — a named list with:
#'
#'   | Element | Description |
#'   |---------|-------------|
#'   | `resampled` | Data frame of interpolated observations (see below). |
#'   | `plot` | A `ggplot` object (or `NULL` if `plot = FALSE`). |
#'   | `interpolation_method` | Interpolation method used (`"spline"` or `"linear"`). |
#'   | `target_frequency` | Requested spacing in minutes between observations. |
#'   | `sampling_rate` | Achieved rate in observations per hour (`60 / target_frequency`). |
#'   | `max_gap` | Gap threshold in hours used for subseries detection. |
#'   | `n_participants` | Number of unique participants in the output. |
#'   | `n_subseries` | Total number of usable subseries. |
#'   | `n_points_raw` | Number of rows in the input `data`. |
#'   | `n_points_resampled` | Number of rows in `resampled`. |
#'   | `avg_points_per_participant` | Mean resampled observations per participant. |
#'
#'   The `resampled` data frame contains:
#'
#'   | Column | Description |
#'   |--------|-------------|
#'   | `ID` | Participant identifier. |
#'   | `Series_ID` | Subseries label `"<ID>_<n>"`. |
#'   | `time` | Regular grid timestamps (`POSIXct`). |
#'   | `<variable columns>` | Interpolated values, one column per numeric input variable. |
#'
#'   Use `print()` for a human-readable summary, or access `$resampled` for
#'   the underlying data frame.
#'
#' @details
#' ### Processing pipeline
#'
#' **Step 1 — Gap detection.**
#' For each participant, consecutive observations are compared. A new subseries
#' is started whenever the gap exceeds `max_gap` (in hours). The subseries index is
#' stored in `Series_ID`.
#'
#' **Step 2 — Spline interpolation.**
#' Within each qualifying subseries (at least `min_points` observations), each
#' variable column is interpolated on a regular grid spaced at
#' `target_frequency` minutes using the chosen `interpolation_method`
#' (natural cubic spline or piecewise linear). Interpolation is strictly
#' bounded to the range of valid observations — no extrapolation is performed.
#' Grid points that fall outside this range, or where too many raw values are
#' `NA`, are silently dropped.
#'
#' **Step 3 — Summary report.**
#' When `verbose = TRUE`, a brief table is printed showing the number of raw
#' observations, interpolated observations, usable subseries, and the output
#' sampling rate.
#'
#' ### Why natural cubic splines?
#' Natural splines minimise the integrated squared second derivative subject to
#' the constraint of passing through the data, which produces smooth
#' interpolants without oscillatory artefacts at the boundaries. They are a
#' standard choice for resampling physiological and psychological signals.
#'
#' @seealso
#' * [stats::spline()] for the underlying interpolation function.
#' * The `nonlinearTseries` package for downstream complexity analyses.
#'
#' @examples
#' # -----------------------------------------------------------------------
#' # Minimal synthetic ESM dataset
#' # Two participants, 15 irregular observations each, two variables
#' # -----------------------------------------------------------------------
#' set.seed(42)
#' n <- 15
#'
#' # Irregular timestamps: gaps of 1, 1.5, 2, or 3 hours chosen at random
#' t_p1 <- as.POSIXct("2025-01-01 09:00:00") +
#'   cumsum(sample(60 * c(60, 90, 120, 180), n, replace = TRUE))
#' t_p2 <- as.POSIXct("2025-01-05 09:00:00") +
#'   cumsum(sample(60 * c(60, 90, 120, 180), n, replace = TRUE))
#'
#' esm_data <- data.frame(
#'   Participant = c(rep("P01", n), rep("P02", n)),
#'   timestamp   = c(t_p1, t_p2),
#'   valence    = rnorm(2 * n),
#'   arousal    = rnorm(2 * n)
#' )
#'
#' # Resample to a regular 60-minute grid — all numeric columns are used
#' result <- resample_timeseries(
#'   esm_data,
#'   time_col             = "timestamp",
#'   participant_col      = "Participant",
#'   target_frequency     = 60,
#'   verbose              = TRUE
#' )
#' result
#' head(result$resampled)
#'
#' # -----------------------------------------------------------------------
#' # Single-participant usage (participant_col = NULL)
#' # -----------------------------------------------------------------------
#' single <- data.frame(
#'   timestamp   = t_p1,
#'   valence     = rnorm(n)
#' )
#' result_single <- resample_timeseries(
#'   single,
#'   time_col            = "timestamp",
#'   target_frequency    = 30,
#'   verbose             = FALSE
#' )
#'
#' # -----------------------------------------------------------------------
#' # Stricter gap threshold: split on gaps > 4 hours
#' # -----------------------------------------------------------------------
#' resampled_strict <- resample_timeseries(
#'   esm_data,
#'   time_col            = "timestamp",
#'   participant_col     = "Participant",
#'   max_gap             = 4,
#'   target_frequency    = 30,
#'   verbose             = FALSE
#' )
#' # Note: more (shorter) subseries are detected
#' length(unique(resampled_strict$resampled$Series_ID))
#'
#' # -----------------------------------------------------------------------
#' # Arbitrary column names — no "Component_" prefix required
#' # -----------------------------------------------------------------------
#' esm_custom <- data.frame(
#'   id        = c(rep("P01", n), rep("P02", n)),
#'   ts        = c(t_p1, t_p2),
#'   stress    = rnorm(2 * n),
#'   energy    = rnorm(2 * n)
#' )
#' result_custom <- resample_timeseries(
#'   esm_custom,
#'   time_col        = "ts",
#'   participant_col = "id",
#'   verbose         = FALSE
#' )
#'
#' @export
resample_timeseries <- function(
    data,
    time_col = "time",
    participant_col = NULL,
    max_gap = 6,
    target_frequency = 60,
    interpolation_method = "spline",
    min_points = 4,
    verbose = TRUE,
    plot = TRUE,
    plot_offset = FALSE
) {
    # --- Input validation ---
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.", call. = FALSE)
    }
    if (!time_col %in% colnames(data)) {
        stop(
            sprintf("`time_col` ('%s') not found in `data`.", time_col),
            call. = FALSE
        )
    }
    if (!is.null(participant_col) && !participant_col %in% colnames(data)) {
        stop(
            sprintf(
                "`participant_col` ('%s') not found in `data`.",
                participant_col
            ),
            call. = FALSE
        )
    }
    interpolation_method <- match.arg(
        interpolation_method,
        choices = c("spline", "linear")
    )

    # Detect variable columns from the raw input: everything numeric except
    # the time column and the participant identifier.
    excluded_raw <- c(
        time_col,
        if (!is.null(participant_col)) participant_col else character(0L)
    )
    variable_cols <- setdiff(
        colnames(data)[vapply(data, is.numeric, logical(1L))],
        excluded_raw
    )
    if (length(variable_cols) == 0L) {
        stop(
            paste0(
                "No numeric variable columns found in `data` after excluding ",
                "the time column ('",
                time_col,
                "')",
                if (!is.null(participant_col)) {
                    paste0(" and participant column ('", participant_col, "')")
                } else {
                    ""
                },
                "."
            ),
            call. = FALSE
        )
    }

    # --- 1. Normalise column names for internal processing ---
    df <- data
    # Rename the time column to the internal name "time"
    if (time_col != "time") {
        df <- dplyr::rename(df, time = dplyr::all_of(time_col))
    }
    # Resolve participant column → always called "ID" internally
    if (is.null(participant_col)) {
        df <- dplyr::mutate(df, ID = "All")
    } else if (participant_col != "ID") {
        df <- dplyr::rename(df, ID = dplyr::all_of(participant_col))
    }

    df <- df |>
        dplyr::mutate(time = as.POSIXct(time)) |>
        dplyr::arrange(ID, time)

    # Convert target frequency (minutes) to step size for seq()
    step_mins <- target_frequency

    # --- 2. Identify subseries (split on gaps > max_gap hours) ---
    df_series <- df |>
        dplyr::group_by(ID) |>
        dplyr::mutate(
            time_diff = as.numeric(difftime(
                time,
                dplyr::lag(time),
                units = "hours"
            )),
            is_gap = dplyr::coalesce(time_diff > max_gap, FALSE),
            subseries_num = cumsum(is_gap) + 1L,
            Series_ID = paste(ID, subseries_num, sep = "_")
        ) |>
        dplyr::ungroup()

    # --- 3. Interpolate within each qualifying subseries ---
    df_interp <- df_series |>
        dplyr::group_by(ID, Series_ID) |>
        dplyr::filter(dplyr::n() >= min_points) |>
        dplyr::reframe(
            time_grid = list(sort(unique(c(
                seq(min(time), max(time), by = paste(step_mins, "mins")),
                max(time)
            )))),
            dplyr::across(
                dplyr::all_of(variable_cols),
                ~ list(
                    if (interpolation_method == "spline") {
                        .safe_spline(
                            x = as.numeric(time),
                            y = .x,
                            xout = as.numeric(time_grid[[1]]),
                            method = "natural",
                            min_pts = min_points
                        )
                    } else {
                        .safe_approx(
                            x = as.numeric(time),
                            y = .x,
                            xout = as.numeric(time_grid[[1]]),
                            min_pts = 2L
                        )
                    }
                )
            )
        ) |>
        tidyr::unnest(
            cols = c(time_grid, dplyr::all_of(variable_cols))
        ) |>
        dplyr::rename(time = time_grid) |>
        tidyr::drop_na(dplyr::all_of(variable_cols)) |>
        dplyr::ungroup()

    # --- 4. Compute summary statistics ---
    sampling_rate <- 60 / target_frequency
    n_participants <- length(unique(df_interp$ID))
    n_subseries <- length(unique(df_interp$Series_ID))
    pts_per_participant <- df_interp |>
        dplyr::count(ID) |>
        dplyr::pull("n")
    avg_pts <- if (length(pts_per_participant) > 0) {
        mean(pts_per_participant)
    } else {
        0L
    }

    # --- 5. Build ggplot ---
    p <- NULL
    if (plot) {
        valid_series <- df_interp |>
            dplyr::distinct(ID, Series_ID)

        # Long-format raw data restricted to usable subseries
        df_raw_long <- df_series |>
            dplyr::semi_join(valid_series, by = c("ID", "Series_ID")) |>
            dplyr::mutate(ID = as.factor(ID)) |>
            tidyr::pivot_longer(
                cols = dplyr::all_of(variable_cols),
                names_to = "Variable",
                values_to = "value"
            )

        # Plotting line: evaluate interpolation on (regular grid + raw timestamps)
        # so the curve passes through observed points and reaches the last one.
        df_plot_line <- df_series |>
            dplyr::semi_join(valid_series, by = c("ID", "Series_ID")) |>
            dplyr::group_by(ID, Series_ID) |>
            dplyr::reframe(
                time_plot = list(sort(unique(c(
                    seq(min(time), max(time), by = paste(step_mins, "mins")),
                    time
                )))),
                dplyr::across(
                    dplyr::all_of(variable_cols),
                    ~ list(
                        if (interpolation_method == "spline") {
                            .safe_spline(
                                x = as.numeric(time),
                                y = .x,
                                xout = as.numeric(time_plot[[1]]),
                                method = "natural",
                                min_pts = min_points
                            )
                        } else {
                            .safe_approx(
                                x = as.numeric(time),
                                y = .x,
                                xout = as.numeric(time_plot[[1]]),
                                min_pts = 2L
                            )
                        }
                    )
                )
            ) |>
            tidyr::unnest(
                cols = c(time_plot, dplyr::all_of(variable_cols))
            ) |>
            dplyr::rename(time = time_plot) |>
            tidyr::drop_na(dplyr::all_of(variable_cols)) |>
            dplyr::ungroup() |>
            dplyr::mutate(ID = as.factor(ID)) |>
            tidyr::pivot_longer(
                cols = dplyr::all_of(variable_cols),
                names_to = "Variable",
                values_to = "value"
            )

        # Offset: align each subseries/raw-series to hour 0
        if (plot_offset) {
            series_origins <- df_series |>
                dplyr::semi_join(valid_series, by = c("ID", "Series_ID")) |>
                dplyr::group_by(ID, Series_ID) |>
                dplyr::summarise(origin = min(time), .groups = "drop") |>
                dplyr::mutate(ID = as.factor(ID))

            df_plot_line <- df_plot_line |>
                dplyr::left_join(series_origins, by = c("ID", "Series_ID")) |>
                dplyr::mutate(
                    time = as.numeric(difftime(time, origin, units = "hours"))
                ) |>
                dplyr::select(-origin)

            df_raw_long <- df_raw_long |>
                dplyr::left_join(series_origins, by = c("ID", "Series_ID")) |>
                dplyr::mutate(
                    time = as.numeric(difftime(time, origin, units = "hours"))
                ) |>
                dplyr::select(-origin)
        }

        x_label <- if (plot_offset) "Hours elapsed" else "Time"

        p <- ggplot2::ggplot(
            mapping = ggplot2::aes(x = time, y = value, colour = ID)
        ) +
            ggplot2::geom_line(
                data = df_plot_line,
                mapping = ggplot2::aes(group = interaction(ID, Series_ID)),
                linewidth = 0.7
            ) +
            ggplot2::geom_point(
                data = df_raw_long,
                shape = 4,
                size = 2.5,
                stroke = 1.2
            ) +
            ggplot2::labs(
                x = x_label,
                y = "Value",
                colour = "Participant"
            ) +
            ggplot2::theme_minimal()

        if (length(variable_cols) > 1) {
            p <- p +
                ggplot2::facet_wrap(
                    ~Variable,
                    scales = "free_y",
                    ncol = 1
                )
        }

        if (n_participants == 1) {
            p <- p + ggplot2::guides(colour = "none")
        }
    }

    # --- 6. Optionally print a verbose summary ---
    if (verbose) {
        cat(sprintf("--- Preprocessing Complete ---\n"))
        cat(sprintf("Original data points    : %d\n", nrow(df)))
        cat(sprintf("Interpolated data points: %d\n", nrow(df_interp)))
        cat(sprintf("Usable subseries        : %d\n", n_subseries))
        cat(sprintf(
            "Output resolution       : 1 observation every %.4g minutes\n",
            step_mins
        ))
        cat(sprintf(
            "Sampling frequency      : %.4g observations per hour\n",
            sampling_rate
        ))
        cat(sprintf("------------------------------\n"))
    }

    # --- 7. Return experdyn_resampled object ---
    structure(
        list(
            resampled = df_interp,
            plot = p,
            interpolation_method = interpolation_method,
            target_frequency = target_frequency,
            sampling_rate = sampling_rate,
            max_gap = max_gap,
            n_participants = n_participants,
            n_subseries = n_subseries,
            n_points_raw = nrow(df),
            n_points_resampled = nrow(df_interp),
            avg_points_per_participant = avg_pts
        ),
        class = "experdyn_resampled"
    )
}

# =============================================================================
# Methods for experdyn_resampled
# =============================================================================

#' Print method for `experdyn_resampled` objects
#'
#' Displays a concise summary of a resampled ESM time series produced by
#' [resample_timeseries()].
#'
#' @param x An `experdyn_resampled` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.experdyn_resampled <- function(x, ...) {
    cat("========= Resampled Time Series =========\n")
    step_mins <- x$target_frequency
    cat(sprintf(
        "  Sampling rate     : %.4g obs/hour  (1 obs every %.4g min)\n",
        x$sampling_rate,
        step_mins
    ))
    cat(sprintf(
        "  Gap threshold     : %.4g hours\n",
        x$max_gap
    ))
    cat(sprintf(
        "  Participants      : %d\n",
        x$n_participants
    ))
    cat(sprintf(
        "  Subseries         : %d\n",
        x$n_subseries
    ))
    cat(sprintf(
        "  Raw obs           : %d\n",
        x$n_points_raw
    ))
    cat(sprintf(
        "  Resampled obs     : %d  (avg %.1f per participant)\n",
        x$n_points_resampled,
        x$avg_points_per_participant
    ))
    cat("  Data frame        : use `$resampled` to extract\n")
    cat("=========================================\n")
    invisible(x)
}
