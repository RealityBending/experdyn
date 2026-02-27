test_that("compute_mvsampen returns experdyn_mvsampen with expected structure", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs, delay = 1, dimensions = 2)
    mse <- compute_mvsampen(ss, r = 0.2)

    expect_s3_class(mse, "experdyn_mvsampen")
    expect_named(
        mse,
        c(
            "participant_entropy",
            "series_entropy",
            "dimensions",
            "r",
            "variables"
        )
    )
    expect_equal(mse$dimensions, 2L)
    expect_equal(mse$r, 0.2)

    pe <- mse$participant_entropy
    expect_s3_class(pe, "data.frame")
    expect_true(all(c("ID", "entropy", "entropy_sd_series") %in% colnames(pe)))
    expect_equal(nrow(pe), 2L) # two participants
})

test_that("compute_mvsampen: entropy_sd_series is numeric (not always NA)", {
    # Build a dataset with 2 well-separated subseries per participant.
    # Use 25 observations per burst so there are enough resampled + embedded
    # points, and r = 1.0 so template matching produces finite (non-Inf) entropy
    # in both bursts. n_series_valid therefore equals 2 for every participant.
    set.seed(7)
    pids <- paste0("P", 1:3)
    rows <- lapply(pids, function(pid) {
        t0 <- as.POSIXct("2025-01-01 09:00:00", tz = "UTC")
        # Burst 1: t0 … t0+24h; burst 2: t0+32h … t0+56h → gap = 8h > max_gap = 6h
        times <- c(
            t0 + (0:24) * 3600,
            t0 + 32 * 3600 + (0:24) * 3600
        )
        data.frame(
            Participant = pid,
            timestamp = times,
            Component_1 = cumsum(rnorm(50)),
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)

    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        max_gap = 6,
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs, delay = 1, dimensions = 2)
    mse <- compute_mvsampen(ss, r = 1.0)

    pe <- mse$participant_entropy

    # All participants must have 2 valid subseries in this dataset
    expect_true(
        all(pe$n_series_valid >= 2),
        label = "all participants should have at least 2 valid series"
    )
    # With >= 2 valid series the SD is always computable, never NaN
    expect_true(
        all(!is.nan(pe$entropy_sd_series)),
        label = "entropy_sd_series should not be NaN for participants with 2+ valid series"
    )
})

test_that("compute_mvsampen: entropy_sd_series is NA for single-series participants", {
    # Single burst per participant → only one series → SD is undefined
    df <- make_esm_df(n = 20, n_participants = 2)
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs, delay = 1, dimensions = 2)
    mse <- compute_mvsampen(ss, r = 0.2)

    # Each participant has a single subseries; sd of one value must be NA
    expect_true(
        all(is.na(mse$participant_entropy$entropy_sd_series)),
        label = "entropy_sd_series should be NA when only one series per participant"
    )
})

test_that("compute_mvsampen: series_entropy has one row per subseries", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    mse <- compute_mvsampen(ss)

    expect_equal(
        nrow(mse$series_entropy),
        length(unique(ss$embedded$Series_ID))
    )
})

test_that("print.experdyn_mvsampen produces output", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    mse <- compute_mvsampen(ss)
    out <- capture.output(print(mse))
    expect_true(any(grepl("Sample Entropy", out)))
})

test_that("compute_mvsampen rejects non-statespace input", {
    df <- make_esm_df()
    expect_error(
        compute_mvsampen(df),
        regexp = "experdyn_statespace"
    )
})
