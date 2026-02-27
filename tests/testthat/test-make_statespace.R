test_that("make_statespace returns experdyn_statespace with expected structure", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs, delay = 1, dimensions = 2)

    expect_s3_class(ss, "experdyn_statespace")
    expect_named(
        ss,
        c(
            "embedded",
            "delay",
            "dimensions",
            "variables",
            "n_series_input",
            "n_series_embedded"
        )
    )
    expect_equal(ss$delay, 1L)
    expect_equal(ss$dimensions, 2L)
    expect_s3_class(ss$embedded, "data.frame")
    expect_gt(nrow(ss$embedded), 0L)
    expect_true(all(c("ID", "Series_ID") %in% colnames(ss$embedded)))

    # lag columns should exist for each variable at each lag
    lag_cols <- grep("_lag\\d+$", colnames(ss$embedded), value = TRUE)
    expect_gte(length(lag_cols), 2L)
})

test_that("make_statespace embedded matrix has correct number of rows", {
    # A simple single-series case: n rows embed to n - (dim-1)*delay rows
    set.seed(1)
    t0 <- as.POSIXct("2025-01-01 09:00:00", tz = "UTC")
    times <- t0 + (0:19) * 3600
    df <- data.frame(
        timestamp = times,
        Component_1 = rnorm(20),
        stringsAsFactors = FALSE
    )
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        verbose = FALSE,
        plot = FALSE
    )
    n_resampled <- nrow(rs$resampled)
    ss <- make_statespace(rs, delay = 1, dimensions = 3)

    # Each series loses (dim - 1) * delay rows
    expect_equal(nrow(ss$embedded), n_resampled - (3 - 1) * 1)
})

test_that("make_statespace rejects non-positive delay/dimensions", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    expect_error(make_statespace(rs, delay = 0), "positive integer")
    expect_error(make_statespace(rs, dimensions = 0), "positive integer")
})

test_that("print.experdyn_statespace produces output", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    out <- capture.output(print(ss))
    expect_true(any(grepl("State Space", out)))
})
