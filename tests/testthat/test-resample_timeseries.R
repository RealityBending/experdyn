test_that("resample_timeseries returns experdyn_resampled with expected structure", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )

    expect_s3_class(rs, "experdyn_resampled")
    expect_named(
        rs,
        c(
            "resampled",
            "plot",
            "interpolation_method",
            "target_frequency",
            "sampling_rate",
            "max_gap",
            "n_participants",
            "n_subseries",
            "n_points_raw",
            "n_points_resampled",
            "avg_points_per_participant"
        )
    )

    expect_s3_class(rs$resampled, "data.frame")
    expect_true(all(c("ID", "Series_ID", "time") %in% colnames(rs$resampled)))
    expect_true(any(grepl("^Component_", colnames(rs$resampled))))

    expect_equal(rs$n_participants, 2L)
    expect_gt(rs$n_points_resampled, 0L)
    expect_null(rs$plot)
})

test_that("resample_timeseries plot is a ggplot when plot = TRUE", {
    df <- make_esm_df(n = 20, n_participants = 1)
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = TRUE
    )
    expect_s3_class(rs$plot, "gg")
})

test_that("resample_timeseries works with a single participant (participant_col = NULL)", {
    df <- make_esm_df(n = 20, n_participants = 1)
    df$Participant <- NULL
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        verbose = FALSE,
        plot = FALSE
    )
    expect_s3_class(rs, "experdyn_resampled")
    expect_equal(rs$n_participants, 1L)
    expect_true(all(rs$resampled$ID == "All"))
})

test_that("resample_timeseries match.arg rejects bad interpolation_method", {
    df <- make_esm_df()
    expect_error(
        resample_timeseries(
            df,
            time_col = "timestamp",
            participant_col = "Participant",
            interpolation_method = "loess",
            verbose = FALSE,
            plot = FALSE
        ),
        regexp = "should be one of"
    )
})

test_that("resample_timeseries accepts partial match for interpolation_method", {
    df <- make_esm_df()
    expect_no_error(
        resample_timeseries(
            df,
            time_col = "timestamp",
            participant_col = "Participant",
            interpolation_method = "lin", # partial match → "linear"
            verbose = FALSE,
            plot = FALSE
        )
    )
})

test_that("resample_timeseries splits across large gaps", {
    set.seed(1)
    t0 <- as.POSIXct("2025-01-01 09:00:00", tz = "UTC")
    # Two bursts: burst 1 ends at t0+4h, burst 2 starts at t0+12h → gap = 8h > max_gap = 6h
    times <- c(
        t0 + (0:4) * 3600,
        t0 + 12 * 3600 + (0:4) * 3600
    )
    df <- data.frame(
        timestamp = times,
        Component_1 = rnorm(10),
        stringsAsFactors = FALSE
    )
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        max_gap = 6,
        verbose = FALSE,
        plot = FALSE
    )
    expect_gte(rs$n_subseries, 2L)
})

test_that("print.experdyn_resampled produces output", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    out <- capture.output(print(rs))
    expect_true(any(grepl("Resampled", out)))
})
