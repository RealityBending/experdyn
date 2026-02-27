test_that("compute_rqa returns experdyn_rqa with expected structure", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs, delay = 1, dimensions = 2)
    rqa <- compute_rqa(ss, r = 1.0, verbose = FALSE)

    expect_s3_class(rqa, "experdyn_rqa")
    expect_named(
        rqa,
        c(
            "participant_rqa",
            "series_rqa",
            "recurrence_data",
            "r",
            "min_diag",
            "min_vert"
        )
    )
    expect_equal(rqa$r, 1.0)
    expect_equal(rqa$min_diag, 2L)
    expect_equal(rqa$min_vert, 2L)

    pr <- rqa$participant_rqa
    expect_s3_class(pr, "data.frame")
    expect_true(all(
        c("ID", "RR", "DET", "L", "Lmax", "ENTR", "LAM", "TT", "Vmax") %in%
            colnames(pr)
    ))
    expect_equal(nrow(pr), 2L)
})

test_that("compute_rqa: RR is in [0, 1]", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs, delay = 1, dimensions = 2)
    rqa <- compute_rqa(ss, r = 1.5, verbose = FALSE)

    rr_vals <- rqa$participant_rqa$RR
    expect_true(
        all(is.na(rr_vals) | (rr_vals >= 0 & rr_vals <= 1)),
        label = "RR must be in [0, 1]"
    )
})

test_that("compute_rqa: DET is in [0, 1]", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs, delay = 1, dimensions = 2)
    rqa <- compute_rqa(ss, r = 1.5, verbose = FALSE)

    det_vals <- rqa$participant_rqa$DET
    expect_true(
        all(is.na(det_vals) | (det_vals >= 0 & det_vals <= 1)),
        label = "DET must be in [0, 1]"
    )
})

test_that("compute_rqa: Lmax/Vmax are NA (not -Inf) when no lines qualify", {
    # Series is large enough for status="ok" but r is tiny so no diagonal/vertical
    # lines of length >= min_diag/min_vert form → Lmax and Vmax must be NA, not -Inf.
    set.seed(3)
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
    ss <- make_statespace(rs, delay = 1, dimensions = 2)
    # r = 0.001 → almost zero recurrences → no lines of length >= 8 qualify
    rqa <- compute_rqa(
        ss,
        r = 0.001,
        min_diag = 8L,
        min_vert = 8L,
        verbose = FALSE
    )

    pr <- rqa$participant_rqa
    expect_gt(
        nrow(pr),
        0L,
        label = "participant_rqa should have at least one row"
    )
    expect_false(
        any(pr$Lmax == -Inf, na.rm = TRUE),
        label = "Lmax should not be -Inf"
    )
    expect_false(
        any(pr$Vmax == -Inf, na.rm = TRUE),
        label = "Vmax should not be -Inf"
    )
})

test_that("compute_rqa: series_rqa has one row per subseries", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    rqa <- compute_rqa(ss, verbose = FALSE)

    expect_equal(
        nrow(rqa$series_rqa),
        length(unique(ss$embedded$Series_ID))
    )
})

test_that("compute_rqa: recurrence_data has correct columns", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    rqa <- compute_rqa(ss, r = 1.0, verbose = FALSE)

    rd <- rqa$recurrence_data
    if (nrow(rd) > 0) {
        expect_true(all(c("i", "j", "ID", "Series_ID") %in% colnames(rd)))
        expect_true(all(rd$i >= 1L))
        expect_true(all(rd$j >= 1L))
    }
})

test_that("compute_rqa rejects bad r", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    expect_error(compute_rqa(ss, r = -1), regexp = "positive")
})

test_that("print.experdyn_rqa produces output", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    rqa <- compute_rqa(ss, verbose = FALSE)
    out <- capture.output(print(rqa))
    expect_true(any(grepl("Recurrence", out)))
})

test_that("plot_rqa returns a ggplot", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    rqa <- compute_rqa(ss, r = 1.5, verbose = FALSE)

    expect_s3_class(plot_rqa(rqa), "gg")
    expect_s3_class(plot_rqa(rqa, facet_by = "series"), "gg")
})

test_that("plot_rqa: y-axis is not reversed (diagonal runs bottom-left to top-right)", {
    df <- make_esm_df()
    rs <- resample_timeseries(
        df,
        time_col = "timestamp",
        participant_col = "Participant",
        verbose = FALSE,
        plot = FALSE
    )
    ss <- make_statespace(rs)
    rqa <- compute_rqa(ss, r = 1.5, verbose = FALSE)
    p <- plot_rqa(rqa)

    # No 'reverse' transform on y — extracting scale should not carry it
    y_scale <- p$scales$get_scales("y")
    has_reverse <- !is.null(y_scale) &&
        !is.null(y_scale$trans) &&
        isTRUE(grepl("reverse", y_scale$trans$name))
    expect_false(
        has_reverse,
        label = "y-axis should not use the reverse transform"
    )
})

test_that("plot_rqa errors on non-experdyn_rqa input", {
    expect_error(plot_rqa(list()), regexp = "experdyn_rqa")
})
