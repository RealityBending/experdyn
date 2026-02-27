# Shared minimal dataset used across several test files.
# Sourced implicitly by testthat via helper-* convention.

make_esm_df <- function(n = 20, n_participants = 2, seed = 42) {
    set.seed(seed)
    pids <- paste0("P", seq_len(n_participants))
    rows <- lapply(pids, function(pid) {
        t0 <- as.POSIXct("2025-01-01 09:00:00", tz = "UTC")
        times <- t0 + cumsum(sample(c(3600, 5400, 7200), n, replace = TRUE))
        data.frame(
            Participant = pid,
            timestamp = times,
            Component_1 = cumsum(rnorm(n)),
            Component_2 = cumsum(rnorm(n)),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}
