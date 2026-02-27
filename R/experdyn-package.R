#' experdyn: Experiential Dynamics for ESM Time Series
#'
#' @description
#' **experdyn** provides a growing toolkit for applying nonlinear dynamics and
#' complexity science to behavioral time series from Experience Sampling Method
#' (ESM) / Ecological Momentary Assessment (EMA) studies.
#'
#' ## What is an ESM study?
#' Participants receive repeated notifications throughout the day (e.g., every
#' 1-3 hours) and respond to short questionnaires about their current mood,
#' activity, or psychological state. The resulting data are inherently
#' **irregularly sampled** multivariate time series — one series per participant.
#'
#' ## What does experdyn offer?
#' The package is organised around three main themes:
#'
#' \describe{
#'   \item{Preprocessing}{Regularise and clean raw ESM time series before
#'     analysis. Key function: `resample_timeseries()`.}
#'   \item{State-space methods}{Reconstruct the attractor of a dynamical system
#'     via delay-based embedding (Takens' theorem).}
#'   \item{Complexity measures}{Compute entropy, fractal dimensions, recurrence
#'     quantification analysis, and related indices that capture the nonlinear
#'     temporal structure of psychological dynamics.}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## All dplyr, ggplot2, tidyr, tidyselect, and nonlinearTseries calls use
## explicit :: notation throughout the package, so no NAMESPACE imports are
## needed beyond listing the packages in DESCRIPTION Imports.
## usethis namespace: end
NULL

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        ".color_group",
        ".traj_group",
        "Dim1",
        "Dim2",
        "Estimate",
        "i",
        "ID",
        "j",
        "n",
        "n_obs",
        "Parameter",
        "Series_ID",
        "Variable",
        "delay",
        "delay_iqr",
        "delay_mean",
        "delay_median",
        "delay_sd",
        "dim",
        "dim_iqr",
        "dim_mean",
        "dim_median",
        "dim_sd",
        "DET",
        "ENTR",
        "L",
        "LAM",
        "Lmax",
        "RR",
        "TT",
        "Vmax",
        "entropy",
        "is_gap",
        "n_series",
        "pairs_m",
        "origin",
        "row_index",
        "status",
        "subseries_num",
        "time",
        "time_diff",
        "time_grid",
        "time_plot",
        "value"
    ))
}
