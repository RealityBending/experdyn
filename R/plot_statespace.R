# =============================================================================
# Exported functions
# =============================================================================

#' Plot Embedded State Space in 2D
#'
#' @description
#' Projects an embedded multivariate state space to two dimensions and returns a
#' ggplot object.
#'
#' By default, the function attempts UMAP projection. If package `uwot` is not
#' installed, it automatically falls back to PCA.
#'
#' Embedded points are linked to visualise trajectories and, when `ID` is
#' present, results are shown in participant-level facets.
#'
#' @param data An `experdyn_statespace` object, or an embedded data frame
#'   (e.g., `x$embedded` from [make_statespace()]).
#' @param method `[character(1): "umap"]`\cr
#'   Projection method: `"umap"` or `"pca"`.
#' @param color_by `[character(1): "ID"]`\cr
#'   Column used for color mapping. If missing from the embedded data, points
#'   are shown in a single colour.
#' @param point_size `[numeric(1): 1.8]`\cr
#'   Point size in the output plot.
#' @param alpha `[numeric(1): 0.7]`\cr
#'   Point transparency in the output plot.
#' @param seed `[integer(1): 123]`\cr
#'   Random seed used for UMAP.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @export
plot_statespace <- function(
    data,
    method = "umap",
    color_by = "ID",
    point_size = 1.8,
    alpha = 0.7,
    seed = 123
) {
    # Resolve input
    embedded <- if (inherits(data, "experdyn_statespace")) {
        data$embedded
    } else {
        data
    }

    if (!is.data.frame(embedded)) {
        stop(
            "`data` must be an `experdyn_statespace` object or an embedded data frame.",
            call. = FALSE
        )
    }
    if (nrow(embedded) == 0) {
        stop("Embedded data is empty; nothing to plot.", call. = FALSE)
    }

    method <- match.arg(method, choices = c("umap", "pca"))

    metadata_cols <- c("ID", "Series_ID", "time", "row_index", ".row_index")
    lag_cols <- grep("_lag\\d+$", colnames(embedded), value = TRUE)

    if (length(lag_cols) == 0) {
        lag_cols <- setdiff(colnames(embedded), metadata_cols)
    }

    if (length(lag_cols) < 2) {
        stop(
            "Need at least 2 embedded dimensions/columns to plot a 2D projection.",
            call. = FALSE
        )
    }

    is_num <- vapply(embedded[lag_cols], is.numeric, logical(1))
    if (!all(is_num)) {
        bad_cols <- lag_cols[!is_num]
        stop(
            sprintf(
                "Non-numeric embedded columns found: %s",
                paste(bad_cols, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    mat <- as.matrix(embedded[, lag_cols, drop = FALSE])
    mat <- scale(mat)
    mat[is.na(mat)] <- 0

    method_used <- method
    if (method == "umap") {
        if (requireNamespace("uwot", quietly = TRUE)) {
            set.seed(seed)
            n_neighbors <- min(15L, max(2L, nrow(mat) - 1L))
            coords <- uwot::umap(
                mat,
                n_components = 2,
                n_neighbors = n_neighbors,
                metric = "euclidean",
                verbose = FALSE,
                ret_model = FALSE
            )
        } else {
            warning(
                "Package `uwot` not installed; falling back to PCA projection.",
                call. = FALSE
            )
            method_used <- "pca"
        }
    }

    if (method_used == "pca") {
        pca <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
        coords <- pca$x[, 1:2, drop = FALSE]
    }

    plot_df <- embedded
    plot_df$Dim1 <- coords[, 1]
    plot_df$Dim2 <- coords[, 2]

    if ("Series_ID" %in% colnames(plot_df)) {
        plot_df$.traj_group <- paste(
            plot_df$ID,
            plot_df$Series_ID,
            sep = "_"
        )
    } else if ("ID" %in% colnames(plot_df)) {
        plot_df$.traj_group <- as.character(plot_df$ID)
    } else {
        plot_df$.traj_group <- "all"
    }

    if (color_by %in% colnames(plot_df)) {
        plot_df$.color_group <- as.factor(plot_df[[color_by]])

        p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(x = Dim1, y = Dim2, colour = .color_group)
        )

        p <- p +
            ggplot2::geom_path(
                ggplot2::aes(group = .traj_group),
                alpha = min(alpha, 0.35),
                linewidth = 0.35
            )

        p <- p +
            ggplot2::geom_point(size = point_size, alpha = alpha) +
            ggplot2::labs(
                x = paste0(toupper(method_used), "-1"),
                y = paste0(toupper(method_used), "-2"),
                colour = color_by
            ) +
            ggplot2::theme_minimal()
    } else {
        p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(x = Dim1, y = Dim2)
        )

        p <- p +
            ggplot2::geom_path(
                ggplot2::aes(group = .traj_group),
                alpha = min(alpha, 0.35),
                linewidth = 0.35,
                colour = "#7F8C8D"
            )

        p <- p +
            ggplot2::geom_point(
                size = point_size,
                alpha = alpha,
                colour = "#2C3E50"
            ) +
            ggplot2::labs(
                x = paste0(toupper(method_used), "-1"),
                y = paste0(toupper(method_used), "-2")
            ) +
            ggplot2::theme_minimal()
    }

    if ("ID" %in% colnames(plot_df)) {
        p <- p + ggplot2::facet_wrap(~ID)
    }

    p
}
