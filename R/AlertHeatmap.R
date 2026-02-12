#==============================================
# NOTE AlertSystemScore Plotting - AlertHeatmap
#==============================================

#' AlertHeatmap: visualize pathway-level perturbation metrics as a heatmap
#'
#' This function visualizes the AlertSystem metrics stored in
#' \code{seurat_obj@misc$AlertSystem$metrics} as a pathways-by-group heatmap.
#' It selects the top N most positively and/or negatively changed pathways
#' based on the mean \code{effect} (summarized logFC) across groups, and can
#' optionally rescale values for visualization.
#'
#' Tiles are filled by the \code{effect} column, while additional metrics
#' such as \code{concordance} and the robustness score
#' (\code{robust_gt_thresh}, fraction of strong, sign-consistent cells among
#' all cells in a group) are displayed as text overlays (labels) if present
#' and non-missing. In the plot, this robustness is shown simply as
#' \code{"robustness"}.
#'
#' @importFrom rlang .data
#'
#' @param seurat_obj A Seurat object that has been processed with
#'   \code{AlertSystemScore()} and contains a metrics table in
#'   \code{seurat_obj@misc$AlertSystem$metrics}. This can be either from
#'   a fresh run of \code{AlertSystemScore()} or after re-attaching a
#'   previously saved AlertSystem state (e.g. via \code{LoadAlertSystem()}
#'   or \code{readRDS("..._AlertSystem_info.rds")}).
#' @param misc_name Name of the misc slot where AlertSystem results are stored.
#'   Default is \code{"AlertSystem"}.
#' @param top_n Integer; number of top pathways to select per side (positive
#'   and/or negative). For example, \code{top_n = 5} and \code{side = "both"}
#'   can yield up to 10 pathways (5 most positive, 5 most negative).
#' @param side Character string indicating which side(s) of the effect
#'   distribution to use. One of \code{"both"}, \code{"positive"}, or
#'   \code{"negative"}. Default is \code{"both"}.
#' @param scale_mode Character; one of:
#'   \describe{
#'     \item{"none"}{Use raw effect values (default).}
#'     \item{"unit"}{Scale effect into approximately \code{c(-1, 1)} by dividing by a
#'       global magnitude (see \code{unit_mode}, \code{unit_q}).}
#'     \item{"zscale"}{Z-score effect values within each pathway across groups:
#'       \eqn{z = (x - \mu)/\sigma}. If there is only one group, this mode
#'       automatically falls back to \code{"none"}.}
#'   }
#' @param unit_mode Character; only used when \code{scale_mode = "unit"}.
#'   One of:
#'   \describe{
#'     \item{"max"}{Divide by the maximum absolute effect across all selected
#'       pathway/group combinations. This makes the most extreme value exactly
#'       ±1.}
#'     \item{"robust"}{Divide by a high quantile (see \code{unit_q}) of the
#'       absolute effect distribution, which reduces the influence of extreme
#'       outliers and improves contrast among moderately changed pathways.}
#'   }
#' @param select_by Character; metric used to rank pathways for inclusion
#'   within each side of the mean \code{effect} distribution. One of
#'   \code{"concordance"} (default), \code{"effect"}, or \code{"abs_effect"}.
#' @param reorder_by Character; if \code{reorder = TRUE}, controls how pathways
#'   are ordered on the y-axis. One of \code{"effect"} (default; negative to
#'   positive mean effect) or \code{"select"} (order by the selection rank).
#' @param unit_q Numeric between 0 and 1; quantile of \code{|effect|} used as
#'   the scaling factor when \code{unit_mode = "robust"}. Default is 0.9.
#' @param reorder Logical; if \code{TRUE} (default), reorder pathways on the
#'   y-axis by their mean effect (from most negative to most positive).
#' @param fill_limits Optional numeric vector of length 2 giving limits for
#'   the fill scale. If \code{NULL}, limits are set automatically depending
#'   on \code{scale_mode}.
#' @param show_concordance Logical; if TRUE (default), show the concordance
#'   metric as text labels when available.
#' @param show_robustness Logical; if TRUE, show the robustness metric
#'   (column \code{robust_gt_thresh}) as text labels when available. Default
#'   is FALSE, so robustness is not shown even if computed, unless
#'   explicitly requested.
#'
#' @return A ggplot2 heatmap object.
#'
#' @export

AlertHeatmap <- function(
  seurat_obj,
  misc_name   = "AlertSystem",
  top_n       = 5,
  side        = c("both", "positive", "negative"),
  scale_mode  = c("none", "unit", "zscale"),
  unit_mode   = c("max", "robust"),
  select_by  = c("concordance", "effect", "abs_effect"),
  reorder_by = c("effect", "select"),
  unit_q      = 0.9,
  reorder     = TRUE,
  fill_limits = NULL,
  show_concordance   = TRUE,
  show_robustness = FALSE
) {

  side       <- match.arg(side)
  scale_mode <- match.arg(scale_mode)
  unit_mode  <- match.arg(unit_mode)
  select_by  <- match.arg(select_by)
  reorder_by <- match.arg(reorder_by)


  ## -------------------------
  ## Fetch metrics
  ## -------------------------
  if (is.null(seurat_obj@misc[[misc_name]]) ||
      is.null(seurat_obj@misc[[misc_name]][["metrics"]])) {
    stop(
      "No metrics table found in seurat_obj@misc[['", misc_name,
      "']][['metrics']]. Run AlertSystemScore() first."
    )
  }

  metrics_df <- seurat_obj@misc[[misc_name]][["metrics"]]

  required_cols <- c("group", "pathway", "effect")
  missing_cols  <- setdiff(required_cols, colnames(metrics_df))
  if (length(missing_cols) > 0L) {
    stop(
      "metrics table must contain columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ", paste(missing_cols, collapse = ", ")
    )
  }

  # Original group label for x-axis: try to recover group.by from params
  group_col_orig <- "group"
  if (!is.null(seurat_obj@misc[[misc_name]][["params"]])) {
    params <- seurat_obj@misc[[misc_name]][["params"]]
    if (!is.null(params$group.by) && !is.na(params$group.by)) {
      group_col_orig <- params$group.by
    }
  }

  # Copy and ensure basic types
  df_long <- metrics_df
  df_long$group   <- as.factor(df_long$group)
  df_long$pathway <- as.character(df_long$pathway)

  ## -------------------------
  ## Pathway selection by mean effect
  ## -------------------------

  # path_stats <- df_long |>
  #   dplyr::group_by(.data$pathway) |>
  #   dplyr::summarise(
  #     mean_effect = mean(.data$effect, na.rm = TRUE),
  #     .groups     = "drop"
  #   )

  has_concord_global <- "concordance" %in% colnames(df_long) &&
  any(!is.na(df_long$concordance))

  path_stats <- df_long |>
    dplyr::group_by(.data$pathway) |>
    dplyr::summarise(
      mean_effect      = mean(.data$effect, na.rm = TRUE),
      mean_abs_effect  = mean(abs(.data$effect), na.rm = TRUE),
      mean_concordance = if (has_concord_global) mean(.data$concordance, na.rm = TRUE) else NA_real_,
      .groups          = "drop"
    )

  # Select “top up / top down” by concordance (default), within effect sign
  # Define a ranking column that’s safe if concordance is missing:
  rank_col <- dplyr::case_when(
    select_by == "concordance" ~ ifelse(has_concord_global, path_stats$mean_concordance, path_stats$mean_abs_effect),
    select_by == "abs_effect"  ~ path_stats$mean_abs_effect,
    TRUE                       ~ abs(path_stats$mean_effect)  # effect: rank by magnitude but keep side split by sign
  )

  path_stats$rank_value <- rank_col


  # Positive and negative sets based on mean_effect
  # pos_paths <- path_stats |>
  #   dplyr::filter(.data$mean_effect > 0) |>
  #   dplyr::arrange(dplyr::desc(.data$mean_effect)) |>
  #   dplyr::slice_head(n = top_n) |>
  #   dplyr::pull(.data$pathway)
  #
  # neg_paths <- path_stats |>
  #   dplyr::filter(.data$mean_effect < 0) |>
  #   dplyr::arrange(.data$mean_effect) |>
  #   dplyr::slice_head(n = top_n) |>
  #   dplyr::pull(.data$pathway)

  pos_paths <- path_stats |>
    dplyr::filter(.data$mean_effect > 0) |>
    dplyr::arrange(dplyr::desc(.data$rank_value)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(.data$pathway)

  neg_paths <- path_stats |>
    dplyr::filter(.data$mean_effect < 0) |>
    dplyr::arrange(dplyr::desc(.data$rank_value)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(.data$pathway)


  # Choose pathways based on "side"
  selected_paths <- character(0)

  if (side %in% c("both", "positive")) {
    if (length(pos_paths) == 0L && side == "both") {
      message("AlertHeatmap: No positively shifted pathways detected (mean effect > 0).")
    }
    selected_paths <- c(selected_paths, pos_paths)
  }

  if (side %in% c("both", "negative")) {
    if (length(neg_paths) == 0L && side == "both") {
      message("AlertHeatmap: No negatively shifted pathways detected (mean effect < 0).")
    }
    selected_paths <- c(selected_paths, neg_paths)
  }

  selected_paths <- unique(selected_paths)

  if (length(selected_paths) == 0L) {
    stop(
      "AlertHeatmap: No pathways selected for side = '", side,
      "'. Check effect values or reduce top_n."
    )
  }

  # Filter to selected pathways only
  df_long <- df_long[df_long$pathway %in% selected_paths, , drop = FALSE]

  # # Optionally reorder pathways: from most negative to most positive mean_effect
  # if (reorder) {
  #   order_df <- path_stats[path_stats$pathway %in% selected_paths, , drop = FALSE]
  #   order_df <- order_df[order(order_df$mean_effect), , drop = FALSE]
  #   df_long$pathway <- factor(df_long$pathway, levels = order_df$pathway)
  # } else {
  #   df_long$pathway <- as.factor(df_long$pathway)
  # }
  if (reorder) {
  order_df <- path_stats[path_stats$pathway %in% selected_paths, , drop = FALSE]

  if (reorder_by == "select") {
    order_df <- order_df[order(order_df$rank_value, decreasing = TRUE), , drop = FALSE]
  } else {
    order_df <- order_df[order(order_df$mean_effect), , drop = FALSE] # neg -> pos
  }

    df_long$pathway <- factor(df_long$pathway, levels = order_df$pathway)
  } else {
    df_long$pathway <- as.factor(df_long$pathway)
  }


  # Number of groups (for zscale sanity)
  n_groups <- dplyr::n_distinct(df_long$group)

  # We will plot using a separate column (plot_value) derived from effect,
  # so the original effect remains untouched in df_long.
  df_long$plot_value <- df_long$effect

  ## -------------------------
  ## Decide what drives the heatmap fill
  ## -------------------------
  legend_title <- "effect (logFC)"

  # if (fill_by == "concordance") {
  # if (!("concordance" %in% colnames(df_long))) {
  #   stop("fill_by = 'concordance' requested but no 'concordance' column is present in metrics.")
  # }
  #   df_long$plot_value <- df_long$concordance
  #   legend_title <- "concordance"
  #
  #   ## IMPORTANT GUARD:
  #   ## Concordance is already bounded [0,1], so scaling makes no sense
  #   scale_mode <- "none"
  #   if (is.null(fill_limits)) {
  #     fill_limits <- c(0, 1)
  #   }
  #
  # } else {
  #   df_long$plot_value <- df_long$effect
  #   legend_title <- "effect (logFC)"
  # }


  if (scale_mode == "unit") {
    legend_title <- "scaled effect"

    if (unit_mode == "max") {
      # scale by the absolute maximum
      scale_val <- max(abs(df_long$plot_value), na.rm = TRUE)
    } else {
      # robust scaling: high quantile of |effect|
      scale_val <- stats::quantile(
        abs(df_long$plot_value),
        probs = unit_q,
        na.rm  = TRUE
      )
    }

    if (is.finite(scale_val) && scale_val > 0) {
      # scale
      df_long$plot_value <- df_long$plot_value / scale_val

      # clamp to [-1, 1] so values outside limits don't become NA (grey)
      df_long$plot_value[df_long$plot_value >  1] <-  1
      df_long$plot_value[df_long$plot_value < -1] <- -1
    } else {
      warning("AlertHeatmap: unit scaling skipped (non-finite scale).")
    }

  } else if (scale_mode == "zscale") {
    if (n_groups < 2) {
      warning(
        "AlertHeatmap: scale_mode = 'zscale' requested but only one group present; ",
        "using scale_mode = 'none' instead."
      )
      scale_mode <- "none"
    } else {
      legend_title <- "z-scaled effect"

      df_long <- df_long |>
        dplyr::group_by(.data$pathway) |>
        dplyr::mutate(
          .mu = mean(.data$plot_value, na.rm = TRUE),
          .sd = stats::sd(.data$plot_value, na.rm = TRUE),
          plot_value = ifelse(is.finite(.sd) & .sd > 0,
                              (plot_value - .mu) / .sd,
                              0)
        ) |>
        dplyr::select(-.mu, -.sd) |>
        dplyr::ungroup()
    }
  }

  ## -------------------------
  ## Fill limits
  ## -------------------------
  if (is.null(fill_limits)) {
    if (scale_mode == "zscale") {
      max_abs <- max(abs(df_long$plot_value), na.rm = TRUE)
      if (!is.finite(max_abs) || max_abs == 0) {
        max_abs <- 1
      }
      fill_limits <- c(-max_abs, max_abs)

    } else if (scale_mode == "unit") {
      min_v <- min(df_long$plot_value, na.rm = TRUE)
      max_v <- max(df_long$plot_value, na.rm = TRUE)
      if (min_v >= 0) {
        fill_limits <- c(0, 1)      # only positive
      } else if (max_v <= 0) {
        fill_limits <- c(-1, 0)     # only negative
      } else {
        fill_limits <- c(-1, 1)     # both sides

      }

    } else { # "none"
      max_abs <- max(abs(df_long$plot_value), na.rm = TRUE)
      if (!is.finite(max_abs) || max_abs == 0) {
        max_abs <- 1
      }
      fill_limits <- c(-max_abs, max_abs)
    }
  }

  # ## -------------------------
  # ## Optional text annotations from metrics
  # ## -------------------------
  # has_concord <- "concordance" %in% colnames(df_long) &&
  #   any(!is.na(df_long$concordance))
  #
  # # has_cred <- "robust_gt_thresh" %in% colnames(df_long) &&
  # #   any(!is.na(df_long$robust_gt_thresh))
  # has_cred <- ("robust_gt_thresh" %in% colnames(df_long)) &&
  #             any(!is.na(df_long$robust_gt_thresh))
  #
  # df_long$label <- NA_character_
  #
  # if (has_concord && has_cred) {
  #   df_long$label <- sprintf("concordance=%.2f\nrobustness=%.2f",
  #                            df_long$concordance,
  #                            df_long$robust_gt_thresh)
  # } else if (has_concord && !has_cred) {
  #   df_long$label <- sprintf("concordance=%.2f", df_long$concordance)
  # } else if (!has_concord && has_cred) {
  #   df_long$label <- sprintf("robustness=%.2f", df_long$robust_gt_thresh)
  # }
  #
  # df_long$label <- NA_character_
  #
  # if (has_concord || has_cred) {
  #   # define a simple "cred" column from robust_gt_thresh
  #   df_long$cred <- if (has_cred) df_long$robust_gt_thresh else NA_real_
  #
  #   if (has_concord && has_cred) {
  #     df_long$label <- sprintf(
  #       "concordance=%.2f\nrobustness=%.2f",
  #       df_long$concordance,
  #       df_long$cred
  #     )
  #   } else if (has_concord && !has_cred) {
  #     df_long$label <- sprintf("concordance=%.2f", df_long$concordance)
  #   } else if (!has_concord && has_cred) {
  #     df_long$label <- sprintf("robustness=%.2f", df_long$cred)
  #   }
  # }

  ## -------------------------
  ## Optional text annotations from metrics
  ## -------------------------
  has_concord_col <- "concordance" %in% colnames(df_long) &&
    any(!is.na(df_long$concordance))

  has_robust_col <- "robust_gt_thresh" %in% colnames(df_long) &&
    any(!is.na(df_long$robust_gt_thresh))

  # # What we actually show depends on both presence *and* user flags
  # use_concord <- show_concordance   && has_concord_col
  # use_cred <- show_robustness && has_robust_col
  #
  # df_long$label <- NA_character_
  # df_long$cred  <- NA_real_

  # if (use_concord || use_cred) {
  #   if (use_cred) {
  #     df_long$cred <- df_long$robust_gt_thresh
  #   }
  #
  #   if (use_concord && use_cred) {
  #     df_long$label <- sprintf(
  #       "concordance=%.2f\nrobustness=%.2f",
  #       df_long$concordance,
  #       df_long$cred
  #     )
  #   } else if (use_concord && !use_cred) {
  #     df_long$label <- sprintf("concordance=%.2f", df_long$concordance)
  #   } else if (!use_concord && use_cred) {
  #     df_long$label <- sprintf("robustness=%.2f", df_long$cred)
  #   }
  # }

  # What we actually show depends on both presence *and* user flags
  use_concord   <- show_concordance   && has_concord_col
  use_robust <- show_robustness && has_robust_col

  df_long$label  <- NA_character_
  df_long$robust <- NA_real_

  if (use_concord || use_robust) {

    if (use_robust) {
      df_long$robust <- df_long$robust_gt_thresh
    }

    if (use_concord && use_robust) {
      df_long$label <- sprintf(
        "concordance=%.2f\nrobustness=%.2f",
        df_long$concordance,
        df_long$robust
      )
    } else if (use_concord && !use_robust) {
      df_long$label <- sprintf("concordance=%.2f", df_long$concordance)
    } else if (!use_concord && use_robust) {
      df_long$label <- sprintf("robustness=%.2f", df_long$robust)
    }
  }


  ## -------------------------
  ## Dynamic title
  ## -------------------------
  plot_title <- dplyr::case_when(
    side == "both"     ~ "Top Positive & Negative Pathways",
    side == "positive" ~ "Top Positively Shifted Pathways",
    side == "negative" ~ "Top Negatively Shifted Pathways",
    TRUE ~ "AlertSystem Pathway Summary"
  )

  ## -------------------------
  ## Build plot
  ## -------------------------
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = group, y = pathway, fill = plot_value)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low       = "#1682C0",
      mid       = "#FFF5CD",
      high      = "#A02248",
      midpoint  = 0,
      name      = legend_title,
      limits    = fill_limits,
      na.value  = "grey90",
      oob       = scales::squish
    ) +
    ggplot2::labs(
      x     = group_col_orig,
      y     = "Pathway",
      title = plot_title
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title  = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      panel.grid  = ggplot2::element_blank()
    )

  # Add labels if we actually have any non-NA annotations
  if (!all(is.na(df_long$label))) {
    p <- p +
      ggplot2::geom_text(
        data = df_long,
        ggplot2::aes(
          label = label,
          color = ifelse(plot_value > 0.1, "white",
                  ifelse(plot_value < -0.1, "white", "black"))
        ),
        size        = 2.6,
        lineheight  = 0.8,
        fontface    = "plain",   # ← not bold
        show.legend = FALSE,
        na.rm       = TRUE
      ) +
      ggplot2::scale_color_identity()
  }

  p
}
