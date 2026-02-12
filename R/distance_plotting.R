# Description: Plotting Computed Distance for each group in a Seurat object

# # Load necessary libraries
# library(ggplot2)
# library(reshape2)
# library(patchwork)  # For arranging multiple plots

# Internal: Get the Upper Triangle of a Matrix
#'
#' This function is used internally to extract the upper triangle of a matrix.
#' It replaces the lower triangle of the input matrix with `NA`, which is useful
#' when working with symmetrical matrices such as distance or correlation matrices
#' where only the upper triangle is needed.
#'
#' @param cormat A numeric matrix (e.g., a correlation or distance matrix).
#' @return A matrix with `NA` values in the lower triangle.
#' @note This function is for internal use and is not exported.
.get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

# Internal: Create a Heatmap from a Matrix
#'
#' This function generates a heatmap from a given matrix. It handles both
#' matrices and data frames (converting them to matrices if needed). It also allows
#' customization of the color palette, axis visibility, and legend display.
#'
#' @param df_matrix A numeric matrix or a data frame to be converted to a matrix.
#' @param title The title for the heatmap.
#' @param min_val The minimum value for the color scale.
#' @param max_val The maximum value for the color scale.
#' @param custom_palette A vector of colors to define the color palette. If `NULL`, a default palette is used.
#' @param show_x_axis Logical; whether to display the x-axis labels. Defaults to `TRUE`.
#' @param show_y_axis Logical; whether to display the y-axis labels. Defaults to `TRUE`.
#' @param show_legend Logical; whether to display the legend. Defaults to `TRUE`.
#' @param custom_order A character vector specifying the order of clusters for the x and y axes.
#' @return A ggplot2 object representing the heatmap.
#' @note This function is for internal use and is not exported.
.create_distance_heatmap <- function(df_matrix, title, min_val, max_val,
                                     custom_palette = NULL,
                                     show_x_axis = TRUE, show_y_axis = TRUE,
                                     show_legend = TRUE,
                                     custom_order = NULL) {

  # Define the color palette, use custom_palette if provided, otherwise use the default palette
  if (is.null(custom_palette)) {
    custom_palette <- c("#F9F3E1", "#F38B60", "#AF3B3B", "#2D1E3E")  # Default palette
  }

  # Check if the input is a dataframe and convert to matrix if needed
  if (is.data.frame(df_matrix)) {
    message("Converting input dataframe to matrix for: '", title, "'")
    df_matrix <- as.matrix(df_matrix)
  } else if (is.matrix(df_matrix)) {
    message("Input is a data matrix: '", title, "'")
  } else {
    stop("Error: The input must be either a data frame or a matrix.")
  }

  if (!isTRUE(all.equal(df_matrix, t(df_matrix), tolerance = 1e-8))) {
    warning("Matrix is not symmetric; plotting upper triangle only.")
  }


  # Get the upper triangle of the matrix
  upper_tri <- .get_upper_tri(df_matrix)

  # Melt the matrix into long format
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

  # Apply custom ordering to the variables on the x and y axes if provided
  if (!is.null(custom_order)) {
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = custom_order)
    melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = custom_order)
  }

  # Create the heatmap with the specified color palette
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = custom_palette, limits = c(min_val, max_val), na.value = "transparent", name = "Distance", guide = if (show_legend) "colourbar" else "none") +
    theme_minimal() +
    theme(axis.text.x = if (show_x_axis) element_text(angle = 45, vjust = 1, size = 12, hjust = 1) else element_blank(),
          axis.text.y = if (show_y_axis) element_text(size = 12) else element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          panel.border = element_blank()) +
    coord_fixed() +
    ggtitle(title)

  return(ggheatmap)
}



#' @title HeatmapDistance: Generate Heatmaps for Original and Perturbed Matrices
#'
#' This function generates two heatmaps from two matrices (original and perturbed)
#' and displays them side by side on the same color scale.
#'
#' @param df_original A numeric matrix representing the original (unperturbed) data.
#' @param df_perturbed A numeric matrix representing the perturbed data.
#' @param custom_palette A vector of colors to define the color palette. Defaults to a red/blue gradient.
#' @param title_original The title for the original heatmap. Defaults to "Original Assay Cluster Similarity Distance".
#' @param title_perturbed The title for the perturbed heatmap. Defaults to "Perturbed Assay Cluster Similarity Distance".
#' @param custom_order A character vector specifying the order of clusters for the x and y axes.
#' @param col_fun Optional color mapping function created by \code{circlize::colorRamp2()}.
#'   If provided, \code{HeatmapDistance()} will derive a shared color scale from the
#'   \code{breaks} and \code{colors} attributes of \code{col_fun} and apply it to both heatmaps.
#'   When \code{col_fun} is supplied, \code{custom_palette} is ignored.
#'
#' @return A patchwork object combining the two heatmaps.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#'   theme_minimal theme element_text element_blank coord_fixed ggtitle
#' @importFrom reshape2 melt
#' @importFrom patchwork wrap_plots
#'
#' @export
#' @examples
#' p <- HeatmapDistance(df_edist_observed, df_edist_perturbed) # , custom_order = custom_order
HeatmapDistance <- function(df_original, df_perturbed, custom_palette = NULL,
                            title_original = "Original Assay Cluster Similarity Distance",
                            title_perturbed = "Perturbed Assay Cluster Similarity Distance",
                            custom_order = NULL,
                            col_fun = NULL) {
  #
  if (is.data.frame(df_original)) df_original <- as.matrix(df_original)
  if (is.data.frame(df_perturbed)) df_perturbed <- as.matrix(df_perturbed)

  # Check if the dimensions of the matrices match
  if (!all(dim(df_original) == dim(df_perturbed))) {
    stop("Error: 'df_original' and 'df_perturbed' must have the same dimensions.")
  } else {
    message("Confirmed: The dimensions of 'df_original' and 'df_perturbed' match.")
  }

  # if custom_order contains a name not present in the matrix, you’ll get a subscript error
  # ---- VALIDATE custom_order (ADD THIS HERE) ----
  if (!is.null(custom_order)) {
    mat <- as.matrix(df_original)

    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop("df_original must have rownames and colnames for custom_order to work.")
    }

    # Check row labels
    missing <- setdiff(custom_order, rownames(mat))
    if (length(missing) > 0) {
      stop(
        "custom_order has labels not found in matrix rows: ",
        paste(missing, collapse = ", ")
      )
    }

    # Check column labels
    missing2 <- setdiff(custom_order, colnames(mat))
    if (length(missing2) > 0) {
      stop(
        "custom_order has labels not found in matrix columns: ",
        paste(missing2, collapse = ", ")
      )
    }

    # Warn if rows/columns are inconsistent
    if (!identical(rownames(mat), colnames(mat))) {
      warning(
        "df_original rownames and colnames differ; ",
        "ensure matrices are labeled consistently."
      )
    }
  }
  # ------------------------------


  # Reorder matrices if custom_order is provided
  if (!is.null(custom_order)) {
    df_original <- df_original[custom_order, custom_order]
    df_perturbed <- df_perturbed[custom_order, custom_order]
  }

  # ---- NEW: if col_fun is provided, derive palette + range from it ----
  # This makes your existing .create_distance_heatmap() work without edits,
  # as long as it uses scale_fill_gradientn(colors=custom_palette, limits=c(min,max)).
  if (!is.null(col_fun)) {
    # Accept either circlize::colorRamp2 function or a list with breaks/colors
    brks <- attr(col_fun, "breaks")
    cols <- attr(col_fun, "colors")

    if (is.null(brks) || is.null(cols)) {
      stop("col_fun must be created by circlize::colorRamp2() (must have 'breaks' and 'colors' attributes).")
    }

    # circlize stores colors with alpha (#RRGGBBAA). ggplot is fine with it, but we can drop AA safely.
    cols <- gsub("^#([0-9A-Fa-f]{6})([0-9A-Fa-f]{2})$", "#\\1", cols)

    # Use these as the palette and shared min/max
    custom_palette <- cols
    combined_min <- min(brks, na.rm = TRUE)
    combined_max <- max(brks, na.rm = TRUE)

  } else {
    # Original behavior: linear min/max over both matrices
    combined_min <- min(min(df_original, na.rm = TRUE), min(df_perturbed, na.rm = TRUE))
    combined_max <- max(max(df_original, na.rm = TRUE), max(df_perturbed, na.rm = TRUE))
  }
  # # Calculate the overall min and max values for the color scale across both matrices
  # combined_min <- min(min(df_original, na.rm = TRUE), min(df_perturbed, na.rm = TRUE))
  # combined_max <- max(max(df_original, na.rm = TRUE), max(df_perturbed, na.rm = TRUE))

  # # Generate heatmaps using the create_distance_heatmap function with a shared color scale
  # heatmap_original <- .create_distance_heatmap(df_original, title_original, combined_min, combined_max, custom_palette, show_x_axis = TRUE, show_y_axis = TRUE, show_legend = FALSE, custom_order)
  # heatmap_perturbed <- .create_distance_heatmap(df_perturbed, title_perturbed, combined_min, combined_max, custom_palette, show_x_axis = TRUE, show_y_axis = FALSE, show_legend = TRUE, custom_order)
  heatmap_original <- .create_distance_heatmap(
  df_original,
  title = title_original,
  min_val = combined_min,
  max_val = combined_max,
  custom_palette = custom_palette,
  show_x_axis = TRUE,
  show_y_axis = TRUE,
  show_legend = FALSE,
  custom_order = custom_order
  )

  heatmap_perturbed <- .create_distance_heatmap(
    df_perturbed,
    title = title_perturbed,
    min_val = combined_min,
    max_val = combined_max,
    custom_palette = custom_palette,
    show_x_axis = TRUE,
    show_y_axis = FALSE,
    show_legend = TRUE,
    custom_order = custom_order
  )

  # Combine heatmaps using patchwork
  combined_plot <- heatmap_original + heatmap_perturbed

  # Message to confirm heatmap generation
  message("Heatmaps have been successfully created for '", title_original, "' and '", title_perturbed, "'.")

  return(combined_plot)
}
