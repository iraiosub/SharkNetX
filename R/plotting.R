#' Plot cumulative mean convergence for network metrics.
#'
#' This function plots the cumulative mean convergence of four network metrics (largest connected component size, clustering coefficient, degree, and edge density) as the number of randomizations increases.
#' The plots are displayed in a 2x2 grid, showing how each metric converges over the randomization iterations. Optionally, the plots can be saved to a PDF file.
#'
#' @param random_lcc_sizes Numeric vector of LCC sizes from the randomizations.
#' @param random_clustering_coeffs Numeric vector of clustering coefficients from the randomizations.
#' @param random_degrees Numeric vector of degrees from the randomizations.
#' @param random_edge_densities Numeric vector of edge densities from the randomizations.
#' @param color Character, the color to use for the lines in the plots (default is "#4682B4").
#' @param output_file Character, optional path to save the plots as a PDF. If NULL, the plots will be displayed but not saved.
#' @return The function generates the plots and, if specified, saves them to a file.
#' @export
plot_metrics_convergence <- function(random_lcc_sizes, random_clustering_coeffs, random_degrees, random_edge_densities, color = "#4682B4", output_file = NULL) {

  # Compute cumulative means
  lcc_conv <- cumsum(random_lcc_sizes) / seq_along(random_lcc_sizes)
  clustering_coeffs_conv <- cumsum(random_clustering_coeffs) / seq_along(random_clustering_coeffs)
  degrees_conv <- cumsum(random_degrees) / seq_along(random_degrees)
  edge_densities_conv <- cumsum(random_edge_densities) / seq_along(random_edge_densities)

  # Check if an output file is specified, and if so, open PDF device
  if (!is.null(output_file)) {
    grDevices::pdf(output_file)
  }

  # Set up the layout for 2x2 plots
  graphics::par(mfrow = c(2, 2))

  # Plot the cumulative mean of LCC as randomizations increase
  graphics::plot(lcc_conv, type = "l", col = color, xlab = "Number of Randomizations", ylab = "LCC", main = "Convergence of LCC")

  # Plot the cumulative mean of clustering coefficients
  graphics::plot(clustering_coeffs_conv, type = "l", col = color, xlab = "Number of Randomizations", ylab = "Clustering coeff", main = "Convergence of clustering coeffs")

  # Plot the cumulative mean of degrees
  graphics::plot(degrees_conv, type = "l", col = color, xlab = "Number of Randomizations", ylab = "Degree", main = "Convergence of degrees")

  # Plot the cumulative mean of edge densities
  graphics::plot(edge_densities_conv, type = "l", col = color, xlab = "Number of Randomizations", ylab = "Edge density", main = "Convergence of edge densities")

  # If output file is specified, close the PDF device
  if (!is.null(output_file)) {
    grDevices::dev.off()
  }
}


#' Plot observed vs. random network metrics for a subset of nodes.
#'
#' This function compares several observed network metrics (largest connected component size, average degree, edge density, and global clustering coefficient) against random network metrics from degree-matched randomizations.
#' The function creates histograms of the random metrics with the observed values overlaid and combines them into a single grid plot for comparison. Optionally, the plot can be saved to a file.
#'
#' @param observed_lcc_size Numeric, the size of the largest connected component (LCC) in the observed subnetwork.
#' @param observed_avg_degree Numeric, the average degree of the observed subnetwork.
#' @param observed_edge_density Numeric, the edge density of the observed subnetwork.
#' @param observed_clustering_coeff Numeric, the global clustering coefficient of the observed subnetwork.
#' @param random_lcc_sizes Numeric vector, LCC sizes from random subnetworks.
#' @param random_avg_degrees Numeric vector, average degrees from random subnetworks.
#' @param random_edge_densities Numeric vector, edge densities from random subnetworks.
#' @param random_clustering_coeffs Numeric vector, global clustering coefficients from random subnetworks.
#' @param binwidth_lcc Numeric, bin width for the LCC size histogram.
#' @param binwidth_degree Numeric, bin width for the average degree histogram.
#' @param binwidth_density Numeric, bin width for the edge density histogram.
#' @param binwidth_clustering Numeric, bin width for the global clustering coefficient histogram.
#' @param color Character, color for the vertical line indicating the observed value in each plot.
#' @param random_fill Character, fill color for the histograms of random metrics.
#' @param output_file Character, optional path to save the combined plot as a file. If NULL, the plot will be displayed.
#' @return A ggplot object of the combined grid of histograms comparing observed and random network metrics.
#' @export
plot_randomization_metrics <- function(observed_lcc_size, observed_avg_degree, observed_edge_density, observed_clustering_coeff,
                                       random_lcc_sizes, random_avg_degrees, random_edge_densities, random_clustering_coeffs,
                                       binwidth_lcc = 2, binwidth_degree = 0.05, binwidth_density = 0.0001, binwidth_clustering = 0.01,
                                       color = "#4682B4", random_fill = "#bdcae2", output_file = NULL) {
  # Plot the LCC size
  lcc_size_gg <- ggplot2::ggplot(data.frame(lcc_size = random_lcc_sizes), ggplot2::aes(x = lcc_size)) +
    ggplot2::geom_histogram(fill = random_fill, color = "grey40", binwidth = binwidth_lcc, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = observed_lcc_size, col = color, lwd = 2, lty = 2) +
    ggplot2::annotate("text", x = observed_lcc_size, y = max(graphics::hist(random_lcc_sizes, plot = FALSE)$counts),
             label = paste("Observed LCC:", observed_lcc_size),
             color = color, hjust = 1.1) +
    ggplot2::labs(x = "Size of Largest Connected Component (LCC)", y = "Frequency") +
    cowplot::theme_cowplot()

  # Plot the Average Degree
  avg_degree_gg <- ggplot2::ggplot(data.frame(degree = random_avg_degrees), ggplot2::aes(x = degree)) +
    ggplot2::geom_histogram(fill = random_fill, color = "grey40", binwidth = binwidth_degree, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = observed_avg_degree, col = color, lwd = 2, lty = 2) +
    ggplot2::annotate("text", x = observed_avg_degree, y = max(graphics::hist(random_avg_degrees, plot = FALSE)$counts),
             label = paste("Observed Avg Degree:", round(observed_avg_degree, 2)),
             color = color, hjust = 1.1) +
    ggplot2::labs(x = "Average Degree", y = "Frequency") +
    cowplot::theme_cowplot()

  # Plot the Edge Density
  edge_density_gg <- ggplot2::ggplot(data.frame(edge_density = random_edge_densities), ggplot2::aes(x = edge_density)) +
    ggplot2::geom_histogram(fill = random_fill, color = "grey40", binwidth = binwidth_density, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = observed_edge_density, col = color, lwd = 2, lty = 2) +
    ggplot2::annotate("text", x = observed_edge_density, y = max(graphics::hist(random_edge_densities, plot = FALSE)$counts),
             label = paste("Observed Edge Density:", round(observed_edge_density, 5)),
             color = color, hjust = 1.1) +
    ggplot2::labs(x = "Edge Density", y = "Frequency") +
    cowplot::theme_cowplot()

  # Plot the Global Clustering Coefficient
  clustering_coeff_gg <- ggplot2::ggplot(data.frame(clustering_coeff = random_clustering_coeffs), ggplot2::aes(x = clustering_coeff)) +
    ggplot2::geom_histogram(fill = random_fill, color = "grey40", binwidth = binwidth_clustering, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = observed_clustering_coeff, col = color, lwd = 2, lty = 2) +
    ggplot2::annotate("text", x = observed_clustering_coeff, y = max(graphics::hist(random_clustering_coeffs, plot = FALSE)$counts),
             label = paste("Observed GCC:", round(observed_clustering_coeff, 3)),
             color = color, hjust = 1.1) +
    ggplot2::labs(x = "Global Clustering Coefficient", y = "Frequency") +
    cowplot::theme_cowplot() +
    ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides = "l")

  # Combine all plots into a grid for easy comparison
  combined_plot <- cowplot::plot_grid(
    lcc_size_gg, avg_degree_gg, edge_density_gg, clustering_coeff_gg,
    ncol = 2, align = 'v'
  )

  # Print or save the plot
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = combined_plot, width = 10, height = 8)
  } else {
    print(combined_plot)
  }

  return(combined_plot)
}



