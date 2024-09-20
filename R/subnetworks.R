#' Perform randomization and measure network properties for a subnetwork.
#'
#' This function computes various network metrics (largest connected component size, average degree, clustering coefficient, and edge density) for a subnetwork formed by nodes of interest from the input graph.
#' It also performs randomizations by sampling subgraphs of the same size from the network and computes the same metrics for these random subgraphs.
#' The function returns both the observed metrics for the subnetwork and distributions of metrics from the randomizations.
#'
#' @param graph An igraph object representing the full network.
#' @param nodes_of_interest A vector of node names representing the subnetwork of interest.
#' @param num_randomizations Integer, the number of random subgraphs to generate (default is 1000).
#' @param seed Integer, the random seed for reproducibility (default is 42).
#' @return A list containing:
#'   - `observed_metrics`: A data frame with the observed network metrics for the subnetwork, including `lcc_size`, `avg_degree`, `clustering_coeff`, and `edge_density`.
#'   - `random_lcc_sizes`: A numeric vector of LCC sizes from the random subgraphs.
#'   - `random_degrees`: A numeric vector of average degrees from the random subgraphs.
#'   - `random_clustering_coeffs`: A numeric vector of clustering coefficients from the random subgraphs.
#'   - `random_edge_densities`: A numeric vector of edge densities from the random subgraphs.
#' @export
run_randomization <- function(graph, nodes_of_interest, num_randomizations = 1000, seed = 42) {

  # Step 1: Construct the subgraph using the nodes of interest
  observed_subgraph <- igraph::induced_subgraph(graph, vids = base::intersect(nodes_of_interest, igraph::V(graph)$name))
  num_nodes <- length(igraph::V(observed_subgraph))

  # Step 2: Initialize vectors to store the LCC sizes, degree, clustering coefficient, and edge density
  random_lcc_sizes <- numeric(num_randomizations)
  random_degrees <- numeric(num_randomizations)
  random_clustering_coeffs <- numeric(num_randomizations)
  random_edge_densities <- numeric(num_randomizations)

  # Step 3: Perform the random sampling and analysis
  base::set.seed(seed)  # Setting a seed for reproducibility
  for (i in 1:num_randomizations) {
    # Randomly sample nodes from the entire graph
    sampled_nodes <- base::sample(igraph::V(graph), num_nodes, replace = FALSE)

    # Construct the subgraph using the sampled nodes
    sampled_subgraph <- igraph::induced_subgraph(graph, sampled_nodes)

    # Calculate the size of the largest connected component in the sampled subgraph
    random_lcc_sizes[i] <- base::max(igraph::components(sampled_subgraph)$csize)

    # Calculate the average degree
    random_degrees[i] <- base::mean(igraph::degree(sampled_subgraph))

    # Calculate the global clustering coefficient (transitivity)
    random_clustering_coeffs[i] <- igraph::transitivity(sampled_subgraph, type = "global")

    # Calculate edge density
    random_edge_densities[i] <- igraph::edge_density(sampled_subgraph)
  }

  # Step 4: Calculate the metrics for the observed subgraph
  observed_lcc_size <- base::max(igraph::components(observed_subgraph)$csize)
  observed_avg_degree <- base::mean(igraph::degree(observed_subgraph))
  observed_clustering_coeff <- igraph::transitivity(observed_subgraph, type = "global")
  observed_edge_density <- igraph::edge_density(observed_subgraph)

  # Step 5: Return a list containing observed metrics and distributions from random samples
  results <- list(
    observed_metrics = data.frame(
      lcc_size = observed_lcc_size,
      avg_degree = observed_avg_degree,
      clustering_coeff = observed_clustering_coeff,
      edge_density = observed_edge_density
    ),
    random_lcc_sizes = random_lcc_sizes,
    random_degrees = random_degrees,
    random_clustering_coeffs = random_clustering_coeffs,
    random_edge_densities = random_edge_densities
  )

  return(results)
}


#' Perform degree-matched randomization and compute network metrics.
#'
#' This function performs degree-matched randomization on the input graph by sampling subgraphs that match the degree distribution of a set of nodes of interest.
#' For each randomization, network metrics (largest connected component size, average degree, edge density, and global clustering coefficient) are calculated.
#' The function returns both the metrics for the observed subgraph and the distributions of metrics from the randomized subgraphs.
#'
#' @param graph An igraph object representing the full network.
#' @param nodes_of_interest A character vector of node names representing the subnetwork of interest.
#' @param num_randomizations Integer, the number of randomizations to perform (default is 1000).
#' @param seed Integer, the random seed for reproducibility (default is 42).
#' @return A list containing:
#'   - `random_metrics`: A list of randomization metrics, including `lcc_sizes`, `avg_degrees`, `edge_densities`, and `clustering_coeffs`.
#'   - `observed_metrics`: A list of observed metrics for the original subnetwork, including `lcc_size`, `avg_degree`, `edge_density`, and `clustering_coeff`.
#' @export
run_degree_matched_randomization <- function(graph, nodes_of_interest, num_randomizations = 1000, seed = 42) {

  # Set seed for reproducibility
  base::set.seed(seed)

  # Step 1: Get degree of nodes of interest in the context of the entire network
  nodes_of_interest_degrees <- igraph::degree(graph, v = igraph::V(graph)[igraph::V(graph)$name %in% nodes_of_interest])

  # Step 2: Count how many nodes of each degree are present in the nodes of interest
  degree_counts <- base::table(nodes_of_interest_degrees)

  # Step 3: Build the pool of candidate nodes for each degree
  candidate_node_pool <- list()
  for (d in base::names(degree_counts)) {
    # Find nodes with the same degree or within a small range
    candidate_nodes <- igraph::V(graph)[igraph::degree(graph) == base::as.numeric(d)]

    # If no exact match, allow some range in degree matching
    if (length(candidate_nodes) == 0) {
      candidate_nodes <- igraph::V(graph)[igraph::degree(graph) >= (base::as.numeric(d) - 1) & igraph::degree(graph) <= (base::as.numeric(d) + 1)]
    }

    # Store the candidate nodes for this degree in a list
    candidate_node_pool[[d]] <- candidate_nodes
  }

  # Initialize vectors to store randomization metrics
  random_lcc_sizes <- base::numeric(num_randomizations)
  random_avg_degrees <- base::numeric(num_randomizations)
  random_edge_densities <- base::numeric(num_randomizations)
  random_clustering_coeffs <- base::numeric(num_randomizations)

  # To track how many smOOPs nodes are present in the sampled networks
  num_smoops_in_sampled_networks <- base::numeric(num_randomizations)

  # Store hashes of sampled subgraphs to ensure uniqueness
  sampled_subgraph_hashes <- list()

  # Step 4: Perform randomizations
  for (i in 1:num_randomizations) {
    repeat {
      sampled_nodes <- c()

      # For each degree in nodes_of_interest, sample the same number of nodes from the precomputed pool
      for (d in base::names(degree_counts)) {
        num_to_sample <- degree_counts[d]
        sampled_nodes <- c(sampled_nodes, base::sample(candidate_node_pool[[d]], num_to_sample, replace = FALSE))
      }

      # Create the subgraph with the sampled nodes
      sampled_subgraph <- igraph::induced_subgraph(graph, sampled_nodes)

      # Check how many of the sampled nodes are in the original smOOPs group
      num_smoops_in_sampled_networks[i] <- base::sum(igraph::V(sampled_subgraph)$name %in% nodes_of_interest)

      # Hash the subgraph to ensure uniqueness
      subgraph_hash <- digest::digest(sampled_subgraph)

      # Check if the subgraph has already been sampled
      if (!(subgraph_hash %in% sampled_subgraph_hashes)) {
        sampled_subgraph_hashes[[base::length(sampled_subgraph_hashes) + 1]] <- subgraph_hash
        break  # Exit the repeat loop if subgraph is unique
      }
    }

    # Calculate metrics for the random subgraph
    random_lcc_sizes[i] <- base::max(igraph::components(sampled_subgraph)$csize)
    random_avg_degrees[i] <- base::mean(igraph::degree(sampled_subgraph))
    random_edge_densities[i] <- igraph::edge_density(sampled_subgraph)
    random_clustering_coeffs[i] <- igraph::transitivity(sampled_subgraph, type = "global")
  }

  # Step 5: Calculate observed metrics for the nodes of interest
  observed_subgraph <- igraph::induced_subgraph(graph, igraph::V(graph)[igraph::V(graph)$name %in% nodes_of_interest])

  observed_lcc_size <- base::max(igraph::components(observed_subgraph)$csize)
  observed_avg_degree <- base::mean(igraph::degree(observed_subgraph))
  observed_edge_density <- igraph::edge_density(observed_subgraph)
  observed_clustering_coeff <- igraph::transitivity(observed_subgraph, type = "global")

  # Return both the randomization and observed metrics
  return(list(
    random_metrics = list(
      lcc_sizes = random_lcc_sizes,
      avg_degrees = random_avg_degrees,
      edge_densities = random_edge_densities,
      clustering_coeffs = random_clustering_coeffs,
      num_smoops_in_sampled_networks = num_smoops_in_sampled_networks
    ),
    observed_metrics = list(
      lcc_size = observed_lcc_size,
      avg_degree = observed_avg_degree,
      edge_density = observed_edge_density,
      clustering_coeff = observed_clustering_coeff
    )
  ))
}



#' Perform global degree-preserving randomization and compute network metrics.
#'
#' This function performs a degree-preserving randomization of the input graph by rewiring edges while maintaining the degree distribution.
#' It then compares various network metrics (largest connected component size, average degree, edge density, and global clustering coefficient)
#' between the observed smOOPs subgraph and the subgraphs generated from the randomized networks. The results of the randomizations and the observed metrics are returned.
#'
#' @param graph An igraph object representing the full network.
#' @param smoops_nodes A character vector of node names corresponding to the smOOPs subgraph.
#' @param num_randomizations Integer, the number of randomizations to perform (default is 1000).
#' @param niter Integer, the number of iterations for the degree-preserving rewiring algorithm (default is 1000).
#' @param seed Integer, the random seed for reproducibility (default is 42).
#' @return A list containing two elements:
#'   - `random_metrics`: A list of randomization metrics, including `lcc_sizes`, `avg_degrees`, `edge_densities`, and `clustering_coeffs`.
#'   - `observed_metrics`: A list of observed metrics for the original smOOPs subgraph, including `lcc_size`, `avg_degree`, `edge_density`, and `clustering_coeff`.
#' @export
run_degree_preserving_randomization <- function(graph, smoops_nodes, num_randomizations = 1000, niter = 1000, seed = 42) {

  # Set seed for reproducibility
  base::set.seed(seed)

  # Initialize vectors to store randomization metrics
  random_lcc_sizes <- base::numeric(num_randomizations)
  random_avg_degrees <- base::numeric(num_randomizations)
  random_edge_densities <- base::numeric(num_randomizations)
  random_clustering_coeffs <- base::numeric(num_randomizations)

  # Step 1: Get the original smOOPs subgraph
  observed_subgraph <- igraph::induced_subgraph(graph, igraph::V(graph)[igraph::V(graph)$name %in% smoops_nodes])

  # Step 2: Calculate observed metrics for the smOOPs subgraph
  observed_lcc_size <- base::max(igraph::components(observed_subgraph)$csize)
  observed_avg_degree <- base::mean(igraph::degree(observed_subgraph))
  observed_edge_density <- igraph::edge_density(observed_subgraph)
  observed_clustering_coeff <- igraph::transitivity(observed_subgraph, type = "global")

  # Step 3: Perform randomizations by rewiring the entire network
  for (i in 1:num_randomizations) {
    # Rewire the entire network while preserving degree
    rewired_graph <- igraph::rewire(graph, with = igraph::keeping_degseq(niter = niter))

    # Create a subgraph using the original smOOPs nodes
    rewired_subgraph <- igraph::induced_subgraph(rewired_graph, igraph::V(rewired_graph)[igraph::V(rewired_graph)$name %in% smoops_nodes])

    # Calculate metrics for the rewired subgraph
    random_lcc_sizes[i] <- base::max(igraph::components(rewired_subgraph)$csize)  # Largest connected component size
    random_avg_degrees[i] <- base::mean(igraph::degree(rewired_subgraph))  # Average degree
    random_edge_densities[i] <- igraph::edge_density(rewired_subgraph)  # Edge density
    random_clustering_coeffs[i] <- igraph::transitivity(rewired_subgraph, type = "global")  # Global clustering coefficient
  }

  # Step 4: Return both the randomization and observed metrics
  return(list(
    random_metrics = list(
      lcc_sizes = random_lcc_sizes,
      avg_degrees = random_avg_degrees,
      edge_densities = random_edge_densities,
      clustering_coeffs = random_clustering_coeffs
    ),
    observed_metrics = list(
      lcc_size = observed_lcc_size,
      avg_degree = observed_avg_degree,
      edge_density = observed_edge_density,
      clustering_coeff = observed_clustering_coeff
    )
  ))
}
