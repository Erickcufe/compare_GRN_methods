process_grn_results <- function(grn_results) {
  # Initialize results list
  graph_results <- list()

  for (cell_type in names(grn_results)) {
    graph_results[[cell_type]] <- list()
    for (method in names(grn_results[[cell_type]])) {
      # Convert GRN to igraph object
      grn <- grn_results[[cell_type]][[method]]
      grn_graph <- graph_from_data_frame(grn, directed = TRUE)

      # Calculate graph properties
      graph_results[[cell_type]][[method]]$degree <- degree(grn_graph)
      graph_results[[cell_type]][[method]]$degree_in <- degree(grn_graph, mode = "in")
      graph_results[[cell_type]][[method]]$degree_out <- degree(grn_graph, mode = "out")
      graph_results[[cell_type]][[method]]$edge_betweenness <- edge_betweenness(grn_graph)
    }
  }

  return(graph_results)
}
