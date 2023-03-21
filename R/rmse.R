calculate_rmse <- function(graph_results) {
  # Initialize RMSE results list
  rmse_results <- list()

  for (cell_type in names(graph_results)) {
    rmse_results[[cell_type]] <- data.frame(
      method1 = character(),
      method2 = character(),
      rmse = numeric(),
      stringsAsFactors = FALSE
    )

    methods <- names(graph_results[[cell_type]])
    num_methods <- length(methods)

    for (i in 1:(num_methods - 1)) {
      for (j in (i + 1):num_methods) {
        method1 <- methods[i]
        method2 <- methods[j]

        degree1 <- graph_results[[cell_type]][[method1]]$degree
        degree2 <- graph_results[[cell_type]][[method2]]$degree

        rmse <- sqrt(mean((degree1 - degree2) ^ 2))

        rmse_results[[cell_type]] <- rbind(
          rmse_results[[cell_type]],
          data.frame(method1 = method1, method2 = method2, rmse = rmse, stringsAsFactors = FALSE)
        )
      }
    }
  }

  return(rmse_results)
}
