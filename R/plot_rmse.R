plot_rmse <- function(rmse_results) {
  par(mfrow = c(2, 3))
  for (cell_type in names(rmse_results)) {
    rmse_df <- rmse_results[[cell_type]]

    barplot(rmse_df$rmse,
            names.arg = paste(rmse_df$method1, rmse_df$method2, sep = " - "),
            main = paste("RMSE for cell type:", cell_type),
            xlab = "Method pairs",
            ylab = "RMSE",
            col = "skyblue",
            ylim = c(0, max(sapply(rmse_results, function(x) max(x$rmse))) * 1.1)
    )
  }
}
