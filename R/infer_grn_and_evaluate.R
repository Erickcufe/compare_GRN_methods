infer_grn <- function(seurat_obj) {
  # Get cell type list and gene expression matrix
  cell_types <- levels(seurat_obj)
  expr_matrix <- seurat_obj@assays$RNA@counts

  # Load TF list from JASPAR
  jaspar_db <- JASPAR2022::loadJasparDb()
  tf_list <- JASPAR2022::getJasparIdList(jaspar_db)

  # Initialize results list
  grn_results <- list()

  # Process each cell type
  for (cell_type in cell_types) {
    # Subset the Seurat object by cell type
    seurat_subset <- subset(seurat_obj, idents = cell_type)
    cell_type_expr_matrix <- seurat_subset@assays$RNA@counts

    # Run SCENIC
    scenic_grn <- runSCENIC_1_coexNetwork2modules(expr_matrix, tf_list)
    grn_results[[cell_type]]$SCENIC <- scenic_grn

    # Run GENIE3
    genie3_grn <- GENIE3::GENIE3(cell_type_expr_matrix, regulators = tf_list)
    grn_results[[cell_type]]$GENIE3 <- genie3_grn

    # Run GRNBoost2
    grnboost2_grn <- GRNBoost2::GRNBoost2(cell_type_expr_matrix, tf_list)
    grn_results[[cell_type]]$GRNBoost2 <- grnboost2_grn
  }

  return(grn_results)
}
