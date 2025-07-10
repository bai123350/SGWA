#' Run sgWGCNA Pipeline
#'
#' A function to encapsulate the hdWGCNA analysis pipeline for a single group.
#'
#' @param seurat_obj The input Seurat object.
#' @param gene_select The method for gene selection. Can be "fraction" or "variable".
#' @param fraction The fraction of cells a gene must be expressed in to be included.
#' @param wgcna_name The name for the WGCNA analysis stored in the Seurat object.
#' @param metacell_group.by The metadata columns to group cells for metacell analysis.
#' @param metacell_ident.group1 The primary metadata column for cell identity.
#' @param metacell_ident.group2 The secondary metadata column for grouping (e.g., sample).
#' @param k The number of nearest neighbors for metacell construction.
#' @param max_shared The maximum shared neighbors for metacell construction.
#' @param min_cells The minimum number of cells required for a group to be included.
#' @param target_metacells The target number of metacells to create.
#' @param reduction The dimensionality reduction to use for metacell construction.
#' @param dat_expr_group_name The specific group name from `metacell_ident.group1` to analyze.
#' @param network_type The type of network to construct (e.g., "signed").
#' @param n_hubs The number of hub genes to identify per module.
#' @param n_genes_score The number of genes to use for module scoring.
#' @param assay The assay to use for expression data.
#' @param plot_path The directory path to save the output plots.
#' @param plot_height The height of the saved plots in inches.
#' @param plot_width The width of the saved plots in inches.
#' @param cols The number of columns for arranging feature plots.
#' @return The processed Seurat object with WGCNA results.
#' @keywords WGCNA scRNA-seq
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- run_sgwgcna_pipeline(
#'   seurat_obj,
#'   dat_expr_group_name = "CD4 T-cells",
#'   min_cells = 15,
#'   plot_path = "wgcna_plots_cd4/"
#' )
#' }
run_sgwgcna_pipeline <- function(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "ASC_consensus",
    metacell_group.by = c("SCT_snn_res.0.5", "sample"),
    metacell_ident.group1 = "SCT_snn_res.0.5",
    metacell_ident.group2 = "sample",
    k = 25,
    max_shared = 12,
    min_cells = 50,
    target_metacells = 250,
    reduction = "harmony",
    dat_expr_group_name = "1",
    network_type = "signed",
    n_hubs = 10,
    n_genes_score = 25,
    assay = NULL,
    meta_yazu = "sample",
    plot_path = "./plots/",
    plot_height = 10,
    plot_width = 10,
    cols = 6) {
  # Create plot directory if it doesn't exist
  if (!dir.exists(plot_path)) {
    dir.create(plot_path, recursive = TRUE)
  }

  # Setup for WGCNA
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = gene_select,
    fraction = fraction,
    wgcna_name = wgcna_name
  )

  # Metacells by groups
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = metacell_group.by,
    ident.group = metacell_ident.group1,
    k = k,
    max_shared = max_shared,
    min_cells = min_cells,
    target_metacells = target_metacells,
    reduction = reduction
  )

  seurat_obj <- NormalizeMetacells(seurat_obj)

  # 获取 assay
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }

  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = dat_expr_group_name,
    group.by = metacell_ident.group1,
    assay = assay,
    layer = "data"
  )

  # Test different soft powers
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = network_type
  )

  # plot the results:
  plot_list <- PlotSoftPowers(seurat_obj)
  pdf(file.path(plot_path, "soft_powers.pdf"), height = plot_height, width = plot_width)
  print(wrap_plots(plot_list, ncol = 2))
  dev.off()

  power_table <- GetPowerTable(seurat_obj)
  print(head(power_table))

  # construct co-expression network:
  overwrite_tom <- file.exists(paste0("TOM/", dat_expr_group_name, "_TOM.rda"))
  seurat_obj <- ConstructNetwork(
    seurat_obj,
    tom_name = dat_expr_group_name,
    overwrite_tom = overwrite_tom
  )

  pdf(file.path(plot_path, "dendrogram.pdf"), height = plot_height, width = plot_width)
  print(PlotDendrogram(seurat_obj, main = paste0(dat_expr_group_name, " hdWGCNA Dendrogram")))
  dev.off()

  TOM <- GetTOM(seurat_obj)

  # need to run ScaleData first or else harmony throws an error:
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

  # compute all MEs in the full single-cell dataset
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars = metacell_ident.group2
  )

  # harmonized module eigengenes:
  hMEs <- GetMEs(seurat_obj)

  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized = FALSE)

  # compute eigengene-based connectivity (kME):
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = metacell_ident.group1, group_name = dat_expr_group_name
  )

  # rename the modules
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste0(dat_expr_group_name, "-M")
  )

  # plot genes ranked by kME for each module
  p <- PlotKMEs(seurat_obj, ncol = cols)
  pdf(file.path(plot_path, "kMEs.pdf"), height = plot_height, width = plot_width)
  print(p)
  dev.off()

  # get the module assignment table:
  modules <- GetModules(seurat_obj) %>% subset(module != "grey")
  print(head(modules[, 1:6]))

  # get hub genes
  hub_df <- GetHubGenes(seurat_obj, n_hubs = n_hubs)
  print(head(hub_df))

  # compute gene scoring for the top hub genes by kME for each module
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = n_genes_score,
    method = "UCell"
  )

  # make a featureplot of hMEs for each module
  plot_list_hMEs <- ModuleFeaturePlot(
    seurat_obj,
    features = "hMEs", # plot the hMEs
    order = TRUE # order so the points with highest hMEs are on top
  )
  pdf(file.path(plot_path, "hMEs_featureplot.pdf"), height = plot_height, width = plot_width * 2)
  print(wrap_plots(plot_list_hMEs, ncol = cols))
  dev.off()

  # make a featureplot of hub scores for each module
  plot_list_scores <- ModuleFeaturePlot(
    seurat_obj,
    features = "scores", # plot the hub gene scores
    order = "shuffle", # order so cells are shuffled
    ucell = TRUE # depending on Seurat vs UCell for gene scoring
  )
  pdf(file.path(plot_path, "scores_featureplot.pdf"), height = plot_height, width = plot_width * 2)
  print(wrap_plots(plot_list_scores, ncol = cols))
  dev.off()

  seurat_obj$cluster <- seurat_obj@meta.data[[meta_yazu]]

  barcodes_subset <- rownames(
    seurat_obj@meta.data[
      seurat_obj@meta.data[[metacell_ident.group1]] == dat_expr_group_name,
    ]
  )

  pdf(file.path(plot_path, "radar_plot.pdf"), height = plot_height, width = plot_width)
  print(
    ModuleRadarPlot(
      seurat_obj,
      group.by = "cluster",
      barcodes = barcodes_subset,
      axis.label.size = 4,
      grid.label.size = 4
    )
  )
  dev.off()

  pdf(file.path(plot_path, "ModuleCorre.pdf"), height = plot_height, width = plot_width)
  print(ModuleCorrelogram(seurat_obj))
  dev.off()

  # get hMEs from seurat object
  MEs <- GetMEs(seurat_obj, harmonized = TRUE)
  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]

  # add hMEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

  # plot with Seurat's DotPlot function
  p <- DotPlot(seurat_obj, features = mods, group.by = metacell_ident.group1)

  p <- p +
    RotatedAxis() +
    scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
  ggsave(
    filename = file.path(plot_path, "dotplot_modules.pdf"),
    plot = p,
    height = plot_height,
    width = plot_width
  )

  return(seurat_obj)
}
