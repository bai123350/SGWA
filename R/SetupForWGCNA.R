#' SetupForWGCNA
#'
#' Create a slot in a Seurat object to store hdWGCNA data
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name name of the WGCNA experiment
#' @param features list of features to use for WGCNA
#' @param metacell_location name of the WGCNA experiment to copy the metacell object from
#' @param ... additional parameters to pass to SelectNetworkGenes
#'
#' @details
#' SetupForWGCNA creates a new slot in the Seurat object (seurat_obj) to store an hdWGCNA experiment
#' with a given name (wgcna_name). This function calls on SelectNetworkGenes to specify which features
#' will be used for network analysis. If there is another hdWGCNA experiment already in the seurat_obj,
#' the same metacell/metaspot object can be used by specifying the name of the hdWGCNA experiment
#' to the metacell_location parameter.
#'
#' @export
SetupForWGCNA <- function(
    seurat_obj, wgcna_name,
    features = NULL,
    metacell_location = NULL,
    ...) {
  # 将指定的 wgcna_name 设置为当前活动的 hdWGCNA 实验。
  # 这使得后续的 hdWGCNA 函数能够知道要处理哪个实验。
  seurat_obj <- SetActiveWGCNA(seurat_obj, wgcna_name)

  # 是否提供了自定义的特征（基因）列表。
  if (is.null(features)) {
    # 如果未提供特征，则使用 SelectNetworkGenes 函数自动为 WGCNA 网络选择基因。
    # 额外的参数 (...) 会传递给该函数。
    seurat_obj <- SelectNetworkGenes(seurat_obj, wgcna_name = wgcna_name, ...)
  } else {
    # 如果提供了特征列表，则使用该列表进行分析。
    # 'gene_select' 设置为 "custom" 表示使用的是预定义的列表。
    seurat_obj <- SelectNetworkGenes(
      seurat_obj,
      wgcna_name = wgcna_name,
      gene_select = "custom",
      gene_list = features
    )
  }

  # 检索为 WGCNA 分析选择的基因列表。
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name = wgcna_name)

  # 如果选择的基因数量过少（少于 500 个），则发出警告，
  # 因为这可能不足以进行稳健的网络分析。
  if (length(genes_use) < 500) {
    warning(paste0(length(genes_use), " features selected. You may wish to change the settings to include more features in your analysis."))
  }

  # 如果提供了现有元细胞对象的位置，则重用它。
  if (!is.null(metacell_location)) {
    # 验证指定的 metacell_location 是否存在于 Seurat 对象中。
    if (!(metacell_location %in% names(seurat_obj@misc))) {
      stop("metacell_location not found in seurat_obj@misc")
    }
    # 将元细胞对象从指定位置复制到当前的 WGCNA 实验中。
    seurat_obj <- SetMetacellObject(seurat_obj, metacell_location, wgcna_name)
  }

  # 返回修改后的 Seurat 对象，其中包含新的 hdWGCNA 实验设置。
  seurat_obj
}
