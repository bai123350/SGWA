#' SelectNetworkGenes
#'
#' 选择将用于共表达网络分析的基因
#'
#' @param seurat_obj 一个 Seurat 对象
#' @param gene_select 如何选择基因？选择 "variable"（变量）、"fraction"（比例）、"all"（全部）或 "custom"（自定义）。
#' @param fraction 一个数值，决定了一个基因必须在多少最小细胞数中表达才能被包含。例如，fraction = 0.05 意味着一个基因必须在 5% 的细胞中表达（计数 > 0）才能被包含。
#' @param group.by Seurat 对象中的一个元数据列，用于对细胞进行分组（例如 seurat_clusters）。
#' @param gene_list 一个基因名称的字符向量，仅在 gene_select = "custom" 时使用。
#' @param assay seurat_obj 中的 Assay，用于计算模块特征基因。默认为 DefaultAssay(seurat_obj)。
#' @param wgcna_name WGCNA 实验的名称。
#'
#' @details
#' SelectNetworkGenes 允许我们指定将用于共表达网络分析的基因。
#' 这个函数由 SetupForWGCNA 调用。默认情况下，使用 VariableFeatures(seurat_obj) 中的可变特征。
#' 也可以通过 gene_list 参数和设置 gene_select="custom" 来使用自定义基因列表。
#'
#' 我们还可以通过设置 gene_select='fraction' 来识别在数据集中特定比例的细胞中表达量大于 0 的基因。
#' 例如，通过设置 fraction=0.1 和 group.by='seurat_clusters'，此函数将识别出
#' 在至少一个聚类中 10% 的细胞中表达的基因集合。
#'
#' @export
SelectNetworkGenes <- function(
    seurat_obj,
    gene_select = "variable",
    fraction = 0.05,
    group.by = NULL, # 应该是 Seurat 对象中的一列，例如 clusters
    gene_list = NULL,
    assay = NULL,
    wgcna_name = NULL) {
  # 如果未提供 wgcna_name，则使用活动的 wgcna
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }

  # 验证输入：
  if (!(gene_select %in% c("variable", "fraction", "all", "custom"))) {
    stop(paste0("无效的选择 gene_select: ", gene_select, ". 有效的 gene_select 是 variable, fraction, all, 或 custom。"))
  }

  # 获取 assay
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }

  # 获取表达矩阵
  if (CheckSeurat5()) {
    expr_mat <- SeuratObject::LayerData(seurat_obj, layer = "counts", assay = assay)
  } else {
    expr_mat <- Seurat::GetAssayData(seurat_obj, slot = "counts", assay = assay)
  }

  # 处理不同的选择策略
  if (gene_select == "fraction") {
    # 为了节省内存，分块对计数��阵进行二值化
    n_chunks <- ceiling(ncol(expr_mat) / 10000)

    if (n_chunks == 1) {
      chunks <- factor(rep(1), levels = 1)
    } else {
      chunks <- cut(1:nrow(expr_mat), n_chunks)
    }
    expr_mat <- do.call(rbind, lapply(levels(chunks), function(x) {
      cur <- expr_mat[chunks == x, ]
      cur[cur > 0] <- 1
      cur
    }))

    group_gene_list <- list()
    if (!is.null(group.by)) {
      # 识别表达的基因
      groups <- unique(seurat_obj@meta.data[[group.by]])
      for (cur_group in groups) {
        # 按此组对表达矩阵进行子集划分
        cur_expr <- expr_mat[, seurat_obj@meta.data[[group.by]] == cur_group]

        gene_filter <- Matrix::rowSums(cur_expr) >= round(fraction * ncol(cur_expr))
        group_gene_list[[cur_group]] <- rownames(seurat_obj)[gene_filter]
      }
      gene_list <- unique(unlist(group_gene_list))
    } else {
      # 识别在至少一部分细胞中表达的基因
      gene_filter <- Matrix::rowSums(expr_mat) >= round(fraction * ncol(seurat_obj))
      gene_list <- rownames(seurat_obj)[gene_filter]
    }
  } else if (gene_select == "variable") {
    gene_list <- VariableFeatures(seurat_obj)
  } else if (gene_select == "all") {
    gene_list <- rownames(seurat_obj)
  } else if (gene_select == "custom") {
    # 确保没有重复项
    gene_list <- unique(gene_list)

    # 所有选定的特征都应存在于 Seurat 对象中：
    if (!all(gene_list %in% rownames(seurat_obj))) {
      stop("一些选定的特征在 rownames(seurat_obj) 中未找到。")
    }

    # 检查基因列表的数据类型：
    if (!is.null(gene_list)) {
      if (class(gene_list) != "character") {
        stop(paste0("gene_list 的类型无效，必须是字符向量。"))
      }
    }
  }

  # 确保基因数量大于 0：
  if (length(gene_list) == 0) {
    stop("未找到基因")
  }

  # 如果基因数量非常少，则发出警告：
  if (length(gene_list) <= 100) {
    warning(paste0("选择的基因非常少 (", length(gene_list), "), 也许应该使用不同的方法来选择基因。"))
  }

  # 将基因添加到基因列表
  seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)

  # 返回更新后的 seurat 对象
  seurat_obj
}
