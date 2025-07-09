############################
# 活动 WGCNA
###########################

#' SetActiveWGCNA
#'
#' @param seurat_obj 一个 Seurat 对象
#' @param wgcna_name seurat_obj@misc 插槽中 hdWGCNA 实验的名称
#' @keywords scRNA-seq
#' @export
SetActiveWGCNA <- function(seurat_obj, wgcna_name) {
  # 设置 active_wgcna 变量
  seurat_obj@misc$active_wgcna <- wgcna_name

  # 如果该 WGCNA 的空列表尚不存在，则初始化它
  if (!(wgcna_name %in% names(seurat_obj@misc))) {
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]] <- list()
  }
  seurat_obj
}

#' GetActiveWGCNA
#'
#' @param seurat_obj 一个 Seurat 对象
#' @keywords scRNA-seq
#' @export
GetActiveWGCNA <- function(seurat_obj) {
  seurat_obj@misc[[seurat_obj@misc$active_wgcna]]
}

#' GetActiveWGCNAName
#'
#' @param seurat_obj 一个 Seurat 对象
#' @keywords scRNA-seq
#' @export
GetActiveWGCNAName <- function(seurat_obj) {
  seurat_obj@misc$active_wgcna
}

#' CheckWGCNAName
#'
#' @param seurat_obj 一个 Seurat 对象
#' @keywords scRNA-seq
#' @export
CheckWGCNAName <- function(seurat_obj, wgcna_name) {
  check <- wgcna_name %in% names(seurat_obj@misc)
  if (!check) {
    stop(paste0("提供的 wgcna_name 无效: ", wgcna_name))
  }
}


# # 获取任何 WGCNA 数据，但默认获取活动数据
# GetWGCNA <- function(seurat_obj, wgcna_name=NULL){

#   # 测试 wgcna_name 是否有效 (TODO)

#   # 如果未提供 wgcna_name，则从活动分析中获取数据
#   if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

#   seurat_obj@misc[[wgcna_name]]
# }

############################
# 元细胞对象
###########################

#' SetMetacellObject
#'
#' @param seurat_obj 一个 Seurat 对象
#' @param metacell_obj 元细胞 Seurat 对象
#' @param wgcna_name seurat_obj@misc 插槽中 hdWGCNA 实验的名称
#' @keywords scRNA-seq
#' @export
SetMetacellObject <- function(seurat_obj, metacell_obj, wgcna_name = NULL) {
  # 如果未提供 wgcna_name，则从活动分析中获取数据
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # 将元细胞对象添加到 Seurat 对象
  seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj <- metacell_obj
  seurat_obj
}

#' GetMetacellObject
#'
#' @param seurat_obj 一个 Seurat 对象
#' @param wgcna_name seurat_obj@misc 插槽中 hdWGCNA 实验的名称
#' @keywords scRNA-seq
#' @export
GetMetacellObject <- function(seurat_obj, wgcna_name = NULL) {
  # 如果未提供 wgcna_name，则从活动分析中获取数据
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  input_class <- class(seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj)
  if (input_class == "Seurat") {
    return(seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj)
  } else if (input_class == "character") {
    metacell_location <- seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj
    return(seurat_obj@misc[[metacell_location]]$wgcna_metacell_obj)
  } else {
    return(NULL)
  }
}

############################
# WGCNA 基因
###########################

#' SetWGCNAGenes
#'
#' @param seurat_obj 一个 Seurat 对象
#' @param gene_list 用于 WGCNA 的基因向量
#' @param wgcna_name seurat_obj@misc 插槽中 hdWGCNA 实验的名称
#' @keywords scRNA-seq
#' @export
SetWGCNAGenes <- function(seurat_obj, gene_list, wgcna_name = NULL) {
  # 如果未提供 wgcna_name，则从活动分析中获取数据
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # 将基因列表添加到 Seurat 对象
  seurat_obj@misc[[wgcna_name]]$wgcna_genes <- gene_list
  seurat_obj
}

#' GetWGCNAGenes
#'
#' @param seurat_obj 一个 Seurat 
#' @param wgcna_name seurat_obj@misc 插槽中 hdWGCNA 实验的名称
#' @keywords scRNA-seq
#' @export
GetWGCNAGenes <- function(seurat_obj, wgcna_name = NULL) {
  # 如果未提供 wgcna_name，则从活动分析中获取数据
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_genes
}

############################
# datExpr
###########################


#' SetDatExpr
#'
#' 此函数指定用于共表达网络分析的基因表达矩阵。
#'
#' @param seurat_obj 一个 Seurat 对象
#' @param group_name 一个字符串，包含提供的 group.by 列中或 Seurat Idents 中的一个组。可以提供一个字符向量以一次选择多个组。
#' @param use_metacells 一个逻辑值，确定我们是使用元细胞 (TRUE) 还是完整的表达矩阵 (FALSE)
#' @param group.by 一个字符串，包含 Seurat 对象中具有细胞组（簇、细胞类型等）的列的名称。如果为 NULL（默认值），hdWGCNA 使用 Seurat Idents 作为组。
#' @param multi.group.by 一个字符串，包含 Seurat 对象中具有用于共识 WGCNA 的组（数据集、样本、条件等）的列的名称
#' @param multi_group_name 一个字符串，包含 multi.group.by 列中存在的组的名称。
#' @param assay Seurat 对象中分析的名称
#' @param slot 用于提取聚合数据的插槽。默认 = 'counts'。插槽与 Seurat v4 一起使用，而不是层。
#' @param layer 用于提取聚合数据的层。默认 = 'counts'。层与 Seurat v5 一起使用，而不是插槽。
#' @param mat 一个包含基因表达数据的矩阵。使用此参数提供矩阵会忽略其他选项。这几乎专门用于伪批量分析。
#' @param features 要用于覆盖先前设置的特征的特征列表。
#' @param wgcna_name 一个字符串，包含 seurat_obj@misc 中 WGCNA 插槽的名称。默认 = NULL，它检索当前活动的 WGCNA 数据
#' @details
#' SetDatExpr 是 hdWGCNA 管道的一个关键函数，它确定将用于网络分析的基因表达矩阵。
#' 我们通常使用此函数使用 group.by 参数选择一个细胞类型或一组细胞类型进行网络分析，但我们提供
#' 用于进一步定制的其他参数。
#' @keywords scRNA-seq
#' @export
SetDatExpr <- function(
    seurat_obj,
    group_name,
    use_metacells = TRUE,
    group.by = NULL,
    multi.group.by = NULL,
    multi_group_name = NULL,
    return_seurat = TRUE,
    assay = NULL,
    slot = "data",
    layer = "data",
    mat = NULL,
    features = NULL,
    wgcna_name = NULL,
    ...) {
  # 如果未提供 wgcna_name，则从活动分析中取数据
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
    warning(paste0("未指定分析，尝试使用分析 ", assay))
  }

  # 检查所选分析是否在 seurat 对象中
  if (!(assay %in% names(seurat_obj))) {
    stop(paste0("无效的分析选择: ", assay, " 在 Assays(seurat_obj) 中未找到。"))
  }

  # 检查插槽是否有效
  if (!(slot %in% c("counts", "data", "scale.data"))) {
    stop("无效的插槽选择。有效选择是 counts、data 或 scale.data。")
  }

  # 从 seurat 对象获取参数
  params <- GetWGCNAParams(seurat_obj, wgcna_name)

  if (is.null(features)) {
    genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  } else {
    if (all(features %in% rownames(seurat_obj))) {
      genes_use <- features
    } else {
      stop("在 rownames(seurat_obj) 中未找到某些特征。")
    }
  }

  # 是否提供了矩阵？
  if (is.null(mat)) {
    # 获取元细胞对象
    m_obj <- GetMetacellObject(seurat_obj, wgcna_name)

    # 使用元细胞还是整个 seurat 对象？
    if (use_metacells & !is.null(m_obj)) {
      s_obj <- m_obj
    } else {
      if (is.null(m_obj)) {
        warning("未找到元细胞 Seurat 对象。改用完整的 Seurat 对象。")
      }
      s_obj <- seurat_obj
    }

    # 从 seurat 对象获取元数据：
    seurat_meta <- s_obj@meta.data

    # 检查 group.by 参数
    if (!is.null(group.by)) {
      # 检查 group.by 是否在 Seurat 对象和元细胞对象中：
      if (!(group.by %in% colnames(s_obj@meta.data))) {
        m_cell_message <- ""
        if (use_metacells) {
          m_cell_message <- "metacell"
        }
        stop(paste0(group.by, " 在 ", m_cell_message, " Seurat 对象的元数据中未找到"))
      }

      # 检查所选组是否在 Seurat 对象中：
      if (!all(group_name %in% s_obj@meta.data[[group.by]])) {
        groups_not_found <- group_name[!(group_name %in% s_obj@meta.data[[group.by]])]
        stop(
          paste0("group_name 中的某些组在 seurat_obj 中未找到: ", paste(groups_not_found, collapse = ", "))
        )
      }
    }

    # 检查 multi.group.by 参数
    if (!is.null(multi.group.by)) {
      # 检查 group.by 是否在 Seurat 对象和元细胞对象中：
      if (!(multi.group.by %in% colnames(s_obj@meta.data))) {
        m_cell_message <- ""
        if (use_metacells) {
          m_cell_message <- "metacell"
        }
        stop(paste0(multi.group.by, " 在 ", m_cell_message, " Seurat 对象的元数据中未找到"))
      }

      # 检查所选组是否在 Seurat 对象中：
      if (!all(multi_group_name %in% s_obj@meta.data[[multi.group.by]])) {
        groups_not_found <- multi_group_name[!(multi_group_name %in% s_obj@meta.data[[multi.group.by]])]
        stop(
          paste0("group_name 中的某些组在 seurat_obj 中未找到: ", paste(groups_not_found, collapse = ", "))
        )
      }
    }

    # 用于按簇/细胞类型分组的列 # nolint
    if (!is.null(group.by)) {
      seurat_meta <- seurat_meta %>% subset(get(group.by) %in% group_name)
    }

    # 如果是 multiExpr，则进一步子集化：
    if (!is.null(multi.group.by)) {
      seurat_meta <- seurat_meta %>% subset(get(multi.group.by) %in% multi_group_name)
    }

    # 获取要使用的细胞列表
    cells <- rownames(seurat_meta)

    # 从 seurat 对象获取表达数据
    if (CheckSeurat5()) {
      exp <- SeuratObject::LayerData(s_obj, assay = assay, layer = layer)
    } else {
      exp <- Seurat::GetAssayData(s_obj, assay = assay, slot = slot)
    }
    datExpr <- as.data.frame(exp)[genes_use, cells]

    # 转置数据
    datExpr <- as.data.frame(t(datExpr))
  } else {
    datExpr <- mat

    # 将其转换为数据框
    if (any(class(datExpr) != "data.frame")) {
      datExpr <- as.data.frame(datExpr)
    }

    # colnames 是基因吗？
    if (!all(colnames(datExpr) %in% rownames(seurat_obj))) {
      stop("提供的矩阵的 colnames 无效。请确保 colnames 是特征（基因），并且所有这些特征都在 seurat_obj 中")
    }

    # 按 WGCNA 基因对 datExpr 进行子集化：
    genes_use <- colnames(datExpr)
    # datExpr <- datExpr[,genes_use]
  }

  if (return_seurat) {
    gene_list <- genes_use[WGCNA::goodGenes(datExpr, ...)]
    datExpr <- datExpr[, gene_list]

    # 更新 WGCNA 基因列表：
    seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)

    # 在 Seurat 对象中设置 datExpr
    seurat_obj@misc[[wgcna_name]]$datExpr <- datExpr
    out <- seurat_obj
  } else {
    out <- datExpr
  }
  out
}



#' GetDatExpr
#'
#' 此函数获取 WGCNA 表达矩阵。
#'
#' @param seurat_obj 一个 Seurat 对象
#' @keywords scRNA-seq
#' @export
GetDatExpr <- function(seurat_obj, wgcna_name = NULL) {
  # 如果未提供 wgcna_name，则从活动分析中获取数据
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$datExpr
}


#' SetMultiExpr
#'
#' 此函数根据元细胞表达矩阵、完整表达矩阵或提供的伪批量表达矩阵，为共识 WGCNA 设置表达矩阵输入。
#'
#' @param seurat_obj 一个 Seurat 对象
#' @param group_name 一个字符串，包含提供的 group.by 列中或 Seurat Idents 中的一个组。
#' @param use_metacells 一个逻辑值，确定我们是使用元细胞 (TRUE) 还是完整的表达矩阵 (FALSE)
#' @param group.by 一个字符串，包含 Seurat 对象中具有细胞组（簇、细胞类型等）的列的名称。如果为 NULL（默认值），hdWGCNA 使用 Seurat Idents 作为组。
#' @param multi.group.by 一个字符串，包含 Seurat 对象中具有用于共识 WGCNA 的组（数据集、样本、条件等）的列的名称
#' @param multi_groups 一个字符向量，包含要选择的组的名称
#' @param assay Seurat 对象中分析的名称
#' @param slot Seurat 对象中插槽的名称（counts、data）
#' @param layer 用于提取聚合数据的层。默认 = 'counts'。层与 Seurat v5 一起使用，而不是插槽。
#' @param mat 一个包含基因表达数据的矩阵。使用此参数提供矩阵会忽略其他选项。这几乎专门用于伪批量分析。
#' @param wgcna_name 一个字符串，包含 seurat_obj@misc 中 WGCNA 插槽的名称。默认 = NULL，它检索当前活动的 WGCNA 数据
#' @keywords scRNA-seq
#' @export
SetMultiExpr <- function(
    seurat_obj,
    group_name,
    use_metacells = TRUE,
    group.by = NULL,
    multi.group.by = NULL,
    multi_groups = NULL,
    assay = NULL,
    slot = "data",
    layer = "data",
    mat = NULL,
    mat_group_delim = 3,
    wgcna_name = NULL,
    ...) {
  # 如果未提供 wgcna_name，则从活动分析中获取数据
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # 获取 WGCNA 基因：
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  gene_names <- GetWGCNAGenes(seurat_obj, wgcna_name)

  s_obj <- seurat_obj

  # 获取分析
  if (is.null(assay)) {
    assay <- DefaultAssay(s_obj)
    warning(paste0("未指定分析，尝试使用分析 ", assay))
  }

  # 如果用户未指定，则获取存在的不同组：
  if (is.null(multi_groups)) {
    multi_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
  } else {
    seurat_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
    if (sum(multi_groups %in% seurat_groups) != length(multi_groups)) {
      stop("在 seurat_obj@meta.data[,multi.group.by] 中未找到 multi_groups 中指定的某些或所有组")
    }
  }

  # 是否提供了矩阵？
  if (is.null(mat)) {
    # 使用元细胞还是整个 seurat 对象？
    if (use_metacells) {
      s_obj <- GetMetacellObject(seurat_obj, wgcna_name)
    } else {
      s_obj <- seurat_obj
    }

    # 获取每个组的 datExpr
    datExpr_list <- lapply(multi_groups, function(cur_group) {
      cur_datExpr <- SetDatExpr(
        seurat_obj,
        group_name = group_name,
        group.by = group.by,
        multi.group.by = multi.group.by,
        multi_group_name = cur_group,
        return_seurat = FALSE,
        use_metacells = use_metacells,
        wgcna_name = wgcna_name,
        assay = assay,
        slot = slot,
        layer = layer
      )
      as.matrix(cur_datExpr)
    })
  } else {
    sample_groups <- do.call(rbind, strsplit(rownames(mat), ":"))[, mat_group_delim]
    datExpr_list <- list()
    for (cur_group in multi_groups) {
      cur_datExpr <- as.data.frame(mat[which(sample_groups == cur_group), ])
      datExpr_list[[cur_group]] <- cur_datExpr
    }
  }

  # 转换为 multiExpr，获取好基因：
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  genes_use <- WGCNA::goodGenesMS(multiExpr)
  gene_names <- gene_names[genes_use]

  # 按好基因对 multiExpr 进行子集化：
  datExpr_list <- lapply(1:length(multiExpr), function(i) {
    multiExpr[[i]]$data[, genes_use]
  })
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  names(multiExpr) <- multi_groups

  # 更新 WGCNA 基因列表：
  seurat_obj <- SetWGCNAGenes(seurat_obj, gene_names, wgcna_name)

  # 在 Seurat 对象中设置 multiExpr
  seurat_obj@misc[[wgcna_name]]$multiExpr <- multiExpr
  seurat_obj
}


#' GetMultiExpr
#'
#' 此函数从元细胞对象获取表达矩阵。
#'
#' @param seurat_obj 一个 Seurat 对象
#' @keywords scRNA-seq
#' @export
GetMultiExpr <- function(seurat_obj, wgcna_name = NULL) {
  # 如果未提供 wgcna_name，则从活动分析中获取数据
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$multiExpr
}

############################
# WGCNA params
###########################


#' SetMetacellParams
#'
#' @param seurat_obj A Seurat object
#' @param params list of WGCNA parameters
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMetacellParams <- function(seurat_obj, params, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$metacell_params <- params
  seurat_obj
}

#' GetMetacellParams
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMetacellParams <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$metacell_params
}

#' SetWGCNAParams
#'
#' @param seurat_obj A Seurat object
#' @param params list of WGCNA parameters
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetWGCNAParams <- function(seurat_obj, params, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_params <- params
  seurat_obj
}

#' GetWGCNAParams
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetWGCNAParams <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_params
}

############################
# SoftPower Table
###########################

#' SetPowerTable
#'
#' @param seurat_obj A Seurat object
#' @param power_table a dataframe containing the results of the soft power test
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetPowerTable <- function(seurat_obj, power_table, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # add power table to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable <- power_table
  seurat_obj
}

#' GetPowerTable
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetPowerTable <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable
}

############################
# WGCNA Network
###########################

#' SetNetworkData
#'
#' @param seurat_obj A Seurat object
#' @param net list of network data from WGCNA
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetNetworkData <- function(seurat_obj, net, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # add network data to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_net <- net
  seurat_obj
}


#' GetNetworkData
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetNetworkData <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$wgcna_net
}

############################
# WGCNA modules dataframe
###########################


#' SetModules
#'
#' @param seurat_obj A Seurat object
#' @param modules dataframe containing gene module assignments
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModules <- function(seurat_obj, modules, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set module df
  seurat_obj@misc[[wgcna_name]]$wgcna_modules <- modules
  seurat_obj
}

#' GetModules
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
GetModules <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_modules
}





#' SetDegrees
#'
#' @param seurat_obj A Seurat object
#' @param degree_df dataframe containing gene module assignments
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
SetDegrees <- function(seurat_obj, degree_df, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set module df
  seurat_obj@misc[[wgcna_name]]$wgcna_degrees <- degree_df
  seurat_obj
}

#' GetDegrees
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
GetDegrees <- function(seurat_obj, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$wgcna_degrees
}


#' GetHubGenes
#'
#' Extract the top N hub genes for a given set of modules. This function outputs
#' a table with the gene name, the module, and the kME for that module for the
#' top N hub genes.
#'
#' @param seurat_obj A Seurat object
#' @param n_hubs the number of hub genes to select for each module
#' @param mods list of modules, selects all modules by default
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetHubGenes <- function(
    seurat_obj,
    n_hubs = 10,
    mods = NULL,
    wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the modules table
  modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != "grey")

  if (is.null(mods)) {
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
  } else {
    if (!all(mods %in% modules$module)) {
      stop("Invalid selection for mods.")
    }
  }

  # get hub genes:
  hub_df <- do.call(rbind, lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur <- cur[, c("gene_name", "module", paste0("kME_", cur_mod))]
    names(cur)[3] <- "kME"
    cur <- dplyr::arrange(cur, desc(kME))
    cur %>% dplyr::slice_max(n = n_hubs, order_by = kME)
  }))
  rownames(hub_df) <- 1:nrow(hub_df)
  hub_df
}


############################
# Module Eigengenes
###########################

#' SetMEs
#'
#' @param seurat_obj A Seurat object
#' @param MEs dataframe or matrix containing module eigengenes
#' @param harmonized logical indicating whether MEs have been harmonized
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMEs <- function(seurat_obj, MEs, harmonized = TRUE, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # harmonized MEs?
  if (harmonized) {
    seurat_obj@misc[[wgcna_name]]$hMEs <- MEs
  } else {
    seurat_obj@misc[[wgcna_name]]$MEs <- MEs
  }
  seurat_obj
}

#' GetMEs
#'
#' Function to retrieve module eigengens from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param harmonized logical indicating whether MEs have been harmonized
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMEs <- function(seurat_obj, harmonized = TRUE, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get harmonized MEs?
  if (harmonized == TRUE && !is.null(seurat_obj@misc[[wgcna_name]]$hMEs)) {
    MEs <- seurat_obj@misc[[wgcna_name]]$hMEs
  } else {
    MEs <- seurat_obj@misc[[wgcna_name]]$MEs
  }
  MEs
}


#' SetMELoadings
#'
#' @param seurat_obj A Seurat object
#' @param loadings named numeric vector with eigengene loadings
#' @param harmonized logical indicating whether MEs have been harmonized
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMELoadings <- function(seurat_obj, loadings, harmonized = TRUE, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # harmonized MEs?
  if (harmonized) {
    seurat_obj@misc[[wgcna_name]]$hME_loadings <- c(seurat_obj@misc[[wgcna_name]]$hME_loadings, loadings)
  } else {
    seurat_obj@misc[[wgcna_name]]$ME_loadings <- c(seurat_obj@misc[[wgcna_name]]$ME_loadings, loadings)
  }
  seurat_obj
}

#' GetMELoadings
#'
#' Function to retrieve module eigengens from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param harmonized logical indicating whether MEs have been harmonized
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMELoadings <- function(seurat_obj, harmonized = TRUE, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get harmonized MEs?
  if (harmonized == TRUE && !is.null(seurat_obj@misc[[wgcna_name]]$hME_loadings)) {
    MEs <- seurat_obj@misc[[wgcna_name]]$hME_loadings
  } else {
    MEs <- seurat_obj@misc[[wgcna_name]]$ME_loadings
  }
  MEs
}


############################
# GO term table
###########################

#' SetEnrichrTable
#'
#' @param seurat_obj A Seurat object
#' @param enrich_table dataframe storing the results of running enrichr
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetEnrichrTable <- function(seurat_obj, enrich_table, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set enrichr table
  seurat_obj@misc[[wgcna_name]]$enrichr_table <- enrich_table
  seurat_obj
}



#' GetEnrichrTable
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetEnrichrTable <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$enrichr_table
}


#' SetEnrichRegulonTable
#'
#' @param seurat_obj A Seurat object
#' @param enrich_table dataframe storing the results of running enrichr
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetEnrichrRegulonTable <- function(seurat_obj, enrich_table, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # set enrichr table
  seurat_obj@misc[[wgcna_name]]$enrichr_regulon_table <- enrich_table
  seurat_obj
}



#' GetEnrichrRegulonTable
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetEnrichrRegulonTable <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$enrichr_regulon_table
}



############################
# Module Scores
###########################


#' SetModuleScores
#'
#' @param seurat_obj A Seurat object
#' @param mod_scores dataframe storing the module expression scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
SetModuleScores <- function(seurat_obj, mod_scores, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_scores <- mod_scores
  seurat_obj
}

#' GetModuleScores
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
GetModuleScores <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_scores
}

############################
# Average Module Expression
###########################

#' SetAvgModuleExpr
#'
#' @param seurat_obj A Seurat object
#' @param avg_mods dataframe storing the average expression of all genes in the same module
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @export
SetAvgModuleExpr <- function(seurat_obj, avg_mods, wgcna_name = NULL) {
  # get data from active assay if wgcna_name is not given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$avg_modules <- avg_mods
  seurat_obj
}

#' GetAvgModuleExpr
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetAvgModuleExpr <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$avg_modules
}

############################
# TOM
###########################

#' GetTOM
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTOM <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get modules
  modules <- GetModules(seurat_obj, wgcna_name)
  gene_names <- modules$gene_name

  # load TOM
  tom_files <- GetNetworkData(seurat_obj, wgcna_name)$TOMFiles

  if (!file.exists(tom_files[[1]])) {
    stop(paste0("TOM file ", tom_files[[1]], " not found. Please update path to TOM file."))
  }

  load(tom_files[[1]])

  TOM <- as.matrix(consTomDS)
  rownames(TOM) <- gene_names
  colnames(TOM) <- gene_names
  TOM
}


############################
# TF Match Matrix (not stored within the WGCNA slot)
###########################

#' SetMotifMatrix
#'
#' @param seurat_obj A Seurat object
#' @param tf_match matrix containing tf-promoter matches
#' @keywords scRNA-seq
#' @export
#' @examples SetMotifMatrix
SetMotifMatrix <- function(seurat_obj, tf_match) {
  # make a spot for the motif info if it's not already there:
  if (is.null(seurat_obj@misc$motifs)) {
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$tf_match_matrix <- tf_match
  seurat_obj
}

#' GetMotifMatrix
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples GetMotifMatrix
GetMotifMatrix <- function(seurat_obj) {
  seurat_obj@misc$motifs$tf_match_matrix
}


############################
# Motif table
###########################

#' SetMotifs
#'
#' @param seurat_obj A Seurat object
#' @param motif_df dataframe containing info about the motifs being analyzed
#' @keywords scRNA-seq
#' @export
#' @examples SetMotifs
SetMotifs <- function(seurat_obj, motif_df) {
  # make a spot for the motif info if it's not already there:
  if (is.null(seurat_obj@misc$motifs)) {
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$motif_df <- motif_df
  seurat_obj
}



#' GetMotifs
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples GetMotifs
GetMotifs <- function(seurat_obj) {
  seurat_obj@misc$motifs$motif_df
}

############################
# PFM List
###########################


#' SetPFMList
#'
#' @param seurat_obj A Seurat object
#' @param pfm_list list of pfm objects
#' @keywords scRNA-seq
#' @export
#' @examples SetPFMList
SetPFMList <- function(seurat_obj, pfm_list) {
  # make a spot for the motif info if it's not already there:
  if (is.null(seurat_obj@misc$motifs)) {
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$pfm_list <- pfm_list
  seurat_obj
}

#' GetPFMList
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples GetPFMList
GetPFMList <- function(seurat_obj) {
  seurat_obj@misc$motifs$pfm_list
}


############################
# TF Target Genes:
###########################

#' SetMotifTargets
#'
#' @param seurat_obj A Seurat object
#' @param motif_targets list of motifs and their target genes
#' @keywords scRNA-seq
#' @export
#' @examples SetMotifTargets
SetMotifTargets <- function(seurat_obj, motif_targets) {
  # make a spot for the motif info if it's not already there:
  if (is.null(seurat_obj@misc$motifs)) {
    seurat_obj@misc$motifs <- list()
  }
  seurat_obj@misc$motifs$motif_targets <- motif_targets
  seurat_obj
}


#' GetMotifTargets
#'
#' @param seurat_obj A Seurat object
#' @keywords scRNA-seq
#' @export
#' @examples GetMotifTargets
GetMotifTargets <- function(seurat_obj) {
  seurat_obj@misc$motifs$motif_targets
}


############################
# motif overlap
###########################


#' SetMotifOverlap
#'
#' @param seurat_obj A Seurat object
#' @param overlap_df dataframe containing motif-module overlap info
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMotifOverlap <- function(seurat_obj, overlap_df, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$motif_module_overlaps <- overlap_df
  seurat_obj
}


#' GetMotifOverlap
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMotifOverlap <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$motif_module_overlaps
}


############################
# motif scores
###########################


#' SetMotifScores
#'
#' @param seurat_obj A Seurat object
#' @param tf_scores dataframe of tf motif target scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetMotifScores <- function(seurat_obj, tf_scores, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$motif_target_scores <- tf_scores
  seurat_obj
}


#' GetMotifScores
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetMotifScores <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$motif_target_scores
}

############################
# ModuleUMAP
###########################

#' SetModuleUMAP
#'
#' @param seurat_obj A Seurat object
#' @param umap_df dataframe of UMAP coordinates
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModuleUMAP <- function(seurat_obj, umap_df, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_umap <- umap_df
  seurat_obj
}

#' GetModuleUMAP
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetModuleUMAP <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$module_umap
}

############################
# ModuleTraitCorrelation
###########################


#' SetModuleTraitCorrelation
#'
#' @param seurat_obj A Seurat object
#' @param mt_cor matrix of module-trait correlation results
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModuleTraitCorrelation <- function(seurat_obj, mt_cor, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$mt_cor <- mt_cor
  seurat_obj
}

#' GetModuleTraitCorrelation
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetModuleTraitCorrelation <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$mt_cor
}

############################
# ModulePreservation
###########################


#' SetModulePreservation
#'
#' @param seurat_obj A Seurat object
#' @param mt_cor matrix of module-trait correlation results
#' @param mod_name name of the module preservation test to store
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetModulePreservation <- function(seurat_obj, mod_pres, mod_name, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # make an empty list if module preservation hasn't been called yet
  if (is.null(seurat_obj@misc[[wgcna_name]]$module_preservation)) {
    seurat_obj@misc[[wgcna_name]]$module_preservation <- list()
  }

  seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]] <- mod_pres
  seurat_obj
}



#' GetModulePreservation
#'
#' @param seurat_obj A Seurat object
#' @param mod_name name of the module preservation test to store
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetModulePreservation <- function(seurat_obj, mod_name, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  if (is.null(seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]])) {
    stop("Invalid module preservation name.")
  }
  seurat_obj@misc[[wgcna_name]]$module_preservation[[mod_name]]
}




#' SetRegulonScores
#'
#' @param seurat_obj A Seurat object
#' @param regulon_scores dataframe storing the TF regulon scores
#' @param target_type dataframe storing the TF regulon scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetRegulonScores <- function(seurat_obj, regulon_scores, target_type, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)


  # if regulon scores have not been set, make a list to store them
  if (is.null(seurat_obj@misc[[wgcna_name]]$regulon_scores)) {
    tmp <- list(regulon_scores)
    names(tmp) <- target_type
    seurat_obj@misc[[wgcna_name]]$regulon_scores <- tmp
  } else {
    seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]] <- regulon_scores
  }
  seurat_obj
}

#' GetRegulonScores
#'
#' @param seurat_obj A Seurat object
#' @param target_type dataframe storing the TF regulon scores
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetRegulonScores <- function(seurat_obj, target_type, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # get the regulon scores
  seurat_obj@misc[[wgcna_name]]$regulon_scores[[target_type]]
}


############################
# getters and setters
###########################

#' SetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param tf_net dataframe storing the TF network info in ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFNetwork <- function(seurat_obj, tf_net, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_net <- tf_net
  seurat_obj
}

#' GetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFNetwork <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_net
}


#' SetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param tf_eval dataframe storing the TF network evaluation info from ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFEval <- function(seurat_obj, tf_eval, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_eval <- tf_eval
  seurat_obj
}

#' GetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFEval <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_eval
}

#' SetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param tf_regulons dataframe storing the TF regulon info from AssignTFRegulons
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFRegulons <- function(seurat_obj, tf_regulons, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_regulons <- tf_regulons
  seurat_obj
}

#' GetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFRegulons <- function(seurat_obj, wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_regulons
}


############################
# Reset module names:
###########################
