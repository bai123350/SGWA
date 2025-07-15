# SGWA: Single-cell Gene Weighted-coexpression Analysis

SGWA 是一个R包，旨在简化对Seurat对象进行的加权基因共表达网络分析 (WGCNA)。它将WGCNA的强大功能与单细胞分析工具Seurat无缝集成，提供了一个从元细胞构建(metacell construction)、网络构建、模块识别到下游功能富集和可视化的完整工作流程。

## 主要功能

*   **与Seurat无缝集成**: 直接在Seurat对象上操作，无需繁琐的数据转换。
*   **元细胞(Metacell)构建**: 通过将相似的细胞分组来处理单细胞数据的稀疏性。
*   **自动化WGCNA流程**: 使用`run_sgwgcna_pipeline`函数一键运行核心分析流程。
*   **丰富的功能模块**:
    *   差异模块特征基因 (DME) 分析
    *   核心基因 (Hub Gene) 识别与网络可视化
    *   模块与细胞表型/分组的关联分析
*   **多样化的可视化**: 提供一系列绘图功能，如树状图、热图、网络图等，以直观展示分析结果。

## 安装

您可以使用 `devtools` 从GitHub安装此包的开发版本：

```r
# 如果您尚未安装devtools，请先安装
# install.packages("devtools")

devtools::install_github("bai123350/SGWA")
```

## 快速入门

这是一个在Seurat对象上运行SGWA主要流程的最小示例。

```r
# 加载库
library(SGWA)
library(Seurat)

# 假设 'seurat_obj' 是您已经预处理过的Seurat对象
# (例如，已完成标准化、降维和聚类)

# 运行完整的WGCNA分析流程
# 该函数将会在seurat_obj中计算并储存WGCNA结果
seurat_obj <- run_sgwgcna_pipeline(seurat_obj,gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "ASC_consensus",
    metacell_group.by = c("celltype", "sample"),
    metacell_ident.group1 = "celltype",
    metacell_ident.group2 = "sample",
    k = 25,
    max_shared = 12,
    min_cells = 50,
    target_metacells = 250,
    reduction = "harmony",
    dat_expr_group_name = "B cell",
    network_type = "signed",
    n_hubs = 20,
    n_genes_score = 25,
    assay = NULL,
    meta_yazu = "sample",
    dot_plot_group = "sample",
    plot_path = "./plots/",
    plot_height = 10,
    plot_width = 10,
    cols = 6)
```

## 依赖

本包依赖于以下核心R包：
- [Seurat](https://satijalab.org/seurat/)
- [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/)
- [UCell](https://github.com/carmonalab/UCell)
- [harmony](https://github.com/immunogenomics/harmony)
- [enrichR](https://github.com/wiflish/enrichR)
- [ggplot2](https://ggplot2.tidyverse.org)
- [dplyr](https.dplyr.tidyverse.org)

## 许可证

本项目采用 MIT 许可证。详情请见 [LICENSE](LICENSE) 文件。
