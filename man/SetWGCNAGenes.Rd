% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters_and_setters.R
\name{SetWGCNAGenes}
\alias{SetWGCNAGenes}
\title{SetWGCNAGenes}
\usage{
SetWGCNAGenes(seurat_obj, gene_list, wgcna_name = NULL)
}
\arguments{
\item{seurat_obj}{一个 Seurat 对象}

\item{gene_list}{用于 WGCNA 的基因向量}

\item{wgcna_name}{seurat_obj@misc 插槽中 hdWGCNA 实验的名称}
}
\description{
SetWGCNAGenes
}
\keyword{scRNA-seq}
