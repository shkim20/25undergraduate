# 1. 대상 gene set 이름 목록
target_gs_list <- c(
  "GOBP_T_HELPER_17_CELL_DIFFERENTIATION",
  "GOBP_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION",
  "GOBP_NEGATIVE_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION"
)

# 2. .gmt 파일 경로
gmt_files <- list(
  HSC = "/home/welcome1/gene_sets(HSC).gmt",
  GMP = "/home/welcome1/gene_sets(GMP).gmt",
  ERP = "/home/welcome1/gene_sets(ERP).gmt"
)

# 3. DotPlot을 저장할 리스트
dotplots <- list()

# 4. 각 gene set에 대해 반복
all_genes <- c()

for (gs_name in target_gs_list) {
  gene_lists <- lapply(gmt_files, function(f) {
    gmt <- getGmt(f)
    geneIds(gmt[[gs_name]])
  })
  common_genes <- Reduce(intersect, gene_lists)
  genes_to_use <- intersect(common_genes, rownames(MDS_Progen_LNLD))
  
  # 핵심: gene 이름에 gene set 이름 prefix 붙이기
  genes_to_use_renamed <- setNames(genes_to_use, paste0(gs_name, "_", genes_to_use))
  all_genes <- c(all_genes, genes_to_use_renamed)
}

DotPlot(MDS_ERP, features = Tcell_gene_GMP, group.by = "SampleGroup") +
  scale_color_gradientn(
    colors = c("skyblue", "white", "red"),
    name = "Avg. Expression"
  ) +
  scale_size(range = c(1, 7), name = "Percent Expressed") +
  ggtitle("T cell-specific genes expression in ERP") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

all_genes
