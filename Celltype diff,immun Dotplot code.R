# 1. 엑셀 파일에서 모든 시트 읽기
sheets <- excel_sheets("/home/welcome1/GSEA files/GSEA_results_(both)DXvsLNLD_Res.xlsx")
sheets <- c("HSC", "GMP", "ERP", "LMPP", "Early Erythroid", "Late Erythroid 1", "Late Erythroid 2", "T cell", "CD4 T cell", "CD8 T cell", "Regulatory T cell", "B cell", "Plasma cell", "NK cell", "CD14 Monocyte", "Neutrophil", "Macrophage")
go_data <- lapply(sheets, function(sheet) {
  read_excel("/home/welcome1/GSEA files/GSEA_results_(both)DXvsLNLD_Res.xlsx", sheet = sheet, skip = 1) %>%
    mutate(Celltype = sheet)
})

combined <- bind_rows(go_data)

combined$NES <- abs(combined$NES)
combined$`-log2(FDR q-val)` <- -log2(combined$`FDR q-val`)
combined$`-log2(FDR q-val)`[is.infinite(combined$`-log2(FDR q-val)`)] <- 10

combined$`-log2(FDR q-val)`[combined$VS == "LNLD_ResvsDX"] <- 
  -abs(combined$`-log2(FDR q-val)`[combined$VS == "LNLD_ResvsDX"])

# 2. 관심 있는 GO term만 필터링
selected_terms <- c("GOBP_ALPHA_BETA_T_CELL_DIFFERENTIATION", "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION","GOBP_REGULATION_OF_T_HELPER_17_CELL_DIFFERENTIATION", "GOBP_NATURAL_KILLER_CELL_DIFFERENTIATON", "GOBP_ALPHA_BETA_T_CELL_PROLIFERATION", "GOBP_B_CELL_DIFFERENTIATION", "GOBP_MONOCYTE_DIFFERENTIATION", "GOBP_MYELOID_DENDRITIC_CELL_DIFFERENTIATION",  "GOBP_NEUTROPHIL_DIFFERENTIATION", "GOBP_MYELOID_CELL_DIFFERENTIATION", "GOBP_REGULATION_OF_T_HELPER_17_TYPE_IMMUNE_RESPONSE", "GOBP_T_CELL_ACTIVATION", "GOBP_ALPHA_BETA_T_CELL_ACTIVATION" , "GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY", "GOBP_REGULATION_OF_NATURAL_KILLER_CELL_ACTIVATION", "GOBP_POSITIVE_REGULATION_OF_B_CELL_ACTIVATION", "GOBP_B_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_MACROPHAGE_ACTIVATION", "GOBP_MACROPHAGE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION")
filtered <- combined %>%
  filter(GS %in% selected_terms)

# 4. factor로 변환하면서 순서 지정
filtered$GS <- factor(filtered$GS, levels = rev(selected_terms))
filtered$Celltype <- factor(filtered$Celltype, levels = sheets)

# 3. DotPlot
ggplot(filtered, aes(x = GS, y = Celltype)) +
  geom_point(aes(size = NES, color = `-log2(FDR q-val)`)) +
  scale_color_gradientn(colors = c("red", "skyblue", "blue"),
                        name = "expression score") +
  scale_size(range = c(1, 7), name = "NES") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1)
  ) +
  coord_flip()

  