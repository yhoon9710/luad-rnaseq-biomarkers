# Hallmark gene set enrichment (fgsea) on the paired DE statistic
suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
})
set.seed(42)

de <- read_tsv("results/tables/de_tumor_vs_normal.tsv", show_col_types = FALSE)

ranks <- de %>%
  filter(!is.na(stat), !is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
ranks <- deframe(ranks[, c("gene_name", "stat")])

# msigdbr >= 10 uses `collection`, older versions use `category`
hm <- if ("collection" %in% names(formals(msigdbr))) {
  msigdbr(species = "Homo sapiens", collection = "H")
} else {
  msigdbr(species = "Homo sapiens", category = "H")
}
pathways <- split(hm$gene_symbol, hm$gs_name)

gsea <- fgsea(pathways, ranks, minSize = 15, maxSize = 500) %>%
  as_tibble() %>% arrange(padj) %>%
  mutate(leadingEdge = map_chr(leadingEdge, ~ paste(head(.x, 25), collapse = ",")),
         pathway = str_remove(pathway, "^HALLMARK_"))

write_tsv(gsea, "results/tables/gsea_hallmarks.tsv")
message("Hallmarks at FDR<0.05: ", sum(gsea$padj < .05, na.rm = TRUE))

plot_df <- gsea %>% filter(padj < .05) %>%
  group_by(direction = sign(NES)) %>% slice_max(abs(NES), n = 10) %>% ungroup()
p <- ggplot(plot_df, aes(NES, reorder(pathway, NES), fill = NES > 0)) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#2c6fbb"),
                    labels = c(`TRUE` = "Up in tumor", `FALSE` = "Down in tumor"),
                    name = NULL) +
  labs(x = "Normalized enrichment score", y = NULL,
       title = "MSigDB Hallmarks: LUAD tumor vs matched normal") +
  theme_bw()
ggsave("results/figures/gsea_hallmarks.png", p, width = 7.5, height = 5.5, dpi = 300)
