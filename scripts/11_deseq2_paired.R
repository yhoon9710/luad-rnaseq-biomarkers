# Paired tumor vs. normal differential expression on TCGA-LUAD
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
})
set.seed(42)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

se <- readRDS("data/tcga_luad_se.rds")

# ---- Sample selection: patients with both a primary tumor and a normal ----
cd <- as.data.frame(colData(se))
cd$barcode   <- colnames(se)
cd$patient   <- substr(cd$barcode, 1, 12)
cd$condition <- ifelse(cd$sample_type == "Primary Tumor", "Tumor", "Normal")

paired_patients <- cd %>%
  distinct(patient, condition) %>%
  count(patient) %>% filter(n == 2) %>% pull(patient)

keep <- cd %>%
  filter(patient %in% paired_patients) %>%
  arrange(barcode) %>%
  distinct(patient, condition, .keep_all = TRUE)   # one aliquot per patient/tissue

se_p <- se[, keep$barcode]
# protein-coding genes only (keeps the multiple-testing burden reasonable)
se_p <- se_p[rowData(se_p)$gene_type == "protein_coding", ]

coldata <- data.frame(
  row.names = keep$barcode,
  patient   = factor(keep$patient),
  condition = factor(keep$condition, levels = c("Normal", "Tumor"))
)
counts <- assay(se_p, "unstranded")
message("Matched pairs: ", length(paired_patients), " | genes: ", nrow(counts))

# ---- DESeq2, paired design ----
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ patient + condition)
dds <- dds[rowSums(counts(dds) >= 10) >= length(paired_patients), ]
dds <- DESeq(dds)

res  <- results(dds, name = "condition_Tumor_vs_Normal", alpha = 0.05)
shr  <- lfcShrink(dds, coef = "condition_Tumor_vs_Normal", type = "apeglm")

gene_ids   <- sub("\\..*$", "", rownames(dds))
gene_names <- rowData(se_p)[rownames(dds), "gene_name"]

de <- tibble(
  gene_id    = gene_ids,
  gene_name  = gene_names,
  baseMean   = res$baseMean,
  log2FC     = shr$log2FoldChange,   # shrunken, for ranking/plots
  log2FC_mle = res$log2FoldChange,
  stat       = res$stat,             # for GSEA ranking
  pvalue     = res$pvalue,
  padj       = res$padj
) %>% arrange(padj)

write_tsv(de, "results/tables/de_tumor_vs_normal.tsv")
sig <- de %>% filter(padj < 0.05, abs(log2FC) > 1)
message("Significant (FDR<0.05, |log2FC|>1): ", nrow(sig),
        " (up ", sum(sig$log2FC > 0), ", down ", sum(sig$log2FC < 0), ")")

# ---- Figures ----
vsd <- vst(dds, blind = TRUE)

pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 2, alpha = .8) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
       title = "TCGA-LUAD matched tumor / normal") +
  theme_bw()
ggsave("results/figures/pca.png", p, width = 6, height = 4.5, dpi = 300)

lab <- bind_rows(slice_max(sig, log2FC, n = 10), slice_min(sig, log2FC, n = 10))
p <- de %>% filter(!is.na(padj)) %>%
  mutate(class = case_when(padj < .05 & log2FC > 1 ~ "Up in tumor",
                           padj < .05 & log2FC < -1 ~ "Down in tumor",
                           TRUE ~ "NS")) %>%
  ggplot(aes(log2FC, -log10(pmax(padj, 1e-300)), color = class)) +
  geom_point(size = .6, alpha = .6) +
  geom_text_repel(data = lab, aes(label = gene_name), color = "black",
                  size = 3, max.overlaps = 30) +
  scale_color_manual(values = c("Up in tumor" = "#c0392b",
                                "Down in tumor" = "#2c6fbb", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(.05), linetype = 2) +
  labs(title = "Paired DESeq2: tumor vs normal", color = NULL) +
  theme_bw()
ggsave("results/figures/volcano.png", p, width = 7, height = 5.5, dpi = 300)

top <- sig %>% slice_min(padj, n = 50, with_ties = FALSE)
idx <- match(top$gene_id, gene_ids)
mat <- assay(vsd)[idx, ]
rownames(mat) <- top$gene_name
mat <- t(scale(t(mat)))
ann <- data.frame(row.names = colnames(mat), condition = coldata$condition)
png("results/figures/heatmap_top50.png", width = 2400, height = 2600, res = 300)
pheatmap(mat, annotation_col = ann, show_colnames = FALSE, fontsize_row = 6,
         main = "Top 50 DE genes (z-scored VST)")
dev.off()

# ---- Expression of top 500 DE genes, long format, for the Shiny app ----
top500 <- de %>% filter(!is.na(padj)) %>% slice_min(padj, n = 500, with_ties = FALSE)
idx <- match(top500$gene_id, gene_ids)
expr <- as.data.frame(assay(vsd)[idx, ])
expr$gene_name <- top500$gene_name
expr %>%
  pivot_longer(-gene_name, names_to = "barcode", values_to = "vst") %>%
  left_join(tibble(barcode = rownames(coldata),
                   patient = as.character(coldata$patient),
                   condition = as.character(coldata$condition)), by = "barcode") %>%
  write_tsv("results/tables/expr_top500_pairs.tsv")

saveRDS(dds, "results/dds_paired.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_deseq2.txt")
