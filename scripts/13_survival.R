# Overall-survival association of the top DE genes across ALL TCGA-LUAD primary tumors
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(tidyverse)
  library(survival)
  library(survminer)
})

se <- readRDS("data/tcga_luad_se.rds")
de <- read_tsv("results/tables/de_tumor_vs_normal.tsv", show_col_types = FALSE)

# one primary tumor per patient
cd <- as.data.frame(colData(se)) %>%
  mutate(barcode = colnames(se), patient = substr(barcode, 1, 12)) %>%
  filter(sample_type == "Primary Tumor") %>%
  arrange(barcode) %>% distinct(patient, .keep_all = TRUE)

# survival data from the GDC patient table (run 05_get_clinical.R first)
clin_all <- read_tsv("data/clinical_luad.tsv", show_col_types = FALSE)
clin <- cd %>%
  select(barcode, patient) %>%
  inner_join(clin_all %>% select(patient, event, time = time_months, age, stage),
             by = "patient")
message("Tumors with survival data: ", nrow(clin), " | deaths: ", sum(clin$event),
        " | alive (censored): ", sum(clin$event == 0))

# variance-stabilized expression of tumors
se_t <- se[rowData(se)$gene_type == "protein_coding", clin$barcode]
dds  <- DESeqDataSetFromMatrix(assay(se_t, "unstranded"),
                               data.frame(row.names = clin$barcode, x = rep(1, nrow(clin))), ~ 1)
vsd  <- vst(dds, blind = TRUE)
gid  <- sub("\\..*$", "", rownames(vsd))

top <- de %>% filter(padj < .05, abs(log2FC) > 1) %>%
  slice_min(padj, n = 25, with_ties = FALSE)

cox <- map_dfr(seq_len(nrow(top)), function(i) {
  x <- as.numeric(scale(assay(vsd)[match(top$gene_id[i], gid), ]))
  d <- mutate(clin, expr = x)
  uni <- summary(coxph(Surv(time, event) ~ expr, data = d))
  adj <- summary(coxph(Surv(time, event) ~ expr + age + stage, data = d))
  tibble(gene_name = top$gene_name[i], gene_id = top$gene_id[i],
         de_log2FC = top$log2FC[i],
         HR = uni$conf.int["expr", 1], HR_low = uni$conf.int["expr", 3],
         HR_high = uni$conf.int["expr", 4], p = uni$coefficients["expr", 5],
         HR_adj = adj$conf.int["expr", 1], p_adj_model = adj$coefficients["expr", 5])
}) %>%
  mutate(fdr = p.adjust(p, "BH"), fdr_adj_model = p.adjust(p_adj_model, "BH")) %>%
  arrange(p)

write_tsv(cox, "results/tables/survival_cox_top25.tsv")
message("Genes with FDR<0.05 (univariate): ", sum(cox$fdr < .05),
        " | adjusted for age+stage: ", sum(cox$fdr_adj_model < .05))

# KM plots (median split, visualization only) for the 3 strongest genes
for (g in head(cox$gene_name, 3)) {
  x <- assay(vsd)[match(cox$gene_id[cox$gene_name == g], gid), ]
  d <- mutate(clin, group = factor(ifelse(x > median(x), "High", "Low"), c("Low", "High")))
  fit <- survfit(Surv(time, event) ~ group, data = d)
  km <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE,
                   xlab = "Months", title = paste0(g, " expression (median split)"),
                   palette = c("#2c6fbb", "#c0392b"))
  png(paste0("results/figures/km_", g, ".png"), width = 2000, height = 2000, res = 300)
  print(km)
  dev.off()
}
