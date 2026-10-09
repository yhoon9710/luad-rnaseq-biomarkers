# Personal follow-up: genes from my undergraduate antiviral-immunity research
# (Yang Lab, Biola) -- TRIM25, ZAP (gene symbol ZC3HAV1), IGF2BP3 -- in TCGA-LUAD.
#   1. Are they differentially expressed in tumor vs. matched normal? (paired DESeq2 results)
#   2. Per-patient paired expression plot
#   3. Are they associated with overall survival across all primary tumors? (Cox, adjusted for age + stage)
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(tidyverse)
  library(survival)
  library(survminer)
})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

genes <- c("TRIM25", "ZC3HAV1", "IGF2BP3")   # ZC3HAV1 = ZAP
pretty <- c(TRIM25 = "TRIM25", ZC3HAV1 = "ZAP (ZC3HAV1)", IGF2BP3 = "IGF2BP3")

# ---- 1. Paired DE results -------------------------------------------------
de <- read_tsv("results/tables/de_tumor_vs_normal.tsv", show_col_types = FALSE)
de_sub <- de %>%
  filter(gene_name %in% genes) %>%
  mutate(fold_change = round(2^log2FC, 2)) %>%
  select(gene_name, gene_id, baseMean, log2FC, fold_change, padj)

cat("\n== Tumor vs matched normal (58 pairs, paired DESeq2) ==\n")
print(as.data.frame(de_sub %>% select(-gene_id)), digits = 3)
write_tsv(de_sub, "results/tables/yanglab_genes_de.tsv")

# ---- 2. Paired expression plot ---------------------------------------------
dds <- readRDS("results/dds_paired.rds")
nc  <- counts(dds, normalized = TRUE)
gid <- sub("\\..*$", "", rownames(dds))
cd  <- as.data.frame(colData(dds))

plot_df <- map_dfr(seq_len(nrow(de_sub)), function(i) {
  idx <- match(de_sub$gene_id[i], gid)
  tibble(gene      = de_sub$gene_name[i],
         patient   = as.character(cd$patient),
         condition = factor(cd$condition, c("Normal", "Tumor")),
         expr      = log2(nc[idx, ] + 1))
}) %>%
  left_join(de_sub %>%
              transmute(gene = gene_name,
                        label = sprintf("%s\nlog2FC = %.2f, FDR = %.1e",
                                        pretty[gene_name], log2FC, padj)),
            by = "gene")

p <- ggplot(plot_df, aes(condition, expr)) +
  geom_line(aes(group = patient), color = "grey75", linewidth = .3) +
  geom_boxplot(aes(fill = condition), width = .4, alpha = .6, outlier.shape = NA) +
  geom_point(size = .8) +
  facet_wrap(~ label, scales = "free_y") +
  scale_fill_manual(values = c(Normal = "#2c6fbb", Tumor = "#c0392b"), guide = "none") +
  labs(x = NULL, y = "log2(normalized counts + 1)",
       title = "Yang Lab antiviral genes in LUAD: tumor vs matched normal",
       subtitle = "Grey lines connect the tumor and normal sample from the same patient") +
  theme_bw()
ggsave("results/figures/yanglab_genes_paired.png", p, width = 10, height = 4.5, dpi = 300)

# ---- 3. Survival across all primary tumors ---------------------------------
se <- readRDS("data/tcga_luad_se.rds")
tum <- as.data.frame(colData(se)) %>%
  mutate(barcode = colnames(se), patient = substr(barcode, 1, 12)) %>%
  filter(sample_type == "Primary Tumor") %>%
  arrange(barcode) %>% distinct(patient, .keep_all = TRUE)

# survival data from the GDC patient table (run 05_get_clinical.R first)
clin_all <- read_tsv("data/clinical_luad.tsv", show_col_types = FALSE)
clin <- tum %>%
  select(barcode, patient) %>%
  inner_join(clin_all %>% select(patient, event, time = time_months, age, stage),
             by = "patient")
cat("\nTumors with survival data:", nrow(clin), "| deaths:", sum(clin$event),
    "| alive (censored):", sum(clin$event == 0), "\n")

tpm <- assay(se, "tpm_unstrand")
rn  <- rowData(se)$gene_name

cox <- map_dfr(genes, function(g) {
  x <- log2(tpm[match(g, rn), clin$barcode] + 1)
  d <- mutate(clin, expr = as.numeric(scale(x)))
  uni <- summary(coxph(Surv(time, event) ~ expr, data = d))
  adj <- summary(coxph(Surv(time, event) ~ expr + age + stage, data = d))

  d$group <- factor(ifelse(x > median(x), "High", "Low"), c("Low", "High"))
  fit <- survfit(Surv(time, event) ~ group, data = d)
  km  <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, xlab = "Months",
                    title = paste0(pretty[g], " expression in LUAD tumors (median split)"),
                    palette = c("#2c6fbb", "#c0392b"))
  png(paste0("results/figures/km_yanglab_", g, ".png"), width = 2000, height = 2000, res = 300)
  print(km)
  dev.off()

  tibble(gene = pretty[g], n = nrow(d), deaths = sum(d$event),
         HR_per_SD = uni$conf.int["expr", 1],
         CI_low = uni$conf.int["expr", 3], CI_high = uni$conf.int["expr", 4],
         p = uni$coefficients["expr", 5],
         HR_adj_age_stage = adj$conf.int["expr", 1],
         p_adj_age_stage  = adj$coefficients["expr", 5])
})

cat("\n== Overall survival, all primary tumors (Cox; HR per 1 SD higher expression) ==\n")
print(as.data.frame(cox), digits = 3)
write_tsv(cox, "results/tables/yanglab_genes_survival.tsv")
cat("\nHR > 1: higher expression = worse survival; HR < 1: better survival.\n",
    "Exploratory: 3 genes chosen in advance, no independent validation cohort.\n", sep = "")
