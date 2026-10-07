# Salmon transcript quantifications -> gene-level counts (tximport)
suppressPackageStartupMessages({
  library(tximport)
  library(readr)
})

samples <- read.delim("config/samples.tsv", comment.char = "#")
files <- file.path("results/salmon", samples$sample, "quant.sf")
names(files) <- samples$sample
stopifnot(all(file.exists(files)))

# GENCODE headers: ENST|ENSG|OTTHUMG|OTTHUMT|tx_name|gene_name|length|biotype|
tx_ids <- read_tsv(files[1], show_col_types = FALSE)$Name
parts <- strsplit(tx_ids, "|", fixed = TRUE)
tx2gene <- data.frame(
  tx        = tx_ids,
  gene_id   = sub("\\..*$", "", vapply(parts, `[`, "", 2)),
  gene_name = vapply(parts, `[`, "", 6)
)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("tx", "gene_id")])

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
counts <- data.frame(gene_id = rownames(txi$counts), round(txi$counts), check.names = FALSE)
genes  <- unique(tx2gene[, c("gene_id", "gene_name")])
counts <- merge(genes, counts, by = "gene_id")
write_tsv(counts, "results/tables/fastq_gene_counts.tsv")
saveRDS(txi, "results/txi.rds")

message("Genes: ", nrow(counts), " | samples: ", ncol(txi$counts))
message("Mapping rates are in results/qc/multiqc_report.html")
