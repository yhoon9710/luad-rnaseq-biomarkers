# Download TCGA-LUAD STAR gene counts (primary tumors + solid tissue normals)
# ~600 files, a few GB. Run once; result cached in data/tcga_luad_se.rds
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

out <- "data/tcga_luad_se.rds"
if (file.exists(out)) { message("Already downloaded: ", out); quit(save = "no") }

query <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)

res <- getResults(query)
message("Files: ", nrow(res), " | ", paste(names(table(res$sample_type)),
        table(res$sample_type), sep = "=", collapse = ", "))

GDCdownload(query, directory = "data/GDCdata", files.per.chunk = 20)
se <- GDCprepare(query, directory = "data/GDCdata")

saveRDS(se, out)
message("Saved ", out, ": ", nrow(se), " genes x ", ncol(se), " samples")
