# LUAD tumor vs. normal explorer (reads results/luad.sqlite)
# Run from repo root: Rscript -e 'shiny::runApp("app")'
library(shiny)
library(DBI)
library(RSQLite)
library(ggplot2)

db_path <- if (file.exists("results/luad.sqlite")) "results/luad.sqlite" else "../results/luad.sqlite"
con <- dbConnect(SQLite(), db_path)
onStop(function() dbDisconnect(con))

has <- function(t) t %in% dbListTables(con)
expr_genes <- if (has("expression")) dbGetQuery(con, "SELECT DISTINCT gene_name FROM expression ORDER BY gene_name")$gene_name else character()

ui <- fluidPage(
  titlePanel("TCGA-LUAD: tumor vs. matched normal"),
  sidebarLayout(
    sidebarPanel(width = 3,
      sliderInput("fdr", "FDR cutoff", min = 0.001, max = 0.1, value = 0.05, step = 0.001),
      sliderInput("lfc", "|log2 fold change| >", min = 0, max = 5, value = 1, step = 0.25),
      selectizeInput("gene", "Gene", choices = expr_genes,
                     selected = if (length(expr_genes)) expr_genes[1] else NULL)
    ),
    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Differential expression",
                 textOutput("n_sig"), plotOutput("volcano", height = 420),
                 tableOutput("de_table")),
        tabPanel("Gene view", plotOutput("paired", height = 420), tableOutput("gene_stats")),
        tabPanel("Pathways", tableOutput("gsea")),
        tabPanel("Survival", tableOutput("surv"))
      )
    )
  )
)

server <- function(input, output, session) {
  de <- reactive({
    dbGetQuery(con, "SELECT gene_name, baseMean, log2FC, padj FROM de_results WHERE padj IS NOT NULL")
  })
  sig <- reactive({
    dbGetQuery(con,
      "SELECT gene_name, ROUND(log2FC, 2) AS log2FC, padj
       FROM de_results WHERE padj < ? AND ABS(log2FC) > ?
       ORDER BY padj LIMIT 100",
      params = list(input$fdr, input$lfc))
  })

  output$n_sig <- renderText({
    n <- dbGetQuery(con, "SELECT COUNT(*) n FROM de_results WHERE padj < ? AND ABS(log2FC) > ?",
                    params = list(input$fdr, input$lfc))$n
    paste(n, "genes pass the current cutoffs (top 100 shown below)")
  })

  output$volcano <- renderPlot({
    d <- de()
    d$sig <- d$padj < input$fdr & abs(d$log2FC) > input$lfc
    ggplot(d, aes(log2FC, -log10(padj), color = sig)) +
      geom_point(size = .6, alpha = .6) +
      scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#c0392b"), guide = "none") +
      geom_point(data = subset(d, gene_name == input$gene), color = "black", size = 3) +
      theme_bw()
  })

  output$de_table <- renderTable(sig(), digits = 4)

  output$paired <- renderPlot({
    req(input$gene)
    d <- dbGetQuery(con, "SELECT patient, condition, vst FROM expression WHERE gene_name = ?",
                    params = list(input$gene))
    d$condition <- factor(d$condition, c("Normal", "Tumor"))
    ggplot(d, aes(condition, vst)) +
      geom_line(aes(group = patient), color = "grey75") +
      geom_boxplot(aes(fill = condition), width = .35, alpha = .6, outlier.shape = NA) +
      geom_point(size = 1) +
      scale_fill_manual(values = c(Normal = "#2c6fbb", Tumor = "#c0392b"), guide = "none") +
      labs(y = "VST expression", x = NULL, title = paste(input$gene, "— lines join the same patient")) +
      theme_bw(base_size = 14)
  })

  output$gene_stats <- renderTable({
    req(input$gene)
    dbGetQuery(con, "SELECT gene_name, gene_id, baseMean, log2FC, padj FROM de_results WHERE gene_name = ?",
               params = list(input$gene))
  }, digits = 4)

  output$gsea <- renderTable({
    if (!has("gsea_hallmarks")) return(NULL)
    dbGetQuery(con, "SELECT pathway, ROUND(NES, 2) AS NES, padj, size FROM gsea_hallmarks ORDER BY padj")
  }, digits = 4)

  output$surv <- renderTable({
    if (!has("survival_cox")) return(NULL)
    dbGetQuery(con, "SELECT gene_name, ROUND(HR, 2) AS HR, ROUND(HR_low, 2) AS HR_low,
                            ROUND(HR_high, 2) AS HR_high, p, fdr, ROUND(HR_adj, 2) AS HR_age_stage, fdr_adj_model
                     FROM survival_cox ORDER BY p")
  }, digits = 4)
}

shinyApp(ui, server)
