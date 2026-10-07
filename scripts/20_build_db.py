"""Load analysis results into a SQLite database used by the Shiny explorer."""
import sqlite3
from pathlib import Path

import pandas as pd

TABLES = Path("results/tables")
DB = Path("results/luad.sqlite")

sources = {
    "de_results": "de_tumor_vs_normal.tsv",
    "gsea_hallmarks": "gsea_hallmarks.tsv",
    "survival_cox": "survival_cox_top25.tsv",
    "expression": "expr_top500_pairs.tsv",
}

DB.unlink(missing_ok=True)
with sqlite3.connect(DB) as con:
    for table, fname in sources.items():
        path = TABLES / fname
        if not path.exists():
            print(f"[skip] {fname} not found")
            continue
        df = pd.read_csv(path, sep="\t")
        df.to_sql(table, con, index=False)
        print(f"[load] {table}: {len(df):,} rows")

    con.executescript(
        """
        CREATE INDEX IF NOT EXISTS ix_de_gene   ON de_results(gene_name);
        CREATE INDEX IF NOT EXISTS ix_de_padj   ON de_results(padj);
        CREATE INDEX IF NOT EXISTS ix_expr_gene ON expression(gene_name);
        """
    )

    # sanity check: a query you might be asked to write in an interview
    q = """
        SELECT CASE WHEN log2FC > 0 THEN 'up' ELSE 'down' END AS direction,
               COUNT(*) AS n_genes
        FROM de_results
        WHERE padj < 0.05 AND ABS(log2FC) > 1
        GROUP BY direction
    """
    print(pd.read_sql(q, con).to_string(index=False))

print(f"Wrote {DB}")
