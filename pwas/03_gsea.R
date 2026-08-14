################################################################################
#
#  Pathway enrichment of the PWAS results — panel-restricted GSEA
#
#  The Olink Explore panel assays a curated subset of the proteome, so a pathway
#  that happens to be well covered by the panel would look enriched relative to
#  one that is barely represented. Enrichment is therefore evaluated against the
#  assayed proteins rather than against the whole annotation database.
#
#  fgsea takes no `universe` argument, so the restriction is applied on the
#  pathway side: each gene set is intersected with the assayed protein set
#  BEFORE fgsea runs, and the minSize/maxSize filter is re-applied afterwards so
#  that sets left too small by the intersection are dropped rather than tested
#  on a handful of genes. `restrict_to_panel()` implements this.
#
#  Input   : input/assoc_results/pwas_primary.csv
#  Output  : output/pwas_gsea/panel_restricted/
#            tables/panel_universe_entrez.csv is consumed by scripts 02 and 03.
#
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(fgsea)
})

cat("================================================================================\n")
cat("Panel-restricted GSEA of the primary PWAS model\n")
cat("================================================================================\n\n")

# =============================================================================
# 0. Paths
# =============================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

pwas_file  <- "input/assoc_results/pwas_primary.csv"
output_dir <- "output/pwas_gsea/panel_restricted"
dir.create(output_dir,                            recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"),       recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"),      recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Input file : %s\n", pwas_file))
cat(sprintf("Output dir : %s\n", output_dir))

# =============================================================================
# 1. PWAS results -> rank vector
# =============================================================================

cat("\n---- Step 1: PWAS results + rank metric ----\n\n")

results <- fread(pwas_file)
results <- results[converged == TRUE & is.finite(beta) & is.finite(se) & se > 0 & is.finite(pval)]

cat(sprintf("Converged proteins        : %d\n", nrow(results)))
cat(sprintf("Bonferroni-sig (pval_bonf<.05): %d\n", sum(results$pval_bonf < 0.05, na.rm = TRUE)))
cat(sprintf("FDR-sig (pval_fdr<.05)        : %d\n", sum(results$pval_fdr  < 0.05, na.rm = TRUE)))

results[, gene_symbol     := toupper(protein)]
results[, rank_signed_logp := sign(beta) * (-log10(pval))]
results[, rank_t           := beta / se]
results[, rank_z           := qnorm(1 - pval/2) * sign(beta)]
results[is.infinite(rank_signed_logp), rank_signed_logp := sign(beta) * 300]
results[is.infinite(rank_z),           rank_z           := sign(beta) * 10]

cat("\nRank metric distribution (rank_t = beta/se):\n"); print(summary(results$rank_t))

# SYMBOL -> ENTREZID
gene_mapping <- bitr(
  unique(results$gene_symbol),
  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
)
cat(sprintf("\nMapping SYMBOL -> ENTREZID: %d / %d (%.1f%%)\n",
            nrow(gene_mapping), length(unique(results$gene_symbol)),
            100 * nrow(gene_mapping) / length(unique(results$gene_symbol))))

results_mapped <- merge(results, gene_mapping,
                        by.x = "gene_symbol", by.y = "SYMBOL", all.x = FALSE)

# One rank per gene: where several proteins map to the same ENTREZID, keep the
# one with the largest absolute rank metric, ties broken by first occurrence.
make_gene_list <- function(dt, metric_col) {
  gene_df <- data.table(ENTREZID = dt$ENTREZID, rank = dt[[metric_col]])
  gene_df <- gene_df[order(-abs(rank))][, .SD[1L], by = ENTREZID]
  gl <- gene_df$rank
  names(gl) <- gene_df$ENTREZID
  sort(gl, decreasing = TRUE)
}

gene_list_main <- make_gene_list(results_mapped, "rank_t")
gene_list_logp <- make_gene_list(results_mapped, "rank_signed_logp")
gene_list_z    <- make_gene_list(results_mapped, "rank_z")

cat(sprintf("\nFinal gene list size : %d ENTREZIDs\n", length(gene_list_main)))

# =============================================================================
# 2. The Olink panel universe
# =============================================================================

cat("\n---- Step 2: define Olink panel universe ----\n\n")

# Universe = ENTREZIDs of every converged Olink protein that successfully mapped.
# This is the "analysis universe" — the same set fgsea implicitly draws from.
panel_entrez <- unique(as.character(results_mapped$ENTREZID))
cat(sprintf("Panel universe size       : %d ENTREZIDs\n", length(panel_entrez)))

# The panel carries ~2,920 converged proteins; symbol-to-Entrez mapping loses a
# few, so a universe far outside this window means the input or the annotation
# database is not the one the reported analysis used.
if (length(panel_entrez) < 2500 || length(panel_entrez) > 3000) {
  warning(sprintf(
    "Panel universe = %d — outside expected 2500-3000 window. Check bitr() mapping.",
    length(panel_entrez)
  ))
}

# Persist the universe for the ORA script and downstream audits.
fwrite(data.table(ENTREZID = panel_entrez),
       file.path(output_dir, "tables", "panel_universe_entrez.csv"))

# =============================================================================
# 3. Gene-set databases (full)
# =============================================================================

cat("\n---- Step 3: load MSigDB collections ----\n\n")

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark$entrez_gene, hallmark$gs_name)
cat(sprintf("  Hallmark : %d pathways\n", length(hallmark_list)))

gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
gobp_list <- split(gobp$entrez_gene, gobp$gs_name)
cat(sprintf("  GO BP    : %d pathways\n", length(gobp_list)))

kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY")
kegg_list <- split(kegg$entrez_gene, kegg$gs_name)
cat(sprintf("  KEGG     : %d pathways\n", length(kegg_list)))

reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_list <- split(reactome$entrez_gene, reactome$gs_name)
cat(sprintf("  Reactome : %d pathways\n", length(reactome_list)))

custom_keywords <- c("TASTE","OLFACT","SMELL","CHEMOSENS",
                     "NEURON","SYNAP","NERVE",
                     "INFLAMM","CYTOKINE","INTERLEUKIN",
                     "IMMUNE","COMPLEMENT","COAGUL",
                     "VIRAL","VIRUS")
gobp_custom <- gobp[grepl(paste(custom_keywords, collapse = "|"),
                          gobp$gs_name, ignore.case = TRUE), ]
gobp_custom_list <- split(gobp_custom$entrez_gene, gobp_custom$gs_name)
cat(sprintf("  Custom   : %d pathways\n", length(gobp_custom_list)))

# =============================================================================
# 4. Panel restriction
# =============================================================================

cat("\n---- Step 4: restrict each pathway to Olink panel intersection ----\n\n")

restrict_to_panel <- function(pathway_list, panel, min_size = 15, max_size = 500,
                              db_label = "") {
  before_n <- length(pathway_list)
  intersected <- lapply(pathway_list, function(g) intersect(as.character(g), panel))
  sizes_before_filter <- lengths(intersected)
  kept <- intersected[sizes_before_filter >= min_size & sizes_before_filter <= max_size]

  cat(sprintf(
    "  %-10s : %4d -> %4d after intersect+size filter (drop = %4d; median panel-size = %.0f)\n",
    db_label, before_n, length(kept), before_n - length(kept),
    median(sizes_before_filter[sizes_before_filter > 0])
  ))
  kept
}

hallmark_panel <- restrict_to_panel(hallmark_list, panel_entrez, db_label = "Hallmark")
gobp_panel     <- restrict_to_panel(gobp_list,     panel_entrez, db_label = "GO BP")
kegg_panel     <- restrict_to_panel(kegg_list,     panel_entrez, db_label = "KEGG")
reactome_panel <- restrict_to_panel(reactome_list, panel_entrez, db_label = "Reactome")
custom_panel   <- if (length(gobp_custom_list) > 0) {
  restrict_to_panel(gobp_custom_list, panel_entrez, db_label = "Custom")
} else NULL

# =============================================================================
# 5. fgsea on panel-restricted pathways
# =============================================================================

cat("\n---- Step 5: fgsea on panel-restricted pathways ----\n\n")
set.seed(42)

run_fgsea <- function(pathways, stats, name) {
  if (length(pathways) == 0) {
    cat(sprintf("  %s: no pathways after restriction, skipping.\n", name))
    return(NULL)
  }
  cat(sprintf("Running %s GSEA on %d pathways ...\n", name, length(pathways)))
  res <- fgsea(pathways = pathways, stats = stats,
               minSize = 15, maxSize = 500, nPermSimple = 10000)
  res <- res[order(pval)]
  res[, direction := ifelse(NES > 0, "Activated", "Suppressed")]
  cat(sprintf("  FDR<0.25: %d | <0.10: %d | <0.05: %d\n\n",
              sum(res$padj < 0.25, na.rm = TRUE),
              sum(res$padj < 0.10, na.rm = TRUE),
              sum(res$padj < 0.05, na.rm = TRUE)))
  res
}

gsea_hallmark <- run_fgsea(hallmark_panel, gene_list_main, "Hallmark")
gsea_gobp     <- run_fgsea(gobp_panel,     gene_list_main, "GO BP")
gsea_kegg     <- run_fgsea(kegg_panel,     gene_list_main, "KEGG")
gsea_reactome <- run_fgsea(reactome_panel, gene_list_main, "Reactome")
gsea_custom   <- if (!is.null(custom_panel)) run_fgsea(custom_panel, gene_list_main, "Custom") else NULL

# =============================================================================
# 6. Rank-metric consistency (Hallmark only, panel-restricted)
# =============================================================================

cat("\n---- Step 6: rank-metric consistency (Hallmark, panel-restricted) ----\n\n")

run_fgsea_simple <- function(pathways, stats_vec) {
  res <- fgsea(pathways = pathways, stats = stats_vec,
               minSize = 15, maxSize = 500, nPermSimple = 10000)
  res <- res[order(pval)]
  as.data.table(res[, c("pathway","size","NES","pval","padj")])
}

gsea_hallmark_t    <- run_fgsea_simple(hallmark_panel, gene_list_main)
gsea_hallmark_logp <- run_fgsea_simple(hallmark_panel, gene_list_logp)
gsea_hallmark_z    <- run_fgsea_simple(hallmark_panel, gene_list_z)

consistency <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE),
                      list(
                        gsea_hallmark_t[,    .(pathway, NES_t    = NES, FDR_t    = padj)],
                        gsea_hallmark_logp[, .(pathway, NES_logp = NES, FDR_logp = padj)],
                        gsea_hallmark_z[,    .(pathway, NES_z    = NES, FDR_z    = padj)]
                      ))
consistency[, same_direction := sign(NES_t) == sign(NES_logp) & sign(NES_t) == sign(NES_z)]

cat(sprintf("Direction agreement       : %d / %d (%.1f%%)\n",
            sum(consistency$same_direction, na.rm = TRUE),
            nrow(consistency),
            100 * mean(consistency$same_direction, na.rm = TRUE)))
cat(sprintf("NES correlation t vs logP : r=%.3f\n",
            cor(consistency$NES_t, consistency$NES_logp, use = "complete.obs")))
cat(sprintf("NES correlation t vs Z    : r=%.3f\n",
            cor(consistency$NES_t, consistency$NES_z,    use = "complete.obs")))

# =============================================================================
# 7. Format + save
# =============================================================================

cat("\n---- Step 7: format + save ----\n\n")

entrez2sym <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(results_mapped$ENTREZID), keytype = "ENTREZID", columns = "SYMBOL"
)

format_gsea_results <- function(gsea_result, db_name, entrez2sym_map) {
  if (is.null(gsea_result) || nrow(gsea_result) == 0) return(NULL)
  le_symbol <- lapply(gsea_result$leadingEdge, function(le) {
    sym <- entrez2sym_map$SYMBOL[match(as.character(le), entrez2sym_map$ENTREZID)]
    unique(sym[!is.na(sym)])
  })
  data.table(
    Database         = db_name,
    Pathway          = gsea_result$pathway,
    Size_panel       = gsea_result$size,
    NES              = round(gsea_result$NES, 3),
    P_value          = gsea_result$pval,
    FDR              = gsea_result$padj,
    Direction        = gsea_result$direction,
    Leading_Edge_N   = vapply(le_symbol, length, integer(1)),
    Leading_Edge     = vapply(le_symbol, function(x) paste(head(x, 10), collapse = ", "), character(1))
  )
}

table_dir <- file.path(output_dir, "tables")

fwrite(format_gsea_results(gsea_hallmark, "Hallmark", entrez2sym),
       file.path(table_dir, "GSEA_Hallmark_panelbg.csv"))
fwrite(format_gsea_results(gsea_gobp,     "GO_BP",    entrez2sym),
       file.path(table_dir, "GSEA_GO_BP_panelbg.csv"))
fwrite(format_gsea_results(gsea_kegg,     "KEGG",     entrez2sym),
       file.path(table_dir, "GSEA_KEGG_panelbg.csv"))
fwrite(format_gsea_results(gsea_reactome, "Reactome", entrez2sym),
       file.path(table_dir, "GSEA_Reactome_panelbg.csv"))
if (!is.null(gsea_custom)) {
  fwrite(format_gsea_results(gsea_custom, "Custom", entrez2sym),
         file.path(table_dir, "GSEA_Custom_panelbg.csv"))
}
fwrite(consistency, file.path(table_dir, "GSEA_Hallmark_rankmetric_consistency_panelbg.csv"))

all_sig <- rbindlist(list(
  if (!is.null(gsea_hallmark)) gsea_hallmark[padj < 0.25, .(Database = "Hallmark", Pathway = pathway, NES, pval, padj, direction)],
  if (!is.null(gsea_gobp))     gsea_gobp    [padj < 0.10, .(Database = "GO_BP",    Pathway = pathway, NES, pval, padj, direction)],
  if (!is.null(gsea_kegg))     gsea_kegg    [padj < 0.25, .(Database = "KEGG",     Pathway = pathway, NES, pval, padj, direction)],
  if (!is.null(gsea_reactome)) gsea_reactome[padj < 0.25, .(Database = "Reactome", Pathway = pathway, NES, pval, padj, direction)]
), fill = TRUE)
fwrite(all_sig, file.path(table_dir, "GSEA_all_significant.csv"))

gsea_results_list <- list(
  hallmark        = gsea_hallmark,
  gobp            = gsea_gobp,
  kegg            = gsea_kegg,
  reactome        = gsea_reactome,
  custom          = gsea_custom,
  gene_list       = gene_list_main,
  hallmark_panel  = hallmark_panel,
  gobp_panel      = gobp_panel,
  kegg_panel      = kegg_panel,
  reactome_panel  = reactome_panel,
  panel_entrez    = panel_entrez,
  entrez2sym      = entrez2sym,
  consistency     = consistency
)
saveRDS(gsea_results_list, file.path(output_dir, "GSEA_results_all_panelbg.rds"))

# =============================================================================
# 8. Summary
# =============================================================================

cat("\n---- SUMMARY (panel-restricted) ----\n")
for (thr in c(0.25, 0.10, 0.05)) {
  cat(sprintf("\nFDR < %.2f:\n", thr))
  cat(sprintf("  Hallmark : %d\n", if (!is.null(gsea_hallmark)) sum(gsea_hallmark$padj < thr, na.rm = TRUE) else 0))
  cat(sprintf("  GO BP    : %d\n", if (!is.null(gsea_gobp))     sum(gsea_gobp$padj     < thr, na.rm = TRUE) else 0))
  cat(sprintf("  KEGG     : %d\n", if (!is.null(gsea_kegg))     sum(gsea_kegg$padj     < thr, na.rm = TRUE) else 0))
  cat(sprintf("  Reactome : %d\n", if (!is.null(gsea_reactome)) sum(gsea_reactome$padj < thr, na.rm = TRUE) else 0))
}

if (!is.null(gsea_hallmark) && nrow(gsea_hallmark) > 0) {
  cat("\nTop Activated (Hallmark, panel-restricted):\n")
  top_act <- head(gsea_hallmark[NES > 0][order(pval)], 5)
  for (i in seq_len(nrow(top_act)))
    cat(sprintf("  %-50s NES=%+.2f  FDR=%.3g\n",
                gsub("HALLMARK_", "", top_act$pathway[i]),
                top_act$NES[i], top_act$padj[i]))

  cat("\nTop Suppressed (Hallmark, panel-restricted):\n")
  top_sup <- head(gsea_hallmark[NES < 0][order(pval)], 5)
  for (i in seq_len(nrow(top_sup)))
    cat(sprintf("  %-50s NES=%+.2f  FDR=%.3g\n",
                gsub("HALLMARK_", "", top_sup$pathway[i]),
                top_sup$NES[i], top_sup$padj[i]))
}

cat(sprintf("\nDone. Results dir: %s\n", output_dir))
cat("Next: source 02_ora_panel_universe.R, then 03_panel_pathway_coverage.R\n")
