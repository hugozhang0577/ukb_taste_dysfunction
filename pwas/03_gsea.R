# Gene set enrichment analysis (fgsea) of the primary PWAS, ranked by beta/SE, with every
# MSigDB gene set restricted to the proteins measured on the Olink panel (panel universe).
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)
library(fgsea)
setwd(Sys.getenv("PROJECT_DIR", "."))

TABLE_DIR <- "output/pwas_gsea/panel_restricted/tables"
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

res <- fread("input/assoc_results/pwas_primary.csv")
res <- res[converged == TRUE & is.finite(beta) & is.finite(se) & se > 0 & is.finite(pval)]

# Ranking statistics: beta/SE (main); signed -log10(P) and signed Z for the rank-metric comparison
res[, gene_symbol      := toupper(protein)]
res[, rank_t           := beta / se]
res[, rank_signed_logp := sign(beta) * (-log10(pval))]
res[, rank_z           := qnorm(1 - pval / 2) * sign(beta)]
res[is.infinite(rank_signed_logp), rank_signed_logp := sign(beta) * 300]
res[is.infinite(rank_z),           rank_z           := sign(beta) * 10]

gene_map <- bitr(unique(res$gene_symbol), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
res <- merge(res, gene_map, by.x = "gene_symbol", by.y = "SYMBOL")

# One value per Entrez ID (largest |rank|), sorted decreasing
make_gene_list <- function(dt, metric) {
  g <- data.table(ENTREZID = dt$ENTREZID, rank = dt[[metric]])
  g <- g[order(-abs(rank))][, .SD[1L], by = ENTREZID]
  sort(setNames(g$rank, g$ENTREZID), decreasing = TRUE)
}
gene_list_t    <- make_gene_list(res, "rank_t")
gene_list_logp <- make_gene_list(res, "rank_signed_logp")
gene_list_z    <- make_gene_list(res, "rank_z")

# Panel universe: every mapped protein in the scan
panel_entrez <- unique(as.character(res$ENTREZID))

# MSigDB collections; custom = GO BP sets matching chemosensory, neural, immune or viral keywords
sets_of <- function(m) split(m$entrez_gene, m$gs_name)
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
CUSTOM_KEYWORDS <- c("TASTE", "OLFACT", "SMELL", "CHEMOSENS", "NEURON", "SYNAP", "NERVE",
                     "INFLAMM", "CYTOKINE", "INTERLEUKIN", "IMMUNE", "COMPLEMENT", "COAGUL",
                     "VIRAL", "VIRUS")
collections <- list(
  Hallmark = sets_of(msigdbr(species = "Homo sapiens", category = "H")),
  GO_BP    = sets_of(gobp),
  KEGG     = sets_of(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY")),
  Reactome = sets_of(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")),
  Custom   = sets_of(gobp[grepl(paste(CUSTOM_KEYWORDS, collapse = "|"), gobp$gs_name, ignore.case = TRUE), ])
)

# Restrict each gene set to the panel universe; keep sets of 15-500 panel genes
restrict_to_panel <- function(sets) {
  s <- lapply(sets, function(g) intersect(as.character(g), panel_entrez))
  s[lengths(s) >= 15 & lengths(s) <= 500]
}
panel_sets <- lapply(collections, restrict_to_panel)

set.seed(42)
run_fgsea <- function(pathways, stats) {
  r <- fgsea(pathways = pathways, stats = stats, minSize = 15, maxSize = 500, nPermSimple = 10000)
  r <- r[order(pval)]
  r[, direction := ifelse(NES > 0, "Activated", "Suppressed")]
}
gsea <- lapply(panel_sets, run_fgsea, stats = gene_list_t)

# Rank-metric comparison (Hallmark): beta/SE vs signed -log10(P) vs signed Z
nes_by_metric <- lapply(list(t = gene_list_t, logp = gene_list_logp, z = gene_list_z), function(s)
  run_fgsea(panel_sets$Hallmark, s)[, .(pathway, NES, padj)])
consistency <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE), list(
  nes_by_metric$t[,    .(pathway, NES_t    = NES, FDR_t    = padj)],
  nes_by_metric$logp[, .(pathway, NES_logp = NES, FDR_logp = padj)],
  nes_by_metric$z[,    .(pathway, NES_z    = NES, FDR_z    = padj)]))
consistency[, same_direction := sign(NES_t) == sign(NES_logp) & sign(NES_t) == sign(NES_z)]

# Result tables; leading edge reported as gene symbols (first 10)
entrez2sym <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(res$ENTREZID),
                                    keytype = "ENTREZID", columns = "SYMBOL")
format_gsea <- function(r, db) {
  le <- lapply(r$leadingEdge, function(x) {
    s <- entrez2sym$SYMBOL[match(as.character(x), entrez2sym$ENTREZID)]
    unique(s[!is.na(s)])
  })
  data.table(Database = db, Pathway = r$pathway, Size_panel = r$size, NES = round(r$NES, 3),
             P_value = r$pval, FDR = r$padj, Direction = r$direction,
             Leading_Edge_N = vapply(le, length, integer(1)),
             Leading_Edge   = vapply(le, function(x) paste(head(x, 10), collapse = ", "), character(1)))
}
for (db in names(gsea))
  fwrite(format_gsea(gsea[[db]], db), file.path(TABLE_DIR, sprintf("GSEA_%s_panelbg.csv", db)))
fwrite(consistency, file.path(TABLE_DIR, "GSEA_Hallmark_rankmetric_consistency_panelbg.csv"))

# Significant sets: FDR < 0.25 (Hallmark, KEGG, Reactome), FDR < 0.10 (GO BP)
FDR_MAX <- c(Hallmark = 0.25, GO_BP = 0.10, KEGG = 0.25, Reactome = 0.25)
all_sig <- rbindlist(lapply(names(FDR_MAX), function(db)
  gsea[[db]][padj < FDR_MAX[[db]], .(Database = db, Pathway = pathway, NES, pval, padj, direction)]),
  fill = TRUE)
fwrite(all_sig, file.path(TABLE_DIR, "GSEA_all_significant.csv"))
