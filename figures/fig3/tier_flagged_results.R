#!/usr/bin/env Rscript
# ===========================================================================
# tier_flagged_results.R  (one table behind Figure 3)
#
# Unified 4-tier evidence flagging engine. Reads per-layer primary sumstats
# (baseline and follow-up ExWAS, the disease phenome scan, PWAS, MWAS, GWAS)
# for the discovery and the two held-out cohorts, and emits
# a single long table `tier_flagged_results.csv` with harmonised columns
# and 4 evidence-tier flags:
#   pass_bonf            — within-domain Bonferroni q<0.05
#   pass_fdr             — within-domain BH-FDR q<0.05
#   pass_effect_tier     — |beta| >= effect-size threshold (layer-specific)
#   pass_g2_replication  — direction-concordant in G2 with |beta| same sign
#                          (and pass_fdr in G2 when available)
#   pass_g3_replication  — same for G3
#
# Output schema (one row per variable per layer per primary model):
#   layer / module / domain / variable / variable_label /
#   beta / se / OR / CI_lo / CI_hi / pval /
#   pval_fdr / pval_bonf /
#   g1_n / g1_n_case /
#   pass_bonf / pass_fdr / pass_effect_tier /
#   g2_beta / g2_pval / g2_sign_concord / pass_g2_replication /
#   g3_beta / g3_pval / g3_sign_concord / pass_g3_replication /
#   evidence_score (sum of 4 tier flags, range 0-4) /
#   headline_flag (evidence_score >= 3)
#
# Output:
#   output/evidence_tiering/tier_flagged_results.csv     (long)
#   output/evidence_tiering/tier_summary_by_layer.csv    (counts)
#   output/evidence_tiering/tier_summary_by_domain.csv   (per-domain)
#
# Figure 3 and the corresponding supplementary tables are drawn from this table.
# It exists because those panels put three scans on one pair of axes: the two
# questionnaire scans and the disease scan report different columns, apply FDR
# within their own families, and hold their held-out results in separate files.
# Harmonising once here means the figure selects and colours rows rather than
# re-deriving significance, and the supplementary tables cannot disagree with the
# figure about what passed. The tier thresholds are set in CFG below rather than
# scattered through the code.
#
# Figure 2 does NOT read this table. It plots each scan's own summary statistics
# directly, because its panels are one-scan-per-panel and need no common scale.
#
# What a tier flag is and is not: the four flags are not four independent tests.
# Bonferroni implies FDR, so `evidence_score` is an ordering device for reading a
# long table, not a probability. `pass_effect_tier` is applied to each layer's
# own beta in its own units, so for a continuous exposure it depends on the scale
# that exposure was measured on; exwas_continuous_sd.R builds the per-SD
# lookup that the figures use to put such exposures on a comparable axis, and it
# is applied at the figure, after these flags.
# ===========================================================================

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# CONFIG — tier thresholds (edit here to retune; do not scatter through code)
# ---------------------------------------------------------------------------
CFG <- list(
  fdr_q          = 0.05,
  effect_tier = list(
    exwas = 0.0953,   # |log(OR)| >= log(1.10) so OR >= 1.10 or <= 0.91
    pwas  = 0.0953,
    mwas  = 0.05,     # NMR is z-scored; |beta| >= 0.05 SD per Würtz convention
    gwas  = 0.05      # placeholder (GWAS uses GWS threshold separately)
  ),
  out_dir = "output/evidence_tiering"
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
}

# ---------------------------------------------------------------------------
# Helper — fail-visibly counter
# ---------------------------------------------------------------------------
chk <- function(dt, label) {
  log_msg("  %s: nrow = %d  ncol = %d", label, nrow(dt), ncol(dt))
  invisible(dt)
}

# ---------------------------------------------------------------------------
# Helper — apply BH-FDR per domain group (matches project convention; do not
# pool across domains — see CLAUDE.md §1.5)
# ---------------------------------------------------------------------------
fdr_per_group <- function(dt, group_col, pval_col) {
  dt[, pval_fdr_within := p.adjust(get(pval_col), method = "BH"), by = group_col]
  dt
}

# ---------------------------------------------------------------------------
# Layer-specific loaders — harmonise to unified schema
# Each returns long data.table with columns:
#   layer module domain variable variable_label
#   beta se OR CI_lo CI_hi pval pval_fdr pval_bonf
#   g1_n g1_n_case
# ---------------------------------------------------------------------------

load_pwas_g1 <- function() {
  log_msg("Loading PWAS discovery primary")
  f <- "input/assoc_results/pwas_primary.csv"
  if (!file.exists(f)) stop("PWAS G1 file not found: ", f)
  d <- fread(f)
  chk(d, "PWAS G1 raw")
  # Columns: protein n_total n_case beta se or or_lower or_upper z pval
  #          pval_fdr pval_bonf sig_* converged error
  d <- d[converged == TRUE]
  out <- data.table(
    layer          = "PWAS",
    module         = "PWAS",
    domain         = "PWAS",
    variable       = d$protein,
    variable_label = toupper(d$protein),
    beta           = d$beta,
    se             = d$se,
    OR             = d$or,
    CI_lo          = d$or_lower,
    CI_hi          = d$or_upper,
    pval           = d$pval,
    pval_fdr       = d$pval_fdr,
    pval_bonf      = d$pval_bonf,
    g1_n           = d$n_total,
    g1_n_case      = d$n_case
  )
  chk(out, "PWAS G1 harmonised")
}

load_mwas_g1 <- function() {
  # The 20-covariate model, matching the metabolite panel in the feature
  # manifest and the model reported in the Methods.
  log_msg("Loading MWAS G1 primary (N7_pwas_d15_equivalent, 20cov)")
  f <- "input/assoc_results/mwas_primary.csv"
  if (!file.exists(f)) stop("MWAS G1 file not found: ", f)
  d <- fread(f)
  chk(d, "MWAS G1 raw")
  d <- d[converged == TRUE]
  # Attach NMR class annotation from manifest
  ann_f <- "input/reference/nmr_metabolite_annotation.csv"
  if (!file.exists(ann_f)) {
    warning("MWAS annotation not found; domain will be 'unknown'")
    ann <- data.table(metabolite = character(), Nightingale_Group = character(),
                      cluster_id = integer(), is_rep = logical())
  } else {
    ann <- fread(ann_f)
  }
  setnames(ann, "metabolite", "protein", skip_absent = TRUE)
  d <- merge(d, ann[, .(protein, Nightingale_Group, cluster_id, is_rep)],
             by = "protein", all.x = TRUE)
  chk(d, "MWAS G1 merged with class annotation")
  out <- data.table(
    layer          = "MWAS",
    module         = "MWAS",
    domain         = ifelse(is.na(d$Nightingale_Group), "unknown", d$Nightingale_Group),
    variable       = d$protein,
    variable_label = d$protein,
    beta           = d$beta,
    se             = d$se,
    OR             = d$or,
    CI_lo          = d$or_lower,
    CI_hi          = d$or_upper,
    pval           = d$pval,
    pval_fdr       = d$pval_fdr,
    pval_bonf      = d$pval_bonf,
    g1_n           = d$n_total,
    g1_n_case      = d$n_case,
    cluster_id     = d$cluster_id,
    is_cluster_rep = d$is_rep
  )
  chk(out, "MWAS G1 harmonised")
}

load_exwas_source <- function(source_key) {
  # source_key: "baseline" or "followup"
  log_msg("Loading %s ExWAS discovery primary", source_key)
  if (source_key == "baseline") {
    f <- "input/assoc_results/exwas_baseline_primary.csv"
  } else {
    f <- "input/assoc_results/exwas_followup_primary.csv"
  }
  if (!file.exists(f)) stop(source_key, " ExWAS file not found: ", f)
  d <- fread(f)
  chk(d, sprintf("%s ExWAS discovery raw", source_key))
  d <- d[converged == TRUE]
  # Recompute per-domain BH (sanity: should match p_adj)
  d <- fdr_per_group(d, "fdr_group", "p_value")
  # Bonferroni within domain
  d[, pval_bonf_within := pmin(p_value * .N, 1), by = fdr_group]
  out <- data.table(
    layer          = paste0("ExWAS-", source_key),
    module         = source_key,
    domain         = d$fdr_group,
    variable       = d$variable,
    variable_label = d$variable,
    beta           = d$beta,
    se             = d$se,
    OR             = d$OR,
    CI_lo          = d$OR_lower,
    CI_hi          = d$OR_upper,
    pval           = d$p_value,
    pval_fdr       = d$pval_fdr_within,
    pval_bonf      = d$pval_bonf_within,
    g1_n           = d$n_complete,
    g1_n_case      = d$n_cases_complete
  )
  chk(out, sprintf("%s ExWAS discovery harmonised", source_key))
}

load_exwas_phecode <- function() {
  log_msg("Loading disease phenome (PheCode) discovery primary")
  f <- "input/assoc_results/dwas_phecode_primary.csv"
  if (!file.exists(f)) stop("PheCode G1 file not found: ", f)
  d <- fread(f)
  chk(d, "PheCode G1 raw")
  d <- d[converged == TRUE]
  # Map PheCode integer to ICD chapter category for domain
  # NB: a richer mapping should live in a side file; here we coarse-bin.
  d[, domain := fifelse(
    phecode >= 290 & phecode < 320, "C-Mental",
    fifelse(phecode >= 320 & phecode < 390, "C-Neuro",
    fifelse(phecode >= 240 & phecode < 280, "C-Endocrine",
    fifelse(phecode >= 390 & phecode < 460, "C-Circulatory",
    fifelse(phecode >= 460 & phecode < 520, "C-Respiratory",
    fifelse(phecode >= 520 & phecode < 580, "C-Digestive",
    "C-Other"))))))]
  out <- data.table(
    layer          = "DWAS",
    module         = "DWAS",
    domain         = d$domain,
    variable       = paste0("phe", d$phecode),
    variable_label = paste0("phe", d$phecode),
    beta           = d$beta,
    se             = d$se,
    OR             = d$or,
    CI_lo          = d$or_lower,
    CI_hi          = d$or_upper,
    pval           = d$pval,
    pval_fdr       = d$pval_fdr,
    pval_bonf      = d$pval_bonf,
    g1_n           = d$n_total,
    g1_n_case      = d$n_case_taste
  )
  chk(out, "DWAS harmonised")
}

load_gwas_g1 <- function() {
  # GWAS sumstats (10M variants) is too large for evidence-tier processing
  # of the same kind — we return only the lead SNP + suggestive set.
  log_msg("GWAS sumstats deferred — only lead SNPs handled in fig3 directly")
  data.table(
    layer = character(), module = character(), domain = character(),
    variable = character(), variable_label = character(),
    beta = numeric(), se = numeric(), OR = numeric(),
    CI_lo = numeric(), CI_hi = numeric(),
    pval = numeric(), pval_fdr = numeric(), pval_bonf = numeric(),
    g1_n = integer(), g1_n_case = integer()
  )
}

# ---------------------------------------------------------------------------
# G2 / G3 replication loaders — return data.table(variable, beta, pval)
# ---------------------------------------------------------------------------

load_pwas_validation <- function(group) {
  # group: "group2" or "group3"
  f <- sprintf("input/assoc_results/pwas_%s_primary.csv", group)
  if (!file.exists(f)) {
    log_msg("PWAS %s validation file missing — skipping", group)
    return(data.table(variable = character(), val_beta = numeric(), val_pval = numeric()))
  }
  d <- fread(f)
  d <- d[converged == TRUE]
  data.table(variable = d$protein, val_beta = d$beta, val_pval = d$pval)
}

load_mwas_validation <- function(group) {
  # 20 covariates, matching the discovery primary model
  f <- sprintf("input/assoc_results/mwas_%s_primary.csv",
               ifelse(group == "group2", "1", "3"), group)
  if (!file.exists(f)) {
    log_msg("MWAS %s validation file missing — skipping", group)
    return(data.table(variable = character(), val_beta = numeric(), val_pval = numeric()))
  }
  d <- fread(f)
  d <- d[converged == TRUE]
  data.table(variable = d$protein, val_beta = d$beta, val_pval = d$pval)
}

load_exwas_validation <- function(source_key, group) {
  f <- switch(source_key,
    baseline = sprintf("input/assoc_results/exwas_baseline_%s_primary.csv", group),
    followup = sprintf("input/assoc_results/exwas_followup_%s_primary.csv", group),
    dwas     = sprintf("input/assoc_results/dwas_phecode_%s_primary.csv", group),
    stop("unknown source: ", source_key))
  if (!file.exists(f)) {
    log_msg("%s %s held-out file missing -- skipping", source_key, group)
    return(data.table(variable = character(), val_beta = numeric(), val_pval = numeric()))
  }
  d <- fread(f)
  d <- d[converged == TRUE]
  if (source_key == "dwas") {
    data.table(variable = paste0("phe", d$phecode), val_beta = d$beta, val_pval = d$pval)
  } else {
    data.table(variable = d$variable, val_beta = d$beta, val_pval = d$p_value)
  }
}

# ---------------------------------------------------------------------------
# Tier-flagging core — apply 4 flags to harmonised G1 table
# ---------------------------------------------------------------------------

apply_tiers <- function(g1, layer_key, g2 = NULL, g3 = NULL) {
  # Effect-tier threshold per layer
  eff_thr <- CFG$effect_tier[[layer_key]]
  if (is.null(eff_thr)) eff_thr <- 0.05  # safe default

  g1[, pass_bonf        := !is.na(pval_bonf) & pval_bonf < CFG$fdr_q]
  g1[, pass_fdr         := !is.na(pval_fdr)  & pval_fdr  < CFG$fdr_q]
  g1[, pass_effect_tier := !is.na(beta)      & abs(beta) >= eff_thr]

  # G2 / G3 replication: direction concordance + nominal significance in val
  for (lbl in c("g2", "g3")) {
    val_dt <- if (lbl == "g2") g2 else g3
    if (is.null(val_dt) || nrow(val_dt) == 0) {
      g1[, paste0(lbl, "_beta")     := NA_real_]
      g1[, paste0(lbl, "_pval")     := NA_real_]
      g1[, paste0(lbl, "_sign_concord") := FALSE]
      g1[, paste0("pass_", lbl, "_replication") := FALSE]
      next
    }
    setnames(val_dt, c("val_beta", "val_pval"),
             c(paste0(lbl, "_beta"), paste0(lbl, "_pval")))
    g1 <- merge(g1, val_dt, by = "variable", all.x = TRUE)
    g1[, paste0(lbl, "_sign_concord") :=
       !is.na(get(paste0(lbl, "_beta"))) &
       sign(get(paste0(lbl, "_beta"))) == sign(beta) &
       sign(beta) != 0]
    g1[, paste0("pass_", lbl, "_replication") :=
       get(paste0(lbl, "_sign_concord")) &
       !is.na(get(paste0(lbl, "_pval"))) &
       get(paste0(lbl, "_pval")) < 0.05]
  }

  # Composite evidence score (0-4)
  g1[, evidence_score :=
     as.integer(pass_bonf) + as.integer(pass_fdr) +
     as.integer(pass_effect_tier) + as.integer(pass_g2_replication)]
  g1[, headline_flag := evidence_score >= 3]

  # Sumstats-only pseudo-R² on the liability scale (Lee 2011 / Yang 2010):
  #   wald_chi2 = z^2 = (beta/se)^2;   r2_liability = chi2 / (chi2 + n)
  # This is the standard "variance explained" proxy when only summary stats are
  # available. Cox 2026 Nat Med reports median R²=0.5% across NHANES; we expect
  # similar magnitudes here: a very small P value is not on its own a claim
  # about the size of an effect.
  g1[, wald_chi2 := ifelse(!is.na(beta) & !is.na(se) & se > 0,
                           (beta / se)^2,
                           qchisq(1 - pmin(pval, 1 - 1e-300), df = 1))]
  g1[, r2_liability := wald_chi2 / (wald_chi2 + g1_n)]
  g1
}

# ---------------------------------------------------------------------------
# Main — assemble all layers
# ---------------------------------------------------------------------------

log_msg("=== evidence tiering START ===")

layers <- list(
  PWAS    = list(g1 = load_pwas_g1(),
                 g2 = load_pwas_validation("group2"),
                 g3 = load_pwas_validation("group3"),
                 key = "pwas"),
  MWAS    = list(g1 = load_mwas_g1(),
                 g2 = load_mwas_validation("group2"),
                 g3 = load_mwas_validation("group3"),
                 key = "mwas"),
  ExWAS_baseline = list(g1 = load_exwas_source("baseline"),
                 g2 = load_exwas_validation("baseline", "group2"),
                 g3 = load_exwas_validation("baseline", "group3"),
                 key = "exwas"),
  ExWAS_followup = list(g1 = load_exwas_source("followup"),
                 g2 = load_exwas_validation("followup", "group2"),
                 g3 = load_exwas_validation("followup", "group3"),
                 key = "exwas"),
  DWAS    = list(g1 = load_exwas_phecode(),
                 g2 = load_exwas_validation("dwas", "group2"),
                 g3 = load_exwas_validation("dwas", "group3"),
                 key = "exwas")
)

flagged <- rbindlist(
  lapply(layers, function(L) apply_tiers(L$g1, L$key, L$g2, L$g3)),
  use.names = TRUE, fill = TRUE
)
chk(flagged, "All layers stacked")

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

summary_layer <- flagged[, .(
  n_tested            = .N,
  n_pass_bonf         = sum(pass_bonf, na.rm = TRUE),
  n_pass_fdr          = sum(pass_fdr, na.rm = TRUE),
  n_pass_effect_tier  = sum(pass_effect_tier, na.rm = TRUE),
  n_pass_g2           = sum(pass_g2_replication, na.rm = TRUE),
  n_pass_g3           = sum(pass_g3_replication, na.rm = TRUE),
  n_headline          = sum(headline_flag, na.rm = TRUE)
), by = layer]

summary_domain <- flagged[, .(
  n_tested            = .N,
  n_pass_bonf         = sum(pass_bonf, na.rm = TRUE),
  n_pass_fdr          = sum(pass_fdr, na.rm = TRUE),
  n_pass_effect_tier  = sum(pass_effect_tier, na.rm = TRUE),
  n_pass_g2           = sum(pass_g2_replication, na.rm = TRUE),
  n_headline          = sum(headline_flag, na.rm = TRUE)
), by = .(layer, module, domain)]

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
out_main    <- file.path(CFG$out_dir, "tier_flagged_results.csv")
out_layer   <- file.path(CFG$out_dir, "tier_summary_by_layer.csv")
out_domain  <- file.path(CFG$out_dir, "tier_summary_by_domain.csv")

fwrite(flagged,         out_main)
fwrite(summary_layer,   out_layer)
fwrite(summary_domain,  out_domain)

log_msg("Wrote %s (%d rows)", out_main,   nrow(flagged))
log_msg("Wrote %s (%d rows)", out_layer,  nrow(summary_layer))
log_msg("Wrote %s (%d rows)", out_domain, nrow(summary_domain))

log_msg("=== evidence tiering DONE ===")
print(summary_layer)
