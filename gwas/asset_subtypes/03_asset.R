# Subset-based meta-analysis (ASSET h.traits, two-sided, DLM P values) of the four per-subtype GWAS.
# Variants in all four scans with identical alleles and MAF >= 0.005 are combined; the correlation between scans
# induced by the shared controls is estimated from null variants (|z| < 3 in all four scans).
# Genome-wide significance P < 5e-8, suggestive P < 1e-6; effect allele = SAIGE Allele2.
library(data.table)
library(ASSET)
library(parallel)
library(ggplot2)
setwd(Sys.getenv("PROJECT_DIR", "."))
IN  <- "input/gwas_subtypes"   # SAIGE step-2 output subtype_{A,B,C,D}_chr{1-22}.txt
OUT <- "output/gwas_subtypes/asset"
dir.create(file.path(OUT, "merged"), recursive = TRUE, showWarnings = FALSE)
SUB <- c("A", "B", "C", "D")
N_CORES <- 16
SEG <- 20000   # variants per h.traits call

# ---- merge the four scans per chromosome ------------------------------------
for (chr in 1:22) {
  read1 <- function(s) {
    d <- fread(file.path(IN, sprintf("subtype_%s_chr%d.txt", s, chr)),
               select = c("CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "N_case", "N_ctrl"))
    setnames(d, c("Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "N_case", "N_ctrl"),
             paste0(c("a1", "a2", "af", "beta", "se", "ncase", "nctrl"), "_", s))
  }
  m <- Reduce(function(x, y) merge(x, y, by = c("CHR", "POS", "MarkerID")), lapply(SUB, read1))
  m <- m[a1_A == a1_B & a1_A == a1_C & a1_A == a1_D & a2_A == a2_B & a2_A == a2_C & a2_A == a2_D]
  m[, maf := pmin(af_A, 1 - af_A)]
  keep <- c("MarkerID", "CHR", "POS", "a1_A", "a2_A", "af_A", "maf", paste0(rep(c("beta_", "se_", "ncase_", "nctrl_"), each = 4), SUB))
  fwrite(m[maf >= 0.005, ..keep], file.path(OUT, "merged", sprintf("merged_chr%d.csv", chr)))
}

# ---- cross-scan z correlation from null variants (40,000 sampled per chromosome) ----
files <- Sys.glob(file.path(OUT, "merged", "merged_chr*.csv"))
set.seed(1)
Z <- do.call(rbind, lapply(files, function(f) {
  d <- fread(f)
  z <- as.matrix(d[, paste0("beta_", SUB), with = FALSE]) / as.matrix(d[, paste0("se_", SUB), with = FALSE])
  z[sample(nrow(z), min(40000L, nrow(z))), , drop = FALSE]
}))
cormat <- cor(Z[apply(abs(Z), 1, max) < 3, ], use = "complete.obs")
dimnames(cormat) <- list(SUB, SUB)
d1 <- fread(files[1], nrows = 1)
ncase <- as.numeric(d1[, paste0("ncase_", SUB), with = FALSE])
nctrl <- as.numeric(d1[, paste0("nctrl_", SUB), with = FALSE])

# ---- ASSET ------------------------------------------------------------------
asset_segment <- function(d) {
  bh <- as.matrix(d[, paste0("beta_", SUB), with = FALSE]); sh <- as.matrix(d[, paste0("se_", SUB), with = FALSE])
  colnames(bh) <- colnames(sh) <- SUB
  nc <- matrix(rep(ncase, each = nrow(d)), ncol = 4); nl <- matrix(rep(nctrl, each = nrow(d)), ncol = 4)
  r  <- h.traits(d$MarkerID, SUB, bh, sh, nc, nl, cor = cormat, cor.numr = FALSE, side = 2, meth.pval = "DLM")
  s2 <- as.data.table(h.summary(r)$Subset.2sided)
  for (cc in c("Pvalue", "Pvalue.1", "Pvalue.2", "OR.1", "OR.2")) set(s2, j = cc, value = as.numeric(s2[[cc]]))
  s2 <- merge(d[, .(SNP = MarkerID, chr = CHR, pos = POS, maf, a_eff = a2_A)], s2, by = "SNP")
  s2[, .(SNP, chr, pos, maf, a_eff, asset_p2 = Pvalue, p_pos = Pvalue.1, p_neg = Pvalue.2,
         OR_pos = OR.1, OR_neg = OR.2, pos_subset = Pheno.1, neg_subset = Pheno.2)]
}
res <- rbindlist(lapply(1:22, function(chr) {
  d <- fread(file.path(OUT, "merged", sprintf("merged_chr%d.csv", chr)))
  rbindlist(mclapply(seq(1, nrow(d), by = SEG), function(st) asset_segment(d[st:min(st + SEG - 1, nrow(d))]),
                     mc.cores = N_CORES))
}))
res <- res[is.finite(asset_p2) & asset_p2 > 0][order(chr, pos)]
fwrite(res, file.path(OUT, "asset_genomewide_results.csv"))
fwrite(res[asset_p2 < 1e-6][order(asset_p2)], file.path(OUT, "asset_suggestive_p1e-6.csv"))
lam <- median(qchisq(res$asset_p2, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
fwrite(data.table(n_variants = nrow(res), lambda_GC = lam, n_p_5e8 = sum(res$asset_p2 < 5e-8),
                  n_p_1e6 = sum(res$asset_p2 < 1e-6)), file.path(OUT, "asset_summary.csv"))

# ---- Manhattan and QQ -------------------------------------------------------
res[, nlp := -log10(asset_p2)]
offs <- c(0, cumsum(as.numeric(res[, max(pos), by = chr]$V1)))[1:22]
res[, cum := pos + offs[chr]]
axis <- res[, .(center = mean(cum)), by = chr]
man <- ggplot(res, aes(cum, nlp, color = factor(chr %% 2))) + geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-6), color = "orange", linetype = "dotted") +
  scale_color_manual(values = c("0" = "grey30", "1" = "grey60"), guide = "none") +
  scale_x_continuous(breaks = axis$center, labels = axis$chr) +
  labs(title = sprintf("Genome-wide ASSET subset-based subtype scan | lambda=%.3f", lam), x = "chromosome",
       y = "-log10 ASSET p (2-sided)") +
  theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank())
ggsave(file.path(OUT, "asset_manhattan.png"), man, width = 12, height = 4.5, dpi = 300)
obs <- -log10(sort(res$asset_p2)); ex <- -log10(ppoints(length(obs)))
keep <- obs > 2 | (seq_along(obs) %% 100L == 0L)   # thin the null bulk
qq <- ggplot(data.table(ex = ex[keep], obs = obs[keep]), aes(ex, obs)) + geom_point(size = 0.5, color = "grey40") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = sprintf("ASSET QQ (lambda=%.3f)", lam), x = "expected -log10 p", y = "observed -log10 p") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT, "asset_qq.png"), qq, width = 4.5, height = 4.5, dpi = 300)
