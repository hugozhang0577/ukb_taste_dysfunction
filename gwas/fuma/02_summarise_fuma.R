#!/usr/bin/env Rscript
# =============================================================================
# FUMA — summarise SNP2GENE output into manuscript tables
# =============================================================================
# Reads the FUMA SNP2GENE output directory (P < 1e-5 run), integrates lead SNPs
# with per-SNP annotation, genomic risk loci, mapped genes and GWAS-Catalog
# overlap, and writes the main-text + supplementary tables.
#
# Outputs (under <FUMA_DIR>/FUMA_summary/):
#   Table1_GWS_loci.csv            genome-wide significant lead loci (main text)
#   TableS1_All_LeadSNPs.csv       all lead SNPs (P < 1e-5)
#   TableS2_GenomicRiskLoci.csv    genomic risk loci
#   TableS3_MappedGenes.csv        positional + eQTL mapped genes
#   TableS4_GWASCatalog.csv        GWAS-Catalog overlap (categorised)
#   TableS5_FunctionalAnnotation.csv  candidate-SNP functional categories
#   LeadSNPs.csv                   lead-SNP master list (chr/pos/effect/gene)
# =============================================================================
suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
FUMA_DIR <- Sys.getenv("FUMA_DIR",
  file.path(PROJECT_DIR, "gwas/cohort_primary/SAIGE/fuma"))
fuma_1e5   <- file.path(FUMA_DIR, "FUMA_1e5")
output_dir <- file.path(FUMA_DIR, "FUMA_summary")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) fread(file.path(fuma_1e5, f))
leadSNPs        <- rd("leadSNPs.txt")
IndSigSNPs      <- rd("IndSigSNPs.txt")
GenomicRiskLoci <- rd("GenomicRiskLoci.txt")
snps            <- rd("snps.txt")
annov           <- rd("annov.txt")
genes           <- rd("genes.txt")
gwascatalog     <- rd("gwascatalog.txt")
cat(sprintf("lead=%d indSig=%d loci=%d candSNP=%d genes=%d gwasCat=%d\n",
            nrow(leadSNPs), nrow(IndSigSNPs), nrow(GenomicRiskLoci),
            nrow(snps), nrow(genes), nrow(gwascatalog)))

# ---- integrate lead SNPs with annotation -----------------------------------
snp_info <- unique(snps[, .(uniqID, rsID, non_effect_allele, effect_allele, MAF,
                            beta, se, nearestGene, dist, func, CADD, RDB)], by = "uniqID")
annov_summary <- annov[, .(
  annot_detail   = paste(unique(annot), collapse = ";"),
  annotated_genes = paste(unique(symbol[symbol != ""]), collapse = ";")), by = uniqID]

lead <- merge(leadSNPs, snp_info, by = "uniqID", all.x = TRUE, suffixes = c("", ".snp"))
lead <- merge(lead, annov_summary, by = "uniqID", all.x = TRUE)
lead <- merge(lead, GenomicRiskLoci[, .(GenomicLocus, start, end, nSNPs, nGWASSNPs)],
              by = "GenomicLocus", all.x = TRUE)
# lead SNP rsID lives in rsID.snp after the snp_info merge; fall back if absent
if ("rsID.snp" %in% names(lead)) lead[, rsID := fifelse(is.na(rsID), rsID.snp, rsID)]
setorder(lead, p)

lead[, significance := fcase(p < 5e-8, "Genome-wide significant",
                             p < 1e-5, "Suggestive", default = "Not significant")]
lead_gws <- lead[significance == "Genome-wide significant"]
cat(sprintf("GWS lead loci: %d | suggestive: %d\n",
            nrow(lead_gws), sum(lead$significance == "Suggestive")))

# ---- GWAS-Catalog overlap (categorised) ------------------------------------
gwascatalog[, category := fcase(
  grepl("alzheimer|dementia|cognitive|memory|brain", Trait, ignore.case = TRUE), "Neurological",
  grepl("cholesterol|lipid|ldl|hdl|triglyceride|apolipoprotein", Trait, ignore.case = TRUE), "Lipids",
  grepl("heart|coronary|myocardial|cardiovascular|stroke|atrial", Trait, ignore.case = TRUE), "Cardiovascular",
  grepl("diabetes|glucose|insulin|hba1c|glycemic", Trait, ignore.case = TRUE), "Metabolic/Diabetes",
  grepl("bmi|obesity|body mass|weight|waist|fat", Trait, ignore.case = TRUE), "Anthropometric",
  grepl("blood pressure|hypertension", Trait, ignore.case = TRUE), "Blood Pressure",
  grepl("longevity|lifespan|aging|age", Trait, ignore.case = TRUE), "Longevity/Aging",
  grepl("inflammation|crp|interleukin|immune", Trait, ignore.case = TRUE), "Inflammatory/Immune",
  default = "Other")]

# ---- gene prioritisation + functional summary ------------------------------
# Mapped genes are reported ordered by their strongest GWAS P value, and both
# mapping counts are carried so a reader can see how each gene was reached. No
# composite score is computed: any weighting of positional against eQTL evidence
# would be an arbitrary choice, and the Methods do not specify one.
func_summary <- snps[, .(n = .N), by = func][order(-n)][, pct := round(n / sum(n) * 100, 1)]

# ---- write tables ----------------------------------------------------------
# Export and dx upload to RAP  (Table1 + S1-S5 + LeadSNPs.csv)
fwrite(lead_gws[, .(Locus = GenomicLocus, rsID, Chr = chr, Position = pos,
                    `EA/NEA` = paste0(effect_allele, "/", non_effect_allele),
                    MAF = sprintf("%.3f", MAF), Beta = sprintf("%.3f", beta),
                    SE = sprintf("%.4f", se), P = formatC(p, format = "e", digits = 2),
                    `Nearest Gene` = nearestGene, Function = func,
                    CADD = sprintf("%.1f", CADD))],
       file.path(output_dir, "Table1_GWS_loci.csv"))

fwrite(lead[, .(Locus = GenomicLocus, rsID, Chr = chr, Position = pos,
                `Effect Allele` = effect_allele, `Non-effect Allele` = non_effect_allele,
                MAF, Beta = beta, SE = se, P = p, `Nearest Gene` = nearestGene,
                `Distance to Gene` = dist, Function = func,
                `Detailed Annotation` = annot_detail, CADD, RDB,
                `nSNPs in Locus` = nSNPs, Significance = significance)],
       file.path(output_dir, "TableS1_All_LeadSNPs.csv"))

fwrite(merge(GenomicRiskLoci, unique(lead[, .(GenomicLocus, nearestGene)]),
             by = "GenomicLocus", all.x = TRUE)[
             , .(GenomicLocus, chr, start, end, rsID, p, nearestGene,
                 nSNPs, nGWASSNPs, nIndSigSNPs, nLeadSNPs)],
       file.path(output_dir, "TableS2_GenomicRiskLoci.csv"))

fwrite(genes[order(minGwasP),
       .(symbol, ensg, chr, start, end, type, posMapSNPs, eqtlMapSNPs,
         eqtlMapminP, eqtlMapts, minGwasP, pLI, GenomicLocus)],
       file.path(output_dir, "TableS3_MappedGenes.csv"))

fwrite(gwascatalog[order(GenomicLocus, P),
       .(GenomicLocus, IndSigSNP, snp, Trait, category, PMID, FirstAuth,
         Journal, P, OrBeta, ReportedGene, MappedGene)],
       file.path(output_dir, "TableS4_GWASCatalog.csv"))

fwrite(func_summary, file.path(output_dir, "TableS5_FunctionalAnnotation.csv"))

fwrite(lead[, .(rsID, chr, pos, p, effect_allele, non_effect_allele, MAF, beta, se,
                nearestGene, func, significance, GenomicLocus)],
       file.path(output_dir, "LeadSNPs.csv"))

cat("FUMA summary tables ->", output_dir, "\n")
