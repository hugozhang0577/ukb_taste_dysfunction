# Summary tables from the FUMA SNP2GENE job (P < 1e-5 lead variants): lead variants with annotation,
# genomic risk loci, mapped genes, GWAS Catalog overlap by trait category, and functional consequences.
library(data.table)
setwd(Sys.getenv("PROJECT_DIR", "."))
FUMA_DIR <- "input/fuma/FUMA_1e5"   # job directory downloaded from FUMA
OUT_DIR  <- "output/gwas_fuma/summary"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) fread(file.path(FUMA_DIR, f))
leadSNPs        <- rd("leadSNPs.txt")
GenomicRiskLoci <- rd("GenomicRiskLoci.txt")
snps            <- rd("snps.txt")
annov           <- rd("annov.txt")
genes           <- rd("genes.txt")
gwascatalog     <- rd("gwascatalog.txt")

snp_info <- unique(snps[, .(uniqID, rsID, non_effect_allele, effect_allele, MAF, beta, se,
                            nearestGene, dist, func, CADD, RDB)], by = "uniqID")
annov_summary <- annov[, .(annot_detail = paste(unique(annot), collapse = ";"),
                           annotated_genes = paste(unique(symbol[symbol != ""]), collapse = ";")), by = uniqID]
lead <- merge(leadSNPs, snp_info, by = "uniqID", all.x = TRUE, suffixes = c("", ".snp"))
lead <- merge(lead, annov_summary, by = "uniqID", all.x = TRUE)
lead <- merge(lead, GenomicRiskLoci[, .(GenomicLocus, start, end, nSNPs, nGWASSNPs)], by = "GenomicLocus", all.x = TRUE)
lead[, rsID := fifelse(is.na(rsID), rsID.snp, rsID)]
setorder(lead, p)
lead[, significance := fcase(p < 5e-8, "Genome-wide significant", p < 1e-5, "Suggestive", default = "Not significant")]
lead_gws <- lead[significance == "Genome-wide significant"]

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
func_summary <- snps[, .(n = .N), by = func][order(-n)][, pct := round(n / sum(n) * 100, 1)]

fwrite(lead_gws[, .(Locus = GenomicLocus, rsID, Chr = chr, Position = pos,
                    `EA/NEA` = paste0(effect_allele, "/", non_effect_allele),
                    MAF = sprintf("%.3f", MAF), Beta = sprintf("%.3f", beta), SE = sprintf("%.4f", se),
                    P = formatC(p, format = "e", digits = 2), `Nearest Gene` = nearestGene,
                    Function = func, CADD = sprintf("%.1f", CADD))],
       file.path(OUT_DIR, "Table1_GWS_loci.csv"))
fwrite(lead[, .(Locus = GenomicLocus, rsID, Chr = chr, Position = pos, `Effect Allele` = effect_allele,
                `Non-effect Allele` = non_effect_allele, MAF, Beta = beta, SE = se, P = p,
                `Nearest Gene` = nearestGene, `Distance to Gene` = dist, Function = func,
                `Detailed Annotation` = annot_detail, CADD, RDB, `nSNPs in Locus` = nSNPs, Significance = significance)],
       file.path(OUT_DIR, "TableS1_All_LeadSNPs.csv"))
fwrite(merge(GenomicRiskLoci, unique(lead[, .(GenomicLocus, nearestGene)]), by = "GenomicLocus", all.x = TRUE)[
         , .(GenomicLocus, chr, start, end, rsID, p, nearestGene, nSNPs, nGWASSNPs, nIndSigSNPs, nLeadSNPs)],
       file.path(OUT_DIR, "TableS2_GenomicRiskLoci.csv"))
fwrite(genes[order(minGwasP), .(symbol, ensg, chr, start, end, type, posMapSNPs, eqtlMapSNPs,
                                eqtlMapminP, eqtlMapts, minGwasP, pLI, GenomicLocus)],
       file.path(OUT_DIR, "TableS3_MappedGenes.csv"))
fwrite(gwascatalog[order(GenomicLocus, P), .(GenomicLocus, IndSigSNP, snp, Trait, category, PMID, FirstAuth,
                                             Journal, P, OrBeta, ReportedGene, MappedGene)],
       file.path(OUT_DIR, "TableS4_GWASCatalog.csv"))
fwrite(func_summary, file.path(OUT_DIR, "TableS5_FunctionalAnnotation.csv"))
fwrite(lead[, .(rsID, chr, pos, p, effect_allele, non_effect_allele, MAF, beta, se, nearestGene, func,
                significance, GenomicLocus)],
       file.path(OUT_DIR, "LeadSNPs.csv"))
