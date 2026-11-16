#!/usr/bin/env Rscript
# Follow-up questionnaire exposures, block 4: experience of pain (Category 154): chronic and
# widespread pain, pain conditions, BPI severity and interference, DN4 neuropathic items, pain sites,
# fibromyalgia symptom-severity items and self-rated health. Runs in the R session of 01.
library(data.table)

yn <- c("Yes" = 1L, "yes" = 1L, "1" = 1L, "No" = 0L, "no" = 0L, "0" = 0L)
recode_yn <- function(x) unname(yn[as.character(x)])
in_range <- function(x, hi) { x <- suppressWarnings(as.numeric(x)); x[which(x < 0 | x > hi)] <- NA; x }
# a site counts as painful unless answered "No"
site_present_flag <- function(x) {
  if (is.numeric(x)) return(fifelse(!is.na(x) & x > 0, 1L, fifelse(is.na(x), NA_integer_, 0L)))
  xc <- as.character(x)
  fifelse(is.na(xc) | xc %in% c("", "Prefer not to answer", "Do not know"), NA_integer_,
          fifelse(xc %in% c("No", "no", "0"), 0L, 1L))
}
# symptom-severity items scored 0 (no problem) to 3 (severe problem)
severity <- c("No problem" = 0, "Slight problem" = 1, "Moderate problem" = 2, "Severe problem" = 3,
              "0" = 0, "1" = 1, "2" = 2, "3" = 3)
recode_severity <- function(x) if (is.numeric(x)) x else unname(severity[as.character(x)])
get_fields <- function(ids) sapply(paste0("p", ids), function(cn) raw[[cn]], simplify = FALSE)

pain <- raw[, .(eid)]
binary_items <- c(pain_chronic = "p120019", pain_widespread = "p120021",
                  cond_neuropathy = "p120008", cond_fibromyalgia = "p120009",
                  cond_cfs = "p120010", cond_migraine = "p120016")
for (v in names(binary_items)) pain[, (v) := recode_yn(raw[[binary_items[[v]]]])]

pain[, pain_avg := in_range(raw$p120088, 10)]                       # BPI average pain, 0-10
bpi <- sapply(get_fields(120091:120097), in_range, hi = 10)  # BPI interference items
pain[, pain_interference_mean := fifelse(rowSums(!is.na(bpi)) >= 4,  # at least 4 of 7 answered
                                         rowMeans(bpi, na.rm = TRUE), NA_real_)]
dn4 <- sapply(get_fields(120046:120052), recode_yn)  # 7 DN4 self-report items
pain[, pain_dn4_score := rowSums(dn4)]                              # all 7 required
pain[, pain_dn4_positive := as.integer(pain_dn4_score >= 3)]

site <- sapply(get_fields(120023:120036), site_present_flag)  # 14 body sites
head_face <- site[, c("p120023", "p120024")]
pain[, pain_orofacial := fifelse(rowSums(!is.na(head_face)) == 0, NA_integer_,
                                 as.integer(rowSums(head_face == 1L, na.rm = TRUE) > 0))]
pain[, pain_n_sites := fifelse(rowSums(!is.na(site)) == 0, NA_real_,
                               rowSums(site == 1L, na.rm = TRUE))]

pain[, `:=`(pain_fatigue      = recode_severity(raw$p120040),
            pain_unrefreshed  = recode_severity(raw$p120041),
            pain_cog_symptoms = recode_severity(raw$p120042),
            pain_health_today = in_range(raw$p120103, 100))]

save_section(pain, "Experience of pain", c(
  pain_chronic = "binary", pain_widespread = "binary", pain_avg = "continuous",
  pain_interference_mean = "continuous", pain_dn4_score = "continuous",
  pain_dn4_positive = "binary", pain_orofacial = "binary", pain_n_sites = "continuous",
  cond_neuropathy = "binary", cond_fibromyalgia = "binary", cond_cfs = "binary",
  cond_migraine = "binary", pain_fatigue = "ordinal", pain_unrefreshed = "ordinal",
  pain_cog_symptoms = "ordinal", pain_health_today = "continuous"))
