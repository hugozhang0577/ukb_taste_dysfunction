# Baseline exposures: clean the baseline questionnaire, physical-measurement and blood-assay
# fields, recode them into continuous, ordinal and binary exposures grouped into exposure
# families, and write one exposure matrix per cohort plus the variable dictionary.
library(data.table)
library(stringr)
FIELDS <- setdiff(trimws(readLines("exwas_baseline/field_list.txt")), c("", "eid"))
setwd(Sys.getenv("PROJECT_DIR", "."))
OUT_DIR <- "output/baseline_exwas/derive"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

D_SES   <- "Demographics & SES";               D_LIFE  <- "Lifestyle & sleep"
D_DIET  <- "Dietary intake & preferences";     D_ANTH  <- "Anthropometric & physiological"
D_ORAL  <- "Oral health";                      D_MENT  <- "General & mental health"
D_PAIN  <- "Sensory function & pain";          D_CLIN  <- "Clinical screening & treatment"
D_BLOOD <- "Blood biochemistry & haematology"; D_ENV   <- "Environmental exposures"
D_HAND  <- "Handedness & laterality";          D_FEM   <- "Reproductive (female-only)"
FAMILIES <- c(D_SES, D_LIFE, D_DIET, D_ANTH, D_ORAL, D_MENT, D_PAIN, D_CLIN, D_BLOOD, D_ENV, D_HAND, D_FEM)

raw <- fread("input/raw/baseline_exwas_table_export.csv", na.strings = c("", "NA"))
pheno <- lapply(setNames(nm = c("group1", "group2", "group3")), function(g)
  fread(sprintf("input/analysis_ready/phenotype_%s.csv", g), na.strings = c("", "NA")))

# 1. Refusal, "do not know" and empty answers -> missing
for (cn in intersect(FIELDS, names(raw))) if (is.character(raw[[cn]]))
  raw[get(cn) %in% c("Do not know", "Prefer not to answer", "Not sure - had a hysterectomy", ""),
      (cn) := NA_character_]

# 2. Repeated measurements
# Blood pressure and pulse: mean of the two automated readings
for (nm in c("dbp_mean", "sbp_mean", "pulse_mean")) {
  f <- c(dbp_mean = 4079, sbp_mean = 4080, pulse_mean = 95)[[nm]]
  raw[, (nm) := rowMeans(.SD, na.rm = TRUE), .SDcols = paste0("p", f, "_i0_a", 0:1)]
  raw[is.nan(get(nm)), (nm) := NA_real_]
}
# Spirometry (FEV1 field 3063, FVC field 3062): best of three blows, FEV1/FVC
for (nm in c("fev1_best", "fvc_best")) {
  cols <- intersect(paste0("p", c(fev1_best = 3063, fvc_best = 3062)[[nm]], "_i0_a", 0:2), names(raw))
  if (length(cols)) {
    raw[, (nm) := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = cols]
    raw[is.infinite(get(nm)), (nm) := NA_real_]
  }
}
if (all(c("fev1_best", "fvc_best") %in% names(raw)))
  raw[, fev1_fvc_ratio := fifelse(fvc_best > 0, fev1_best / fvc_best, NA_real_)]
# Grip strength: the stronger hand
raw[, grip_max := pmax(p46_i0, p47_i0, na.rm = TRUE)]
raw[is.infinite(grip_max), grip_max := NA_real_]

# 3. Multi-response fields -> one indicator per answer + number of answers
multi_response_labels <- list(
  p6138_i0 = c("1" = "College or University degree", "2" = "A levels/AS levels or equivalent",
               "3" = "O levels/GCSEs or equivalent", "4" = "CSEs or equivalent",
               "5" = "NVQ or HND or HNC or equivalent",
               "6" = "Other professional qualifications eg: nursing, teaching",
               "-7" = "None of the above", "-3" = "Prefer not to answer"),
  p6142_i0 = c("1" = "In paid employment or self-employed", "2" = "Retired",
               "3" = "Looking after home and/or family",
               "4" = "Unable to work because of sickness or disability", "5" = "Unemployed",
               "6" = "Doing unpaid or voluntary work", "7" = "Full or part-time student",
               "-7" = "None of the above", "-3" = "Prefer not to answer"),
  p6144_i0 = c("1" = "Eggs or foods containing eggs", "2" = "Dairy products", "3" = "Wheat products",
               "4" = "Sugar or foods/drinks containing sugar", "5" = "I eat all of the above",
               "-3" = "Prefer not to answer"),
  p6145_i0 = c("1" = "Serious illness, injury or assault to yourself",
               "2" = "Serious illness, injury or assault of a close relative",
               "3" = "Death of a close relative", "4" = "Death of a spouse or partner",
               "5" = "Marital separation/divorce", "6" = "Financial difficulties",
               "-7" = "None of the above", "-3" = "Prefer not to answer"),
  p6150_i0 = c("1" = "Heart attack", "2" = "Angina", "3" = "Stroke", "4" = "High blood pressure",
               "-7" = "None of the above", "-3" = "Prefer not to answer"),
  p6159_i0 = c("1" = "Headache", "2" = "Facial pain", "3" = "Neck or shoulder pain", "4" = "Back pain",
               "5" = "Stomach or abdominal pain", "6" = "Hip pain", "7" = "Knee pain",
               "8" = "Pain all over the body", "-7" = "None of the above", "-3" = "Prefer not to answer"))

# Fields exported as numeric codes are first mapped to their answer labels
for (cn in names(multi_response_labels)) {
  vals  <- as.character(raw[[cn]])
  first <- vals[!is.na(vals) & vals != ""][1]
  if (is.na(first) || grepl("[A-Za-z]", first)) next
  lab <- multi_response_labels[[cn]]
  new <- vapply(strsplit(vals, "|", fixed = TRUE), function(parts) {
    m <- lab[parts]; m[is.na(m)] <- parts[is.na(m)]; paste(unname(m), collapse = "|")
  }, character(1))
  new[is.na(vals)] <- NA_character_
  raw[, (cn) := new]
}

decompose_multiresponse <- function(dt, colname, prefix, none_label = "None of the above") {
  skip  <- c(none_label, "Do not know", "Prefer not to answer")
  vals  <- dt[[colname]]
  vals[is.na(vals)] <- ""
  items <- trimws(unique(unlist(strsplit(vals, "\\|"))))
  for (item in items[items != "" & !items %in% skip]) {
    nm <- gsub("_$", "", gsub("_+", "_", paste0(prefix, "_", gsub("[^a-zA-Z0-9]", "_", tolower(item)))))
    dt[, (nm) := as.integer(grepl(item, get(colname), fixed = TRUE))]
  }
  cnt <- paste0(prefix, "_count")
  dt[, (cnt) := str_count(get(colname), "\\|") +
       fifelse(get(colname) != "" & !get(colname) %in% skip, 1L, 0L)]
  dt[get(colname) == "" | is.na(get(colname)), (cnt) := NA_integer_]
  dt[get(colname) %in% none_label, (cnt) := 0L]
  invisible(dt)
}
raw <- decompose_multiresponse(raw, "p6145_i0", "stressor")
raw <- decompose_multiresponse(raw, "p6159_i0", "pain_month")
raw <- decompose_multiresponse(raw, "p6144_i0", "never_eat", none_label = "I eat all of the above")
raw <- decompose_multiresponse(raw, "p6150_i0", "vascular")
raw <- decompose_multiresponse(raw, "p6138_i0", "qualification")
raw <- decompose_multiresponse(raw, "p6142_i0", "employment")

# 4. Answer labels -> numeric codes (unlisted labels -> NA)
recode <- function(x, map) unname(map[as.character(x)])
YES_NO <- c(Yes = 1L, No = 0L)

for (cn in c("p1210_i0", "p2020_i0", "p2188_i0", "p2247_i0", "p2443_i0", "p2956_i0", "p3404_i0",
             "p3414_i0", "p3571_i0", "p3741_i0", "p3773_i0", "p3799_i0", "p4067_i0", "p4501_i0",
             "p4598_i0", "p4631_i0", "p1920_i0"))
  raw[, (paste0(cn, "_bin")) := recode(get(cn), YES_NO)]
raw[, p2453_i0_bin := fcase(grepl("^Yes", p2453_i0), 1L, p2453_i0 == "No", 0L, default = NA_integer_)]

FOOD_FREQ <- c("Never" = 0L, "Less than once a week" = 1L, "Once a week" = 2L,
               "2-4 times a week" = 3L, "5-6 times a week" = 4L, "Once or more daily" = 5L)
for (cn in c("p1329_i0", "p1339_i0", "p1349_i0", "p1359_i0", "p1369_i0", "p1379_i0", "p1389_i0", "p1408_i0"))
  raw[, (paste0(cn, "_ord")) := recode(get(cn), FOOD_FREQ)]

# Amounts recorded as text ("Less than one" = 0.5)
for (cn in c("p1289_i0", "p1299_i0", "p1309_i0", "p1319_i0", "p1438_i0", "p1458_i0", "p1488_i0",
             "p1498_i0", "p1528_i0", "p845_i0", "p1160_i0", "p3680_i0", "p2714_i0", "p3581_i0"))
  raw[, (paste0(cn, "_num")) := suppressWarnings(as.numeric(gsub("Less than one", "0.5", get(cn), fixed = TRUE)))]

SLEEP_FREQ <- c("Never/rarely" = 0L, "Sometimes" = 1L, "Often" = 2L, "Usually" = 3L)
for (cn in c("p1190_i0", "p1200_i0", "p1220_i0")) raw[, (paste0(cn, "_ord")) := recode(get(cn), SLEEP_FREQ)]
PHQ <- c("Not at all" = 0L, "Several days" = 1L, "More than half the days" = 2L, "Nearly every day" = 3L)
for (cn in c("p2050_i0", "p2060_i0")) raw[, (paste0(cn, "_ord")) := recode(get(cn), PHQ)]
SOCIAL_FREQ <- c("No friends/family" = 0L, "Never or almost never" = 0L, "Once every few months" = 1L,
                 "About once a month" = 2L, "About once a week" = 3L, "2-4 times a week" = 4L,
                 "Almost daily" = 5L)
for (cn in c("p1031_i0", "p2110_i0")) raw[, (paste0(cn, "_ord")) := recode(get(cn), SOCIAL_FREQ)]

raw[, salt_added_ord := recode(p1478_i0, c("Never/rarely" = 0L, "Sometimes" = 1L, "Usually" = 2L, "Always" = 3L))]
raw[, health_rating_ord := recode(p2178_i0, c("Excellent" = 4L, "Good" = 3L, "Fair" = 2L, "Poor" = 1L))]
raw[, falls_ord := recode(p2296_i0, c("No falls" = 0L, "Only one fall" = 1L, "More than one fall" = 2L))]
raw[, weight_change_cat := recode(p2306_i0, c("No - weigh about the same" = 0L,
                                              "Yes - gained weight" = 1L, "Yes - lost weight" = -1L))]
raw[, getting_up_ord := recode(p1170_i0, c("Very easy" = 4L, "Fairly easy" = 3L,
                                           "Not very easy" = 2L, "Not at all easy" = 1L))]
raw[, chronotype_ord := recode(p1180_i0, c(
  "Definitely a 'morning' person" = 2L, "More a 'morning' than 'evening' person" = 1L,
  "More an 'evening' than a 'morning' person" = -1L, "Definitely an 'evening' person" = -2L))]
raw[, income_ord := recode(p738_i0, c("Less than 18,000" = 1L, "18,000 to 30,999" = 2L,
                                      "31,000 to 51,999" = 3L, "52,000 to 100,000" = 4L,
                                      "Greater than 100,000" = 5L))]
raw[, diet_variation_ord := recode(p1548_i0, c("Never/rarely" = 0L, "Sometimes" = 1L, "Often" = 2L))]
raw[, hot_drink_temp_ord := recode(p1518_i0, c("Do not drink hot drinks" = 0L, "Warm" = 1L,
                                               "Hot" = 2L, "Very hot" = 3L))]
raw[, tinnitus_ord := recode(p4803_i0, c(
  "No, never" = 0L, "Yes, but not now, but have in the past" = 1L, "Yes, now some of the time" = 2L,
  "Yes, now most of the time" = 3L, "Yes, now all of the time" = 4L))]
raw[, menopause_bin := recode(p2724_i0, YES_NO)]
raw[, prospective_memory_ord := recode(p20018_i0, c(
  "Correct recall on first attempt" = 2L, "Correct recall on second attempt" = 1L,
  "Instruction not recalled, either skipped or incorrect" = 0L))]
raw[, body_size_10_ord := recode(p1687_i0, c("Thinner" = -1L, "About average" = 0L, "Plumper" = 1L))]
raw[, handedness_cat := recode(p1707_i0, c("Right-handed" = "right", "Left-handed" = "left",
                                           "Use both right and left hands equally" = "ambidextrous"))]
raw[, diet_change_bin := recode(p1538_i0, c("No" = 0L, "Yes, because of illness" = 1L,
                                            "Yes, because of other reasons" = 1L))]
raw[, diet_change_illness := as.integer(p1538_i0 == "Yes, because of illness")]

# 5. Biologically implausible values -> NA
bio_ranges <- list(
  p50_i0 = c(130, 210), p21002_i0 = c(30, 250), p48_i0 = c(50, 180), p49_i0 = c(50, 180),
  p46_i0 = c(0, 90), p47_i0 = c(0, 90), sbp_mean = c(70, 250), dbp_mean = c(30, 150),
  pulse_mean = c(25, 200), p23099_i0 = c(1, 75),
  p30600_i0 = c(15, 60), p30690_i0 = c(1, 20), p30700_i0 = c(10, 700), p30710_i0 = c(0, 200),
  p30740_i0 = c(1, 40), p30750_i0 = c(10, 200), p30870_i0 = c(0, 30), p30880_i0 = c(50, 900),
  p30890_i0 = c(0, 400),
  fev1_best = c(0.3, 7), fvc_best = c(0.5, 9), fev1_fvc_ratio = c(0.15, 1.1),
  p20161_i0 = c(0, 200), p22040_i0 = c(0, 20000), p24006 = c(0, 100), p20023_i0 = c(50, 5000),
  p20127_i0 = c(0, 12))
for (v in intersect(names(bio_ranges), names(raw)))
  raw[get(v) < bio_ranges[[v]][1] | get(v) > bio_ranges[[v]][2], (v) := NA]

# 6. Derived variables
raw[, whr := fifelse(p49_i0 > 0, p48_i0 / p49_i0, NA_real_)]
raw[, phq2_score := p2050_i0_ord + p2060_i0_ord]
pain_bin <- paste0(c("p2956_i0", "p3404_i0", "p3414_i0", "p3571_i0", "p3741_i0", "p3773_i0",
                     "p3799_i0", "p4067_i0"), "_bin")
raw[, pain_site_count := rowSums(.SD, na.rm = TRUE), .SDcols = pain_bin]
all_pain_na <- rowSums(!is.na(raw[, ..pain_bin])) == 0
raw[all_pain_na, pain_site_count := NA_integer_]
raw[, srt_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("p20019_i0", "p20021_i0")]
raw[is.nan(srt_mean), srt_mean := NA_real_]

# 7. Exposure registry (dictionary order = column order in the exposure matrices)
dict <- list()
reg <- function(var_name, source_field, domain, var_type, description = "")
  dict[[length(dict) + 1]] <<- data.table(var_name, source_field, domain, var_type,
                                          female_only = domain == D_FEM, description)
add_exposure <- function(name, src, domain, type, desc = "") {
  if (!src %in% names(raw)) return(invisible())   # source not in the extract
  raw[, (name) := raw[[src]]]
  reg(name, src, domain, type, desc)
}

continuous <- list(
  c("height", "p50_i0", D_ANTH, "Standing height (cm)"),
  c("weight", "p21002_i0", D_ANTH, "Weight (kg)"),
  c("waist_circ", "p48_i0", D_ANTH, "Waist circumference (cm)"),
  c("hip_circ", "p49_i0", D_ANTH, "Hip circumference (cm)"),
  c("grip_left", "p46_i0", D_ANTH, "Hand grip strength left (kg)"),
  c("grip_right", "p47_i0", D_ANTH, "Hand grip strength right (kg)"),
  c("body_fat_pct", "p23099_i0", D_ANTH, "Body fat percentage"),
  c("whole_body_fat", "p23100_i0", D_ANTH, "Whole body fat mass (kg)"),
  c("fat_free_mass", "p23101_i0", D_ANTH, "Whole body fat-free mass (kg)"),
  c("trunk_fat_pct", "p23127_i0", D_ANTH, "Trunk fat percentage"),
  c("sbp", "sbp_mean", D_ANTH, "Systolic BP mean of 2 readings (mmHg)"),
  c("dbp", "dbp_mean", D_ANTH, "Diastolic BP mean of 2 readings (mmHg)"),
  c("pulse", "pulse_mean", D_ANTH, "Pulse rate mean of 2 readings (bpm)"),
  c("grip_max", "grip_max", D_ANTH, "Max grip strength L/R (kg)"),
  c("whr", "whr", D_ANTH, "Waist-hip ratio"),
  c("fev1", "fev1_best", D_ANTH, "FEV1 best of 3 blows (L)"),
  c("fvc", "fvc_best", D_ANTH, "FVC best of 3 blows (L)"),
  c("fev1_fvc", "fev1_fvc_ratio", D_ANTH, "FEV1/FVC ratio"),
  c("pack_years", "p20161_i0", D_LIFE, "Pack-years of smoking"),
  c("physical_met", "p22040_i0", D_LIFE, "Physical activity MET min/week"),
  c("neuroticism", "p20127_i0", D_MENT, "Neuroticism score (0\u{2013}12)"),
  c("srt_left", "p20019_i0", D_PAIN, "Speech reception threshold left (dB)"),
  c("srt_right", "p20021_i0", D_PAIN, "Speech reception threshold right (dB)"),
  c("srt_mean", "srt_mean", D_PAIN, "Speech reception threshold mean (dB)"),
  c("fluid_intel", "p20016_i0", D_CLIN, "Fluid intelligence score"),
  c("reaction_time", "p20023_i0", D_CLIN, "Mean reaction time (ms)"),
  c("n_medications", "p137_i0", D_MENT, "Number of medications taken"),
  c("pm25", "p24006", D_ENV, "PM2.5 air pollution 2010 (ug/m3)"),
  c("age_menarche", "p2714_i0_num", D_FEM, "Age at menarche (years)"),
  c("age_menopause", "p3581_i0_num", D_FEM, "Age at menopause (years)"))
for (v in continuous) add_exposure(v[1], v[2], v[3], "continuous", v[4])

biochem <- c(albumin = 30600, alp = 30610, alt = 30620, apoa = 30630, apob = 30640, ast = 30650,
             direct_bili = 30660, urea = 30670, calcium = 30680, cholesterol = 30690,
             creatinine = 30700, crp = 30710, cystatin_c = 30720, ggt = 30730, glucose = 30740,
             hba1c = 30750, hdl = 30760, igf1 = 30770, ldl = 30780, lpa = 30790,
             oestradiol = 30800, phosphate = 30810, rheumatoid_f = 30820, shbg = 30830,
             total_bili = 30840, testosterone = 30850, total_protein = 30860,
             triglycerides = 30870, urate = 30880, vitamin_d = 30890)
for (nm in names(biochem))
  add_exposure(nm, paste0("p", biochem[[nm]], "_i0"), D_BLOOD, "continuous", paste("Serum", nm))
haem <- c(wbc = 30000, rbc = 30010, haemoglobin = 30020, haematocrit = 30030, mcv = 30040,
          mch = 30050, mchc = 30060, rdw = 30070, platelet = 30080, platelet_crit = 30090,
          mpv = 30100, platelet_dw = 30110, lymphocyte = 30120, monocyte = 30130,
          neutrophil = 30140, eosinophil = 30150, basophil = 30160, nrbc = 30170,
          lymphocyte_pct = 30180, monocyte_pct = 30190, neutrophil_pct = 30200,
          eosinophil_pct = 30210, basophil_pct = 30220, nrbc_pct = 30230, retic_pct = 30240,
          retic_count = 30250, mean_retic_vol = 30260, mean_spher_vol = 30270,
          immature_retic = 30280, hls_retic_pct = 30290, hls_retic = 30300)
for (nm in names(haem))
  add_exposure(nm, paste0("p", haem[[nm]], "_i0"), D_BLOOD, "continuous", paste("Blood", nm))

coded <- list(
  c("snoring", "p1210_i0_bin", D_LIFE, "binary"),
  c("loneliness", "p2020_i0_bin", D_MENT, "binary"),
  c("longstanding_ill", "p2188_i0_bin", D_MENT, "binary"),
  c("hearing_diff", "p2247_i0_bin", D_PAIN, "binary"),
  c("diabetes_sr", "p2443_i0_bin", D_CLIN, "binary"),
  c("cancer_sr", "p2453_i0_bin", D_CLIN, "binary"),
  c("general_pain_3mo", "p2956_i0_bin", D_PAIN, "binary"),
  c("neck_pain_3mo", "p3404_i0_bin", D_PAIN, "binary"),
  c("hip_pain_3mo", "p3414_i0_bin", D_PAIN, "binary"),
  c("back_pain_3mo", "p3571_i0_bin", D_PAIN, "binary"),
  c("stomach_pain_3mo", "p3741_i0_bin", D_PAIN, "binary"),
  c("knee_pain_3mo", "p3773_i0_bin", D_PAIN, "binary"),
  c("headache_3mo", "p3799_i0_bin", D_PAIN, "binary"),
  c("facial_pain_3mo", "p4067_i0_bin", D_PAIN, "binary"),
  c("fam_nonacc_death", "p4501_i0_bin", D_CLIN, "binary"),
  c("ever_depressed", "p4598_i0_bin", D_MENT, "binary"),
  c("ever_unenthu", "p4631_i0_bin", D_MENT, "binary"),
  c("mood_swings", "p1920_i0_bin", D_MENT, "binary"),
  c("menopause", "menopause_bin", D_FEM, "binary"),
  c("diet_change", "diet_change_bin", D_DIET, "binary"),
  c("diet_change_ill", "diet_change_illness", D_DIET, "binary"),
  c("health_rating", "health_rating_ord", D_MENT, "ordinal"),
  c("falls", "falls_ord", D_MENT, "ordinal"),
  c("weight_change", "weight_change_cat", D_MENT, "ordinal"),
  c("getting_up", "getting_up_ord", D_LIFE, "ordinal"),
  c("chronotype", "chronotype_ord", D_LIFE, "ordinal"),
  c("nap_daytime", "p1190_i0_ord", D_LIFE, "ordinal"),
  c("insomnia", "p1200_i0_ord", D_LIFE, "ordinal"),
  c("daytime_dozing", "p1220_i0_ord", D_LIFE, "ordinal"),
  c("depressed_mood", "p2050_i0_ord", D_MENT, "ordinal"),
  c("unenthusiasm", "p2060_i0_ord", D_MENT, "ordinal"),
  c("phq2_score", "phq2_score", D_MENT, "ordinal"),
  c("friend_visits", "p1031_i0_ord", D_MENT, "ordinal"),
  c("able_confide", "p2110_i0_ord", D_MENT, "ordinal"),
  c("income", "income_ord", D_SES, "ordinal"),
  c("salt_added", "salt_added_ord", D_DIET, "ordinal"),
  c("diet_variation", "diet_variation_ord", D_DIET, "ordinal"),
  c("hot_drink_temp", "hot_drink_temp_ord", D_DIET, "ordinal"),
  c("tinnitus", "tinnitus_ord", D_PAIN, "ordinal"),
  c("prosp_memory", "prospective_memory_ord", D_CLIN, "ordinal"),
  c("body_size_10", "body_size_10_ord", D_ANTH, "ordinal"),
  c("pain_site_count", "pain_site_count", D_PAIN, "ordinal"))
for (v in coded) add_exposure(v[1], v[2], v[3], v[4])

food_freq <- c(oily_fish_freq = "p1329_i0", nonoily_fish_freq = "p1339_i0",
               processed_meat_freq = "p1349_i0", poultry_freq = "p1359_i0", beef_freq = "p1369_i0",
               lamb_freq = "p1379_i0", pork_freq = "p1389_i0", cheese_freq = "p1408_i0")
for (nm in names(food_freq)) add_exposure(nm, paste0(food_freq[[nm]], "_ord"), D_DIET, "ordinal")
amounts <- c(cooked_veg = "p1289_i0", salad_veg = "p1299_i0", fresh_fruit = "p1309_i0",
             dried_fruit = "p1319_i0", bread_intake = "p1438_i0", cereal_intake = "p1458_i0",
             tea_intake = "p1488_i0", coffee_intake = "p1498_i0", water_intake = "p1528_i0",
             age_last_meat = "p3680_i0")
for (nm in names(amounts)) add_exposure(nm, paste0(amounts[[nm]], "_num"), D_DIET, "continuous")
add_exposure("sleep_hours", "p1160_i0_num", D_LIFE, "continuous", "Sleep duration (hours)")
add_exposure("age_education", "p845_i0_num", D_SES, "continuous", "Age completed full-time education (years)")

# Categorical fields: one indicator per observed answer
dummy_code <- function(colname, prefix, domain) {
  x <- raw[[colname]]
  for (lv in sort(unique(x[!is.na(x)]))) {
    nm <- paste0(prefix, "_", gsub("^_|_$", "", gsub("_+", "_", gsub("[^a-zA-Z0-9]", "_", tolower(lv)))))
    raw[, (nm) := fifelse(!is.na(x) & x == lv, 1L, fifelse(is.na(x), NA_integer_, 0L))]
    reg(nm, colname, domain, "binary", paste0(prefix, ": ", lv))
  }
}
dummy_code("p1418_i0", "milk_type", D_DIET)
dummy_code("p1428_i0", "spread_type", D_DIET)
dummy_code("p1448_i0", "bread_type", D_DIET)
dummy_code("p1468_i0", "cereal_type", D_DIET)
dummy_code("p1508_i0", "coffee_type", D_DIET)
dummy_code("p2654_i0", "spread_detail", D_DIET)
dummy_code("handedness_cat", "handedness", D_HAND)

multi_domain <- c(stressor = D_MENT, pain_month = D_PAIN, never_eat = D_DIET, vascular = D_CLIN,
                  qualification = D_SES, employment = D_SES)
for (pfx in names(multi_domain))
  for (dc in grep(paste0("^", pfx, "_"), names(raw), value = TRUE))
    reg(dc, dc, multi_domain[[pfx]], if (grepl("_count$", dc)) "ordinal" else "binary",
        paste("Decomposed from", pfx))

# Oral health and treatment covariates from the cohort phenotype files, tested here as exposures
from_pheno <- list(
  c("oral_bleeding_gums_baseline", D_ORAL, "binary", "Bleeding gums at baseline"),
  c("oral_painful_gums_baseline", D_ORAL, "binary", "Painful gums at baseline"),
  c("oral_loose_teeth_baseline", D_ORAL, "binary", "Loose teeth at baseline"),
  c("oral_toothache_baseline", D_ORAL, "binary", "Toothache at baseline"),
  c("oral_mouth_ulcers_baseline", D_ORAL, "binary", "Mouth ulcers at baseline"),
  c("oral_dentures_baseline", D_ORAL, "binary", "Dentures at baseline"),
  c("periodontal_indicator_baseline", D_ORAL, "binary", "Periodontal composite indicator at baseline"),
  c("any_oral_problem_baseline", D_ORAL, "binary", "Any oral problem composite at baseline"),
  c("surg_taste_affecting_baseline", D_CLIN, "binary", "Taste-affecting surgery before baseline"),
  c("surg_non_taste_affecting_baseline", D_CLIN, "binary",
    "Non-taste-affecting surgery before baseline (negative control)"),
  c("surg_tooth_extraction_baseline", D_CLIN, "binary", "Tooth extraction surgery before baseline"),
  c("taste_drug_user_baseline_numeric", D_CLIN, "binary", "Taste-affecting medication use at baseline"),
  c("egfr", D_BLOOD, "continuous", "Estimated GFR (derived from creatinine + cystatin C)"))
for (v in from_pheno) reg(v[1], "phenotype", v[2], v[3], v[4])

dict <- rbindlist(dict)
dict[, source := domain]
fwrite(dict, file.path(OUT_DIR, "baseline_exwas_variable_dictionary.csv"))

# 8. One exposure matrix per cohort: participants with a defined outcome, families in FAMILIES order
for (g in names(pheno)) {
  ph <- pheno[[g]][!is.na(taste_2w_strict)]
  parts <- lapply(FAMILIES, function(fam) {
    d  <- dict[source == fam]
    rv <- d[source_field != "phenotype", var_name]
    pv <- intersect(d[source_field == "phenotype", var_name], names(ph))
    x  <- if (length(rv)) raw[eid %in% ph$eid, c("eid", rv), with = FALSE] else ph[, .(eid)]
    if (length(pv)) x <- merge(x, ph[, c("eid", pv), with = FALSE], by = "eid", all.x = TRUE)
    x[, c("eid", intersect(d$var_name, names(x))), with = FALSE]
  })
  fwrite(Reduce(function(a, b) merge(a, b, by = "eid", all = TRUE), parts),
         file.path(OUT_DIR, sprintf("exwas_baseline_%s.csv", g)))
}
