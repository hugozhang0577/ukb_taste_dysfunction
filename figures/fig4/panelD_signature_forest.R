# =============================================================================
# Fig 4 Panel D — subtype signature forest
#
# For each subtype, the standardised mean difference of eight variables against
# the whole-cohort mean, with a 95% CI. z is computed WITHIN each sex and then
# pooled, so the sex-dimorphic variables (grip strength, SHBG) are not driven by
# each subtype's sex composition.
#
# Descriptive only: most of these variables entered the clustering, so the panel
# describes the partition rather than testing it. No p-values or significance
# stars are drawn (post-clustering testing on the clustering inputs would be
# circular; Zhang 2019 Cell Systems). Independent characterisation is Figure 5.
#
# The plotted variables are pre-specified, not ranked out of the data — see the
# selection block below.
# =============================================================================
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

P7 <- "output/subtyping"; ML <- "output/ml_ready"
MF <- "output/feature_manifest"; OUT <- "output/figures/fig4"
SUBTYPE_MAP <- list(m=c("1"="D","2"="A","3"="B","4"="C"), f=c("1"="C","2"="B","3"="A","4"="D"))
SUB_COL <- c(A="#4E79A7", B="#E15759", C="#F28E2B", D="#59A14F")
SUB_TTL <- c(A="A   Aging frailty", B="B   Psychosomatic",
             C="C   Cardiometabolic", D="D   Young idiopathic")

manifest <- fread(file.path(MF, "master_feature_manifest_final.csv"))
is_cont <- function(src,vt,phe){
  if(!is.na(src)&&src%in%c("PWAS","MWAS")) return(TRUE)
  if(!is.na(phe)&&nzchar(as.character(phe))) return(FALSE)
  if(tolower(ifelse(is.na(vt),"",vt))%in%c("binary","derived_binary")) return(FALSE)
  TRUE }
manifest[, cont := mapply(is_cont, source_analysis, var_type, phecode)]
src_of <- setNames(manifest$source_analysis, manifest$feature_id)
EXCL <- c("sex","assess_centre_id","eid","age_baseline","years_baseline_to_taste",
          "grip_left","grip_right","whole_body_fat","body_fat_pct","trunk_fat_pct",
          "smoking","drink","taste_2w_strict","taste_any")  # drop collinear adiposity
cont_feats <- setdiff(manifest[cont==TRUE, feature_id], EXCL)
cont_feats <- cont_feats[!grepl("^(PC[0-9]+|smell|taste_(basic|2w|4w|any))", cont_feats)]

# --- within-sex z, then pooled across sex -------------------------------------
dat <- as.data.table(readRDS(file.path(ML,"group1_full.rds")))
cases_all <- dat[taste_2w_strict==1]
zl <- list()
for(sx in c("m","f")){
  cl <- readRDS(file.path(P7,"clusters",sprintf("cluster_assignments_g1_%s_k4.rds",sx)))
  sub <- data.table(eid=as.integer(rownames(cl$Z_active)),
                    subtype=SUBTYPE_MAP[[sx]][as.character(cl$final_labels)])
  cs <- cases_all[sex==(if(sx=="m")1L else 0L)][sub, on="eid", nomatch=0]
  for(fid in intersect(cont_feats,names(cs))){
    x <- suppressWarnings(as.numeric(cs[[fid]]))
    if(sum(!is.na(x))<50 || length(unique(x[!is.na(x)]))<6) next
    zl[[paste(sx,fid)]] <- data.table(feature_id=fid, sex=sx, subtype=cs$subtype,
                                      eid=cs$eid, z=(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
  }
}
Z <- rbindlist(zl)[!is.na(z)&!is.na(subtype)]
pooled <- Z[, .(m=mean(z), se=sd(z)/sqrt(.N), n=.N), by=.(feature_id,subtype)]

# --- view + legibility: whitelist-based (unrecognized -> "Other" -> dropped) --
ADI         <- c("BMI","waist_circ","hip_circ","whr","weight")   # general/central adiposity
FUNC        <- c("grip_max")                                     # physical function
LIPID_CLIN  <- c("hdl","ldl","ldl_direct","hdl_direct","cholesterol","triglycerides","apoa","apob","chol")
LIPID_NMR   <- c("HDL_C","ApoA1","Total_TG","VLDL_C","LDL_size")            # legible NMR lipids
NMR_METAB   <- c("GlycA","Albumin","Citrate","Gln","Val","DHA","DHA_pct","PUFA_pct","MUFA","MUFA_pct","SFA","LA_pct")
BIOCHEM     <- c("urate","glucose","hba1c","creatinine","cystatin_c","igf1","shbg","crp","calcium",
                 "phosphate","vitamin_d","alt","ast","ggt","alp","testosterone","total_protein",
                 "albumin","direct_bilirubin","total_bilirubin","alkaline_phosphatase")
AFFECT      <- c("phq9_total","phq2_score","gad7_total","neuroticism")
view_of <- function(fid){
  s <- src_of[[fid]]
  if(fid=="age") return("Demographic")
  if(fid %in% ADI) return("Adiposity")
  if(fid %in% FUNC) return("Physical function")
  if(fid %in% AFFECT) return("Affective")
  if(grepl("^phq15|^pain", fid)) return("Somatic/pain")
  if(grepl("^cog_", fid)) return("Cognition")
  if(fid %in% LIPID_CLIN || fid %in% LIPID_NMR) return("Lipids")
  if(fid %in% NMR_METAB) return("NMR metabolite")
  if(fid %in% BIOCHEM) return("Clinical biochemistry")
  if(!is.na(s)&&s=="PWAS") return("Protein")
  if(grepl("dx_count", fid)) return("Disease burden")
  if(fid=="townsend") return("Sociodemographic")
  "Other" }                                          # unrecognized -> dropped

# --- labels -------------------------------------------------------------------
LAB <- c(age="Age", BMI="BMI", waist_circ="Waist circ.", hip_circ="Hip circ.",
         whr="Waist-to-hip ratio", body_fat_pct="Body fat %", trunk_fat_pct="Trunk fat %",
         weight="Weight", grip_max="Grip strength", phq9_total="Depression (PHQ-9)",
         phq2_score="Depression (PHQ-2)", gad7_total="Anxiety (GAD-7)", neuroticism="Neuroticism",
         phq15_somatic_total="Somatic symptoms (PHQ-15)", pain_interference_mean="Pain interference",
         pain_health_today="Pain: health today", pain_month_count="Pain sites (count)",
         cog_symbol_digit="Cognition: symbol-digit", cog_tmt_a="Cognition: trail A",
         cog_tmt_b="Cognition: trail B", cog_matrix_reasoning="Cognition: matrix reasoning",
         urate="Urate", glucose="Glucose",
         hba1c="HbA1c", creatinine="Creatinine", cystatin_c="Cystatin C", igf1="IGF-1",
         shbg="SHBG", crp="CRP", total_dx_count="Total diagnoses", dx_count_psych="Psychiatric diagnoses (count)",
         townsend="Deprivation (Townsend)",
         GlycA="GlycA (inflammation)", HDL_C="HDL cholesterol", Total_TG="Triglycerides",
         ApoA1="ApoA1", Albumin="Albumin", VLDL_C="VLDL cholesterol", LDL_size="LDL size",
         hdl="HDL cholesterol", ldl="LDL cholesterol", ldl_direct="LDL cholesterol",
         cholesterol="Total cholesterol", triglycerides="Triglycerides",
         apoa="Apolipoprotein A", apob="Apolipoprotein B",
         Citrate="Citrate", Gln="Glutamine", Val="Valine", DHA="Omega-3 (DHA)",
         DHA_pct="Omega-3 % (DHA)", PUFA_pct="PUFA %", MUFA_pct="MUFA %", MUFA="MUFA",
         SFA="Saturated FA", LA_pct="Linoleic acid %",
         pon3="PON3", cfb="Complement factor B", lep="Leptin", ccl19="CCL19",
         ctsd="Cathepsin D", orm1="Orosomucoid (AGP)", a1bg="A1BG", apcs="Serum amyloid P",
         f9="Coagulation factor IX", oxt="Oxytocin", inhbc="Inhibin C", enpp6="ENPP6",
         ces1="Carboxylesterase 1", cntn5="Contactin-5", wfikkn2="WFIKKN2", igsf9="IGSF9", palm2="PALM2")
mklab <- function(id) ifelse(id %in% names(LAB), LAB[id],
                             ifelse(src_of[id] %in% "PWAS", toupper(id), gsub("_"," ",id)))

# --- the plotted variables ----------------------------------------------------
# Eight per subtype, pre-specified rather than ranked out of the data: their union
# is exactly the variable set of main-text Table 2, whose footnote points back to
# this panel. Age and BMI are shared by all four panels so the rows align; the
# remaining six per subtype span distinct views (at most two adiposity measures,
# collinear body-fat percentages excluded).
CORE    <- c("age","BMI")                                 # universal core (triangles), all panels
ANCHORS <- c("age","BMI","grip_max")                      # pinned to top rows in listed order
                                                          # (grip aligns across A & D where present)
PANEL_VARS <- list(
  A = c("age","BMI","waist_circ","grip_max","hdl","shbg","urate","GlycA"),
  B = c("age","BMI","phq9_total","gad7_total","neuroticism","phq15_somatic_total",
        "pain_interference_mean","dx_count_psych"),
  C = c("age","BMI","GlycA","hba1c","triglycerides","hdl","urate","cystatin_c"),
  D = c("age","BMI","cog_symbol_digit","grip_max","igf1","cystatin_c","GlycA","hba1c"))
fin <- rbindlist(lapply(names(PANEL_VARS), function(s){
  d <- pooled[subtype==s & feature_id %in% PANEL_VARS[[s]]]
  miss <- setdiff(PANEL_VARS[[s]], d$feature_id)
  if(length(miss)) cat("  [not in ml_ready]", s, ":", paste(miss, collapse=", "), "\n")
  d }))
fin[, `:=`(lo=m-1.96*se, hi=m+1.96*se, label=mklab(feature_id),
           view=vapply(feature_id, view_of, ""), core=feature_id %in% CORE)]

cat("\n===== plotted variables (Age+BMI core, then the per-subtype tail) =====\n")
for(s in c("A","B","C","D")){
  cat(sprintf("\n-- %s  (views: %s) --\n", SUB_TTL[[s]],
              paste(sort(unique(fin[subtype==s, view])), collapse=", ")))
  print(fin[subtype==s][order(-abs(m)), .(label, view, z=round(m,2), n)])
}

# --- plot 1x4 strip; core pinned top (Age row1, BMI row2), tail by z below ----
XL <- c(-1.75, 2.4)
mk <- function(dt, s){
  d <- dt[subtype==s]; col <- SUB_COL[[s]]
  anch_fid <- intersect(ANCHORS, d$feature_id)           # ANCHORS order: age, BMI, grip
  anch_lab <- d[match(anch_fid, feature_id), label]      # top block (age, BMI[, grip])
  tail_lab <- d[!(feature_id %in% anch_fid)][order(m), label]  # ascending -> most neg bottom
  d[, label := factor(label, levels = c(tail_lab, rev(anch_lab)))]  # top: age, BMI, grip
  ggplot(d, aes(m,label)) +
    geom_vline(xintercept=0, linewidth=0.35, colour="grey60") +
    geom_errorbarh(aes(xmin=lo,xmax=hi), height=0.14, linewidth=0.5, colour=col) +
    geom_point(aes(shape=core), size=1.9, colour=col) +
    geom_text(aes(label=sprintf("%+.2f",m)), hjust=ifelse(d$m>=0,-0.28,1.28),
              size=2.2, colour="grey30") +
    scale_shape_manual(values=c("FALSE"=16,"TRUE"=17), guide="none") +   # core = triangle
    scale_x_continuous(limits=XL, breaks=c(-1,0,1)) +
    labs(title=SUB_TTL[[s]], x=NULL, y=NULL) +
    theme_bw(base_size=10) +
    theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_line(colour="grey93"),
          panel.border=element_rect(colour=col, linewidth=0.8),
          plot.title=element_text(face="bold", size=10, colour=col, margin=margin(b=1)),
          axis.text.y=element_text(size=7.6, colour="grey15"),
          axis.text.x=element_text(size=7.8), plot.margin=margin(2,5,1,2))
}
P <- (mk(fin,"A")|mk(fin,"B")|mk(fin,"C")|mk(fin,"D"))
xlab <- wrap_elements(full = grid::textGrob("Standardized mean difference vs cohort (z)",
                                            gp = grid::gpar(fontsize = 10, col = "grey20")))
Pf <- P / xlab + plot_layout(heights = c(1, 0.05))
ggsave(file.path(OUT,"fig4_panelD_signature_forest.png"), Pf, width=13.5, height=3.15, dpi=300, bg="white")
ggsave(file.path(OUT,"fig4_panelD_signature_forest.pdf"), Pf, width=13.5, height=3.15)
cat("\nwrote fig4_panelD_signature_forest.{png,pdf}\n")
