# =============================================================================
# _subtype_map.R
#
# Single source of truth for the subtype letter <-> cluster ID mapping.
# source()'d by every subtype-aware downstream script (06_fullpanel_omics /
# _tables_k4 / evidence-layer builders). k-means assigns arbitrary cluster IDs
# (1-4); this file fixes which ID corresponds to which biological subtype per sex.
#
# How to (re-)derive the mapping (three independent crosschecks):
#   (a) cluster profiles (04_cluster_discriminators) per cluster_id:
#         older + low BMI + high vitamin D            -> A (aging frailty)
#         high PHQ-9 / GAD-7 / neuroticism            -> B (psychosomatic)
#         high BMI + urate + waist                    -> C (cardiometabolic)
#         youngest + no extreme feature               -> D (young idiopathic)
#   (b) matched pairs (subtype_reproducibility/01): each M cluster's best-matched F cluster inherits its
#         letter via (a).
#   (c) cluster-size proportions cross-checked between sexes.
# =============================================================================

# cluster_id (string) -> subtype letter
SUBTYPE_MAP <- list(
  m = c("1" = "D", "2" = "A", "3" = "B", "4" = "C"),
  f = c("1" = "C", "2" = "B", "3" = "A", "4" = "D")
)

# subtype letter -> cluster_id (integer)
SUBTYPE_MAP_INV <- list(
  m = c(A = 2, B = 3, C = 4, D = 1),
  f = c(A = 3, B = 2, C = 1, D = 4)
)

# ---- self-consistency assertion (catches typos / wrong edits) ---------------
local({
  for (sx in c("m", "f")) {
    fwd <- SUBTYPE_MAP[[sx]]; inv <- SUBTYPE_MAP_INV[[sx]]
    for (L in names(inv)) {
      cid <- as.character(inv[[L]])
      if (is.na(fwd[cid]) || fwd[cid] != L)
        stop(sprintf("_subtype_map.R inconsistency: sex=%s letter=%s", sx, L))
    }
    stopifnot(setequal(unname(fwd), names(inv)))
    stopifnot(setequal(names(fwd), as.character(unname(inv))))
  }
})
