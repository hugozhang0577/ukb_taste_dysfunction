# Sex-specific k-means cluster id -> subtype letter, subtype names, and the MOFA+ initializations.
SUBTYPE_MAP <- list(
  m = c("1" = "B", "2" = "D", "3" = "A", "4" = "C"),
  f = c("1" = "B", "2" = "A", "3" = "D", "4" = "C")
)
SUBTYPE_MAP_INV <- lapply(SUBTYPE_MAP, function(x) setNames(as.integer(names(x)), x)[c("A", "B", "C", "D")])
SUBTYPE_NAME <- c(A = "Aging frailty", B = "Psychosomatic", C = "Cardiometabolic", D = "Younger idiopathic")

# 21 initializations per sex; primary fit = highest median ARI against the other 20, within each sex
FIT_SEEDS    <- c(20260413L, 10001:10020)
PRIMARY_SEED <- c(m = 10014L, f = 10010L)
other_seeds  <- function(sx) setdiff(FIT_SEEDS, PRIMARY_SEED[[tolower(sx)]])
