#!/usr/bin/env Rscript
# Follow-up questionnaire exposures, block 2: food-preference questionnaire (Category 1039).
# Seven taste-quality items enter individually; the other items are summarized by five principal
# components. Runs in the R session of 01 (see its header).
library(data.table)

fpq_ids <- c(
  20600, 20601, 20602, 20603, 20604, 20605, 20606, 20607, 20608, 20609,
  20610, 20611, 20612, 20613, 20615, 20616, 20617, 20618, 20619, 20620,
  20621, 20622, 20623, 20624, 20625, 20626, 20627, 20628, 20629, 20630,
  20631, 20632, 20633, 20634, 20635, 20636, 20637, 20638, 20639, 20640,
  20642, 20643, 20644, 20645, 20646, 20647, 20648, 20649, 20650, 20651,
  20652, 20653, 20654, 20655, 20658, 20659, 20660, 20661, 20662, 20663,
  20664, 20665, 20666, 20667, 20671, 20672, 20673, 20674, 20675, 20676,
  20677, 20678, 20679, 20680, 20681, 20682, 20683, 20684, 20685, 20686,
  20687, 20688, 20689, 20690, 20691, 20692, 20693, 20694, 20695, 20696,
  20697, 20698, 20699, 20700, 20701, 20702, 20703, 20704, 20705, 20706,
  20707, 20708, 20709, 20710, 20711, 20712, 20713, 20714, 20715, 20716,
  20717, 20718, 20719, 20720, 20721, 20722, 20723, 20724, 20725, 20726,
  20727, 20728, 20729, 20730, 20731, 20732, 20734, 20735, 20736, 20737,
  20738, 20739, 20740, 20742, 20743, 20744, 20745, 20746, 20747, 20748)
fpq_desc <- c(
  "adding salt to foods", "aniseed", "apple juice", "apples", "asparagus",
  "aubergine", "avocados", "bacon", "baked/steamed fish", "bananas",
  "barbequed/grilled meat", "beef steak", "beetroot", "bell pepper",
  "biscuits", "bitter foods", "bitter/ale", "black olives", "black pepper",
  "blue cheese", "bolognese sauce", "broad beans", "broccoli", "brown rice",
  "brussel sprouts", "burgers (meat)", "burn of spicy foods", "butter on bread",
  "butternut squash", "cabbage", "cake", "cake icing", "capers", "cauliflower",
  "cereal/granola bar", "cheesecake", "cherries", "chicken",
  "chilli pepper", "chips/french fries", "cod", "coffee with sugar",
  "coffee without sugar", "coriander", "corn flakes", "cream", "croissant",
  "cucumber", "curry", "dairy products", "dark chocolate",
  "diet fizzy drinks", "dried fruit", "eggs", "extra virgin olive oil",
  "fatty foods", "fresh tomatoes", "fried chicken", "fried/battered fish",
  "fruit", "garlic", "gherkins", "globe artichoke", "goat's cheese",
  "grapefruit", "green olives", "haddock", "ham", "hard cheese",
  "herring", "honey", "horseradish/wasabi", "ice cream", "jam",
  "kiwi fruit", "lager", "lamb", "lemons", "lentils/beans", "liver",
  "mackerel", "marzipan", "mayonnaise", "melon", "milk chocolate",
  "mushrooms", "onions", "orange juice", "oranges", "pasta", "pears",
  "pizza", "plain yogurt", "plums", "pollock", "pork chop", "porridge",
  "potato crisps", "potatoes", "prawns", "raw carrots", "red meat",
  "red wine", "regular fizzy drinks", "roast chicken", "salad dressing",
  "salad leaves", "salami", "salmon", "salty foods", "salty pretzels",
  "sardines", "sausages (meat)", "savoury biscuits", "shellfish",
  "skimmed milk", "smoked fish", "soft cheese", "soy sauce", "soya milk",
  "spicy foods", "spinach", "spirits", "strawberries",
  "sweet coffee house drinks", "sweet foods", "tea with sugar",
  "tea without sugar", "tomato ketchup", "tinned tuna", "turnip",
  "vegetables", "vinegar", "whisky", "white bread", "white rice",
  "white wine", "whole grain breakfast cereal", "whole milk", "wholemeal bread")
taste_items <- c(fpq_taste_bitter = 20616, fpq_taste_sweet = 20732, fpq_taste_salty = 20716,
                 fpq_taste_salt_adding = 20600, fpq_taste_fatty = 20659,
                 fpq_taste_spicy = 20727, fpq_taste_spicy_burn = 20627)

# Data-coding 96: 1 = extremely dislike ... 9 = extremely like; "never tried" -> NA
fpq_map <- c("extremely dislike" = 1L, "2" = 2L, "3" = 3L, "4" = 4L,
             "neither like nor dislike" = 5L, "6" = 6L, "7" = 7L, "8" = 8L, "extremely like" = 9L)
fpq <- raw[, .(eid)]
for (id in fpq_ids) fpq[, (paste0("p", id)) := fpq_map[raw[[paste0("p", id)]]]]
for (v in names(taste_items)) fpq[, (v) := get(paste0("p", taste_items[[v]]))]

# PCA of the correlation matrix (pairwise-complete) over the general food items
X <- as.matrix(fpq[, paste0("p", setdiff(fpq_ids, taste_items)), with = FALSE])
col_var <- apply(X, 2, var, na.rm = TRUE)
X <- X[, !is.na(col_var) & col_var >= 1e-10, drop = FALSE]
R <- cor(X, use = "pairwise.complete.obs")
R[is.na(R)] <- 0
eig <- eigen(R)
N_PC <- 5L
loadings <- eig$vectors[, 1:N_PC]
X_std <- sweep(sweep(X, 2, colMeans(X, na.rm = TRUE), "-"), 2, apply(X, 2, sd, na.rm = TRUE), "/")
scores <- X_std %*% loadings
pc_names <- paste0("fpq_pca_pc", 1:N_PC)
for (k in 1:N_PC) fpq[, (pc_names[k]) := scores[, k]]

pc_lab <- paste0("PC", 1:N_PC)
fwrite(data.table(variable = colnames(X), description = fpq_desc[match(colnames(X), paste0("p", fpq_ids))],
                  setnames(as.data.table(round(loadings, 5)), pc_lab)),
       file.path(OUT_DIR, "followup_exwas_food_preferences_pca_loadings.csv"))
var_pct <- eig$values / sum(eig$values) * 100
fwrite(data.table(PC = pc_lab, eigenvalue = round(eig$values[1:N_PC], 3),
                  var_explained_pct = round(var_pct[1:N_PC], 1),
                  cum_var_pct = round(cumsum(var_pct)[1:N_PC], 1)),
       file.path(OUT_DIR, "followup_exwas_food_preferences_pca_variance.csv"))

save_section(fpq, "Food preferences",
             c(setNames(rep("continuous", length(taste_items)), names(taste_items)),
               setNames(rep("derived_continuous", N_PC), pc_names)))
