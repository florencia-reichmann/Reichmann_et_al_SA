################################################################################
## Reichmann et al.
## Artificial light at night reorganizes habitat connectivity 
## and erodes landscape capacity to sustain biodiversity across South America
################################################################################
# This script is part of the code accompanying the manuscript above.
#
# 02 — Light pollution
# Joins the light pollution input table to M.SA.J and saves the combined
# matrix M.SA.J.LP in the input RData file. Columns are selected by name.

# Join LP by coordinate keys; all column selections use names.
setwd(".")
load("Paisaje_SA_EcoR_M_SA_J_LP.RData")
load("M.SA.J.RData")

LPF <- LP_grid
LPF <- LPF[, c("CENTROID_X", "CENTROID_Y", "magnitud", "_mean", "_median")]
names(LPF)[names(LPF) == "_mean"] <- "LPF_mean"
names(LPF)[names(LPF) == "_median"] <- "LPF_median"

# Unir por coordenadas redondeadas
LPF$CENTROID_X_3 <- round(LPF$CENTROID_X, 3)
LPF$CENTROID_Y_3 <- round(LPF$CENTROID_Y, 3)
M.SA.J$CENTROID_X_3 <- round(M.SA.J$CENTROID_X, 3)
M.SA.J$CENTROID_Y_3 <- round(M.SA.J$CENTROID_Y, 3)
LPF <- LPF[, c("CENTROID_X_3", "CENTROID_Y_3", "magnitud", "LPF_mean", "LPF_median")]
M.SA.J.LP <- merge(M.SA.J, LPF, by = c("CENTROID_X_3", "CENTROID_Y_3"), all.x = TRUE)
M.SA.J.LP <- M.SA.J.LP[, setdiff(names(M.SA.J.LP), c("PageName", "CENTROID_X_3", "CENTROID_Y_3"))]

# Ajustar LP
M.SA.J.LP$LPF_mean[M.SA.J.LP$LPF_mean > 100] <- 97.7555
M.SA.J.LP$LPF_mean_est <- pmin(M.SA.J.LP$LPF_mean, 5) / 5
M.SA.J.LP$magnitud[M.SA.J.LP$magnitud < 10] <- NA
M.SA.J.LP$LP_magnitud_est <- (22 - M.SA.J.LP$magnitud) / (22 - 11.4301)

save(M.SA, limites_ecoregions, LP_grid, M.SA.J, M.SA.J.LP, file = "Paisaje_SA_EcoR_M_SA_J_LP.RData")
