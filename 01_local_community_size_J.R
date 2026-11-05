################################################################################
## Reichmann et al.
## Artificial light at night reorganizes habitat connectivity 
## and erodes landscape capacity to sustain biodiversity across South America
################################################################################
# This script is part of the code accompanying the manuscript above.
#
# 01 — Local community size (J)
# Calculates freshwater and LAND J from habitat areas. 
# Creates M.SA.J with all freshwater and LAND J columns and saves it as M.SA.J.RData.

# Calculate freshwater and LAND J from habitat areas.
setwd(".")
load("Paisaje_SA_EcoR_M_SA_J_LP.RData")

J.max <- 400
area.max.sa <- max(M.SA$FRESHWATER)
b.ef <- 0.5

# Freshwater
M.SA.J <- M.SA
M.SA.J$J.freshwater <- ceiling((J.max / area.max.sa^b.ef) * M.SA$FRESHWATER^b.ef)
M.SA.J$J.efimero <- ceiling((J.max / area.max.sa^b.ef) * M.SA$FRESHWATER.efimeros^b.ef)
M.SA.J$J.temporal <- ceiling((J.max / area.max.sa^b.ef) * M.SA$FRESHWATER.temporales^b.ef)
M.SA.J$J.permanente <- ceiling((J.max / area.max.sa^b.ef) * M.SA$FRESHWATER.permanentes^b.ef)
M.SA.J$CENTROID_Y_3 <- round(M.SA.J$CENTROID_Y, 3)
M.SA.J$CENTROID_X_3 <- round(M.SA.J$CENTROID_X, 3)

# LAND habitat combinations and J are added directly to M.SA.J.
M.SA.J$LAND.efimeros <- M.SA$LAND + M.SA$FRESHWATER.temporales + M.SA$FRESHWATER.permanentes
M.SA.J$LAND.temporales <- M.SA$LAND + M.SA$FRESHWATER.efimeros + M.SA$FRESHWATER.permanentes
M.SA.J$LAND.permanentes <- M.SA$LAND + M.SA$FRESHWATER.efimeros + M.SA$FRESHWATER.temporales
M.SA.J$J.LAND <- (400 / 100^0.5) * sqrt(M.SA$LAND)
M.SA.J$J.LAND.efimeros <- (400 / 100^0.5) * sqrt(M.SA.J$LAND.efimeros)
M.SA.J$J.LAND.temporales <- (400 / 100^0.5) * sqrt(M.SA.J$LAND.temporales)
M.SA.J$J.LAND.permanentes <- (400 / 100^0.5) * sqrt(M.SA.J$LAND.permanentes)

save(M.SA.J, file = "M.SA.J.RData")
