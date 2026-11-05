################################################################################
## Reichmann et al.
## Artificial light at night reorganizes habitat connectivity 
## and erodes landscape capacity to sustain biodiversity across South America
################################################################################
# This script is part of the code accompanying the manuscript above.
#
# 03 — Networks and centrality metrics
# Builds connectivity networks and calculates in-degree, out-degree and
# betweenness centrality, with and without considering the impact of light pollution.

setwd(".")
load("Paisaje_SA_EcoR_M_SA_J_LP.RData")


# Visualiación de estructura inicial de la matriz de datos
head(M.SA.J.LP) 
table(M.SA.J.LP$EcoR_ID)
M.SA.J.LP$LP <- M.SA.J.LP$LPF_mean_est # LP normalizado y con maximo en 5

unique(M.SA.J.LP$EcoR_Name)
unique(M.SA.J.LP$EcoR_ID)

# Carga de librerías
library(geosphere)  # Para cálculo de distancias
library(sp)         # Manejo de datos espaciales
library(stars)      # Manipulación de datos geométricos
library(sf)
library(sna)        # Para medidas de centralidad en redes

columnas_base_red <- c(
  "CENTROID_X", "CENTROID_Y", "UTM_sect", "UTM_x", "UTM_y",
  "COUNTRY", "ISO", "EcoR_ID", "EcoR_Name", "PAGENAME",
  "MAR", "LAND", "LP")


# ids <-350 # one EcoR_ID to test
ids <-301:352
ids <- setdiff(ids, c(348, 351)) # 348 va por separado en Malvinas y Patagonia - 351 Juan Fernandez 1 celda
A_list <- list()


###### FOR _ D50.4 _ Q0 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q0.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q0.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", id), M.final.SA.LP.Q0.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", id), M.final.SA.LP.Q0.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}
#save.image("Paisaje_SA_LP_Q0_D50.4_EcoR_ID_345_352.RData")



###### FOR _ D50.4 _ Q025 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0.25  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q025 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q025 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas de M.LAND.D50.4 y M.Permanentes.D50.4
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q025.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q025.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q025.D50.4.PERM.ID_", id), M.final.SA.LP.Q025.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q025.D50.4.LAND.ID_", id), M.final.SA.LP.Q025.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}
#save.image("Paisaje_SA_LP_Q025_D50.4_EcoR_ID_345_352.RData")
 

###### FOR _ D50.4 _ Q1 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <-1  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q1 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q1 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q1.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q1.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q1.D50.4.PERM.ID_", id), M.final.SA.LP.Q1.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q1.D50.4.LAND.ID_", id), M.final.SA.LP.Q1.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}
#save.image("Paisaje_SA_LP_Q1_D50.4_EcoR_ID_345_352.RData")



###### FOR _ D50.4 _ Q5 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 5  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q5 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q5 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q5.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q5.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_", id), M.final.SA.LP.Q5.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_", id), M.final.SA.LP.Q5.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}
#save.image("Paisaje_SA_LP_Q5_D50.4_EcoR_ID_345_352.RData")






































################################ MALVINAS ##### 
ids <-348
A_list <- list()
###### FOR _ D50.4 _ Q0 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  #A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q0.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q0.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_M", id), M.final.SA.LP.Q0.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_M", id), M.final.SA.LP.Q0.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}
#save.image("Paisaje_SA_LP_Q0_D50.4_EcoR_ID_345.RData")
#save.image("Paisaje_SA_LP_Q0_D50.4_EcoR_ID_345.RData")

###### FOR _ D50.4 _ Q025 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  #A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0.25  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q025 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q025 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas de M.LAND.D50.4 y M.Permanentes.D50.4
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q025.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q025.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q025.D50.4.PERM.ID_M", id), M.final.SA.LP.Q025.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q025.D50.4.LAND.ID_M", id), M.final.SA.LP.Q025.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}
#save.image("Paisaje_SA_LP_Q025_D50.4_EcoR_ID_345.RData")
#save.image("Paisaje_SA_LP_Q025_D50.4_EcoR_ID_345.RData")

###### FOR _ D50.4 _ Q1 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  #A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <-1  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q1 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q1 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q1.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q1.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q1.D50.4.PERM.ID_M", id), M.final.SA.LP.Q1.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q1.D50.4.LAND.ID_M", id), M.final.SA.LP.Q1.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}
#save.image("Paisaje_SA_LP_Q1_D50.4_EcoR_ID_345.RData")
#save.image("Paisaje_SA_LP_Q1_D50.4_EcoR_ID_345.RData")


###### FOR _ D50.4 _ Q5 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  #A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 5  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q5 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q5 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q5.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q5.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_M", id), M.final.SA.LP.Q5.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_M", id), M.final.SA.LP.Q5.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}

#save.image("Paisaje_SA_LP_Q5_D50.4_EcoR_ID_348_M.RData")

#save.image("Paisaje_SA_LP_TODO_D50.4_EcoR_ID_348_M.RData")

################################ PATAGONIA ##### 
ids <-348
A_list <- list()
###### FOR _ D50.4 _ Q0 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  #A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q0.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q0.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_P", id), M.final.SA.LP.Q0.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_P", id), M.final.SA.LP.Q0.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$LP[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q0.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q0.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q0.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q0.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q0.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q0.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}

#save.image("Paisaje_SA_LP_Q0_D50.4_EcoR_ID_348.RData")

###### FOR _ D50.4 _ Q025 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  #A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0.25  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q025 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q025 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas de M.LAND.D50.4 y M.Permanentes.D50.4
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q025.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q025.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q025.D50.4.PERM.ID_P", id), M.final.SA.LP.Q025.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q025.D50.4.LAND.ID_P", id), M.final.SA.LP.Q025.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$LP[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q025.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q025.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q025.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q025.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q025.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q025.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}

#save.image("Paisaje_SA_LP_Q025_D50.4_EcoR_ID_348.RData")

###### FOR _ D50.4 _ Q1 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  #A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <-1  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q1 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q1 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q1.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q1.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q1.D50.4.PERM.ID_P", id), M.final.SA.LP.Q1.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q1.D50.4.LAND.ID_P", id), M.final.SA.LP.Q1.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$LP[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q1.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q1.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q1.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q1.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q1.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q1.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}

#save.image("Paisaje_SA_LP_Q1_D50.4_EcoR_ID_348_P.RData")


###### FOR _ D50.4 _ Q5 ######
for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  #A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Visualiza parte de la matriz de distancias
  dim(aa_en_metros)
  aa_en_metros[1:10,1:10]
  
  # Ploteo inicial de los puntos de la EcoR en un mapa
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1.8)
  par(mfrow = c(1, 2))
  plot(M.SA.J.LP$CENTROID_Y ~ M.SA.J.LP$CENTROID_X, cex = 0.1, col = "light gray", pch = 15)
  points(A$CENTROID_Y ~ A$CENTROID_X, cex = 0.1, col = "orange", pch = 15)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  par(mar = c(4, 4, .5, .5), cex = 1.4, bty = "l")
  plot(exp(-b.mig * (0:200)) ~ I(0:200), ylim = c(0, 1), type = "l", lwd = 4, col = "dark green", 
       ylab = "Dispersal probability", xlab = "Distance (Km)")
  abline(h = 0.5, v = 4)
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 5  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Identifica áreas permanentes y terrestres en la EcoR
  id.land <- which(A$LAND != 0); length(id.land)
  id.permanentes <- which(A$FRESHWATER.permanentes != 0); length(id.permanentes)
  
  # Calculo de grados de la red (in y out) para las áreas permanentes y terrestres
  library(sna)
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "indegree") -> grado.in.permanente.D50.4
  degree(mm.permanente.D50.4[id.permanentes, id.permanentes], gmode = "digraph", cmode = "outdegree") -> grado.out.permanente.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "indegree") -> grado.in.land.D50.4
  degree(mm.LAND.D50.4[id.land, id.land], gmode = "digraph", cmode = "outdegree") -> grado.out.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  detach("package:sna", unload = TRUE)
  library(igraph) 
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  betweenness(g0.permanente.D50.4, directed = TRUE) -> bc.permanente.D50.4
  betweenness(g0.land.D50.4, directed = TRUE) -> bc.land.D50.4
  
  # Cálculo de betweenness centrality para las redes de áreas permanentes y terrestres
  graph_from_adjacency_matrix(1 / (mm.permanente.D50.4[id.permanentes, id.permanentes] + 1e-5), mode = "directed", weighted = TRUE) -> g0.permanente.D50.4
  graph_from_adjacency_matrix(1 / (mm.LAND.D50.4[id.land, id.land] + 1e-5), mode = "directed", weighted = TRUE) -> g0.land.D50.4
  
  # Calcular betweenness sin pasar el argumento directed
  bc.permanente.D50.4 <- betweenness(g0.permanente.D50.4)
  bc.land.D50.4 <- betweenness(g0.land.D50.4)
  
  # Crea matrices finales con los resultados
  M.Permanentes.D50.4 <- data.frame(A[id.permanentes, c(columnas_base_red, "J.freshwater")], `EF(1).T(2).P(3).` = rep(3, length(id.permanentes)), grado.in.D50.4 = grado.in.permanente.D50.4, grado.out.D50.4 = grado.out.permanente.D50.4, bc.D50.4 = bc.permanente.D50.4, check.names = FALSE)
  M.LAND.D50.4 <- data.frame(A[id.land, c(columnas_base_red, "J.LAND")], `EF(1).T(2).P(3).` = rep(100, length(id.land)), grado.in.D50.4 = grado.in.land.D50.4, grado.out.D50.4 = grado.out.land.D50.4, bc.D50.4 = bc.land.D50.4, check.names = FALSE)
  
  # Calcular la nueva columna J
  M.Permanentes.D50.4$J.q5 <- M.Permanentes.D50.4$J.freshwater * (1 - M.Permanentes.D50.4$LP)^q
  M.LAND.D50.4$J.q5 <-M.LAND.D50.4$J.LAND* (1 - M.LAND.D50.4$LP)^q
  
  # Renombra columnas
  names(M.LAND.D50.4)[names(M.LAND.D50.4) == "J.LAND"] <- "J.uk"
  names(M.Permanentes.D50.4)[names(M.Permanentes.D50.4) == "J.freshwater"] <- "J.uk"
  
  # Almacena matrices finales para la eco-región actual
  M.LAND <- A[, columnas_base_red]
  for (columna in setdiff(names(M.LAND.D50.4), names(M.LAND))) M.LAND[[columna]] <- 0
  M.LAND <- M.LAND[, names(M.LAND.D50.4)]
  
  M.final.SA.LP.Q5.D50.4.PERM <- rbind(M.Permanentes.D50.4, M.LAND)
  M.final.SA.LP.Q5.D50.4.LAND <- rbind(M.LAND.D50.4, M.LAND)
  
  
  # Asigna los resultados a variables dinámicas para cada eco-región
  assign(paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_P", id), M.final.SA.LP.Q5.D50.4.PERM)
  assign(paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_P", id), M.final.SA.LP.Q5.D50.4.LAND)
  
  # Genera gráficos para cada métrica calculada (LP, In-degree, Out-degree, Betweenness)
  par(mar = c(2.5, 2.5, 1.5, 1.5), bty = "l", cex = 1)
  par(mfrow = c(1, 7), oma = c(0, 0, 0, 0))
  
  #### Gráfico de contaminación lumínica (LP)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y ~ M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X, col = "gray", cex = 0.4, pch = 15, main = "Light Pollution")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$LP[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)] ~ 
           M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)` == 100)], cex = 0.4,
         col = colorRampPalette(c("#212C52", "#6B74F7", "#DE247B","#FF7900","#FFDE94"))(599)[ii], pch = 15)
  
  #### In_Degree_MassEffect #### D50=4Km
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="In-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$grado.in.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  #### Out_Degree_SourcePatchEffect ####
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Out-Degree Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$grado.out.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  
  #### Betweenness ####
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness LAND (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), 
                         max(M.final.SA.LP.Q5.D50.4.LAND$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.LAND$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)]~M.final.SA.LP.Q5.D50.4.LAND$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.LAND$`EF(1).T(2).P(3)`==100)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 0.4)
  
  rm(ii)
  plot(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Betweenness Permanent (4Km)")
  ii <- cut(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)], 
            breaks = seq(min(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), 
                         max(M.final.SA.LP.Q5.D50.4.PERM$bc.D50.4[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
  points(M.final.SA.LP.Q5.D50.4.PERM$CENTROID_Y[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)]~M.final.SA.LP.Q5.D50.4.PERM$CENTROID_X[which(M.final.SA.LP.Q5.D50.4.PERM$`EF(1).T(2).P(3)`==3)],cex=0.4,
         col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
  
  Sys.sleep(1) 
}
















################################# D50.4 - SUMAR 348 al resto #####
if (!exists("M.final.SA.LP.Q0.D50.4.LAND.ID_301")) load("Paisaje_SA_LP_Q0_D50.4_EcoR_ID_301_352.RData")
if (!exists("M.final.SA.LP.Q025.D50.4.LAND.ID_301")) load("Paisaje_SA_LP_Q025_D50.4_EcoR_ID_301_352.RData")
if (!exists("M.final.SA.LP.Q1.D50.4.LAND.ID_301")) load("Paisaje_SA_LP_Q1_D50.4_EcoR_ID_301_352.RData")
if (!exists("M.final.SA.LP.Q5.D50.4.LAND.ID_301")) load("Paisaje_SA_LP_Q5_D50.4_EcoR_ID_301_352.RData")
if (!exists("M.final.SA.LP.Q5.D50.4.LAND.ID_M348") || !exists("M.final.SA.LP.Q5.D50.4.LAND.ID_P348")) load("Paisaje_SA_LP_Q_TODO_D50.4_EcoR_ID_348_PM.RData")
######## renombrar para plotear Malvinas y Patagonia #####
# Patagonia queda con id 348
# Malvinas le asigno id 353 
M.final.SA.LP.Q0.D50.4.LAND.ID_353 <- M.final.SA.LP.Q0.D50.4.LAND.ID_M348
M.final.SA.LP.Q025.D50.4.LAND.ID_353 <- M.final.SA.LP.Q025.D50.4.LAND.ID_M348
M.final.SA.LP.Q1.D50.4.LAND.ID_353 <- M.final.SA.LP.Q1.D50.4.LAND.ID_M348
M.final.SA.LP.Q5.D50.4.LAND.ID_353 <- M.final.SA.LP.Q5.D50.4.LAND.ID_M348

M.final.SA.LP.Q0.D50.4.LAND.ID_348 <- M.final.SA.LP.Q0.D50.4.LAND.ID_P348
M.final.SA.LP.Q025.D50.4.LAND.ID_348 <- M.final.SA.LP.Q025.D50.4.LAND.ID_P348
M.final.SA.LP.Q1.D50.4.LAND.ID_348 <- M.final.SA.LP.Q1.D50.4.LAND.ID_P348
M.final.SA.LP.Q5.D50.4.LAND.ID_348 <- M.final.SA.LP.Q5.D50.4.LAND.ID_P348

M.final.SA.LP.Q0.D50.4.PERM.ID_353 <- M.final.SA.LP.Q0.D50.4.PERM.ID_M348
M.final.SA.LP.Q025.D50.4.PERM.ID_353 <- M.final.SA.LP.Q025.D50.4.PERM.ID_M348
M.final.SA.LP.Q1.D50.4.PERM.ID_353 <- M.final.SA.LP.Q1.D50.4.PERM.ID_M348
M.final.SA.LP.Q5.D50.4.PERM.ID_353 <- M.final.SA.LP.Q5.D50.4.PERM.ID_M348

M.final.SA.LP.Q0.D50.4.PERM.ID_348 <- M.final.SA.LP.Q0.D50.4.PERM.ID_P348
M.final.SA.LP.Q025.D50.4.PERM.ID_348 <- M.final.SA.LP.Q025.D50.4.PERM.ID_P348
M.final.SA.LP.Q1.D50.4.PERM.ID_348 <- M.final.SA.LP.Q1.D50.4.PERM.ID_P348
M.final.SA.LP.Q5.D50.4.PERM.ID_348 <- M.final.SA.LP.Q5.D50.4.PERM.ID_P348


save(list = ls(pattern = "^M\\.final\\.SA\\.LP\\."), file = "Resultados_redes_D50.4.RData")