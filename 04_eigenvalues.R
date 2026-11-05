################################################################################
## Reichmann et al.
## Artificial light at night reorganizes habitat connectivity 
## and erodes landscape capacity to sustain biodiversity across South America
################################################################################
# This script is part of the code accompanying the manuscript above.
#
# 04 — Eigenvalues
# Calculates landscape carrying capacity by ecoregion using eigenvalues and
# adds the final table to the results file.


# Carga de librerías
library(geosphere)  # Para cálculo de distancias
library(sp)         # Manejo de datos espaciales
library(stars)      # Manipulación de datos geométricos
library(sna)        # Para medidas de centralidad en redes
library(igraph)
library(RSpectra)
library(sf)
library(readr) 

setwd(".")
load("Paisaje_SA_EcoR_M_SA_J_LP.RData")
load("Resultados_redes_D50.4.RData")


# Visualiación de estructura inicial de la matriz de datos
head(M.SA.J.LP) 
table(M.SA.J.LP$EcoR_ID)
M.SA.J.LP$LP <- M.SA.J.LP$LPF_mean_est # LP normalizado entre 0 y 1, con maximo en 5 nw
range(M.SA.J.LP$LP, na.rm = TRUE)

unique(M.SA.J.LP$EcoR_Name)
unique(M.SA.J.LP$EcoR_ID)


################################## EIGENVALUES #####
ids <-301:352
ids <- setdiff(ids, c(348, 351)) # 348 va por separado en Malvinas y Patagonia - 351 Juan Fernandez 1 celda
A_list <- list()

###### FOR _ D50.4 _ Q0 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q0 <- data.frame(EcoR_ID = character(), 
                             Lambda_1_LAND = numeric(), 
                             Lambda_1_PERM = numeric(), 
                             stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q0 <- rbind(eigenvalues_q0, 
                          data.frame(EcoR_ID = as.character(id), 
                                     Lambda_1_LAND = lambda_1_land, 
                                     Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
eigenvalues_q0$Lambda_1_LAND <- Re(eigenvalues_q0$Lambda_1_LAND)
eigenvalues_q0$Lambda_1_PERM <- Re(eigenvalues_q0$Lambda_1_PERM)
print(eigenvalues_q0)
###### FOR _ D50.4 _ Q025 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q025 <- data.frame(EcoR_ID = character(), 
                             Lambda_1_LAND = numeric(), 
                             Lambda_1_PERM = numeric(), 
                             stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0.25  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q025 <- rbind(eigenvalues_q025, 
                          data.frame(EcoR_ID = as.character(id), 
                                     Lambda_1_LAND = lambda_1_land, 
                                     Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
eigenvalues_q025$Lambda_1_LAND <- Re(eigenvalues_q025$Lambda_1_LAND)
eigenvalues_q025$Lambda_1_PERM <- Re(eigenvalues_q025$Lambda_1_PERM)
print(eigenvalues_q025)

###### FOR _ D50.4 _ Q5 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q5 <- data.frame(EcoR_ID = character(), 
                             Lambda_1_LAND = numeric(), 
                             Lambda_1_PERM = numeric(), 
                             stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 5  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q5 <- rbind(eigenvalues_q5, 
                          data.frame(EcoR_ID = as.character(id), 
                                     Lambda_1_LAND = lambda_1_land, 
                                     Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
print(eigenvalues_q5)
eigenvalues_q5$Lambda_1_LAND <- Re(eigenvalues_q5$Lambda_1_LAND)
eigenvalues_q5$Lambda_1_PERM <- Re(eigenvalues_q5$Lambda_1_PERM)

print(eigenvalues_q5)
print(eigenvalues_q0)









## MALVINAS ##### 
ids <-348
A_list <- list()
###### FOR _ D50.4 _ Q0 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q0_M <- data.frame(EcoR_ID = character(), 
                             Lambda_1_LAND = numeric(), 
                             Lambda_1_PERM = numeric(), 
                             stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  #A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q0_M <- rbind(eigenvalues_q0_M, 
                          data.frame(EcoR_ID = as.character(id), 
                                     Lambda_1_LAND = lambda_1_land, 
                                     Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
eigenvalues_q0_M$Lambda_1_LAND <- Re(eigenvalues_q0_M$Lambda_1_LAND)
eigenvalues_q0_M$Lambda_1_PERM <- Re(eigenvalues_q0_M$Lambda_1_PERM)
print(eigenvalues_q0_M)
###### FOR _ D50.4 _ Q025 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q025_M <- data.frame(EcoR_ID = character(), 
                               Lambda_1_LAND = numeric(), 
                               Lambda_1_PERM = numeric(), 
                               stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  #A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q025_M <- rbind(eigenvalues_q025_M, 
                            data.frame(EcoR_ID = as.character(id), 
                                       Lambda_1_LAND = lambda_1_land, 
                                       Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
eigenvalues_q025_M$Lambda_1_LAND <- Re(eigenvalues_q025_M$Lambda_1_LAND)
eigenvalues_q025_M$Lambda_1_PERM <- Re(eigenvalues_q025_M$Lambda_1_PERM)
print(eigenvalues_q025_M)
###### FOR _ D50.4 _ Q5 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q5_M <- data.frame(EcoR_ID = character(), 
                             Lambda_1_LAND = numeric(), 
                             Lambda_1_PERM = numeric(), 
                             stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  #A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 5  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q5_M <- rbind(eigenvalues_q5_M, 
                          data.frame(EcoR_ID = as.character(id), 
                                     Lambda_1_LAND = lambda_1_land, 
                                     Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
print(eigenvalues_q5_M)
eigenvalues_q5_M$Lambda_1_LAND <- Re(eigenvalues_q5_M$Lambda_1_LAND)
eigenvalues_q5_M$Lambda_1_PERM <- Re(eigenvalues_q5_M$Lambda_1_PERM)

print(eigenvalues_q5_M)
print(eigenvalues_q0_M)

## PATAGONIA ##### 
ids <-348
A_list <- list()
###### FOR _ D50.4 _ Q0 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q0_P <- data.frame(EcoR_ID = character(), 
                             Lambda_1_LAND = numeric(), 
                             Lambda_1_PERM = numeric(), 
                             stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  #A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q0_P <- rbind(eigenvalues_q0_P, 
                          data.frame(EcoR_ID = as.character(id), 
                                     Lambda_1_LAND = lambda_1_land, 
                                     Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
eigenvalues_q0_P$Lambda_1_LAND <- Re(eigenvalues_q0_P$Lambda_1_LAND)
eigenvalues_q0_P$Lambda_1_PERM <- Re(eigenvalues_q0_P$Lambda_1_PERM)
print(eigenvalues_q0_P) 
###### FOR _ D50.4 _ Q025 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q025_P <- data.frame(EcoR_ID = character(), 
                               Lambda_1_LAND = numeric(), 
                               Lambda_1_PERM = numeric(), 
                               stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  #A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 0  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q025_P <- rbind(eigenvalues_q025_P, 
                            data.frame(EcoR_ID = as.character(id), 
                                       Lambda_1_LAND = lambda_1_land, 
                                       Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
eigenvalues_q025_P$Lambda_1_LAND <- Re(eigenvalues_q025_P$Lambda_1_LAND)
eigenvalues_q025_P$Lambda_1_PERM <- Re(eigenvalues_q025_P$Lambda_1_PERM)
print(eigenvalues_q025_P) 
###### FOR _ D50.4 _ Q5 ######
# Crear una tabla vacía para almacenar los eigenvalues
eigenvalues_q5_P <- data.frame(EcoR_ID = character(), 
                             Lambda_1_LAND = numeric(), 
                             Lambda_1_PERM = numeric(), 
                             stringsAsFactors = FALSE)

for (id in ids) {
  # Filtra la matriz para cada EcoR
  A <- M.SA.J.LP[which(M.SA.J.LP$EcoR_ID == as.character(id)), ]
  A_list[[as.character(id)]] <- A  # Guarda los datos en la lista
  #A<-A[A$CENTROID_X > -62, ] ### PARA 348 MALVINAS
  A<-A[A$CENTROID_X < -62, ] ### PARA 348 PATAGONIA
  # Si no hay datos, continuar con la siguiente iteración
  if (nrow(A) == 0) next  
  
  # Transforma las coordenadas de la EcoR a formato espacial
  coord.xy <- cbind(A$CENTROID_X, A$CENTROID_Y)
  colnames(coord.xy) <- c("CENTROID_X", "CENTROID_Y")
  LAND_coord <- st_as_sf(as.data.frame(coord.xy), coords = c("CENTROID_X", "CENTROID_Y"))  # A formato 'sf'
  LAND_coord_spatial <- as(LAND_coord, "Spatial")  # Convertir a formato 'sp' para cálculos de distancia
  
  # Cálculo de la matriz de distancias en metros
  aa_en_metros <- distm(LAND_coord_spatial@coords, LAND_coord_spatial@coords, fun = distGeo)
  
  # Cálculo de probabilidad de migración basada en la distancia D50 (4 km)
  D50 <- 4  
  b.mig <- (-log(.5) / D50)  # Parámetro de migración basado en la distancia
  
  b.mig_D50.4 <- b.mig
  
  # Define parámetro q que modula el impacto de la contaminación lumínica
  q <- 5  # Modifica la atenuación de J por LP
  
  # Cálculo de la probabilidad de migración
  mm.permanente.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.permanente * (1 - A$LP)^q
  mm.LAND.D50.4 <- exp(-b.mig_D50.4 * (aa_en_metros / 1000)) * A$J.LAND * (1 - A$LP)^q
  
  # Convertir la matriz en un grafo dirigido con pesos
  g_land <- graph_from_adjacency_matrix(mm.LAND.D50.4, mode = "directed", weighted = TRUE)
  g_perm <- graph_from_adjacency_matrix(mm.permanente.D50.4, mode = "directed", weighted = TRUE)
  
  # Extraer la matriz de adyacencia como un objeto numérico
  adj_land <- as.matrix(as_adjacency_matrix(g_land, attr = "weight"))
  adj_perm <- as.matrix(as_adjacency_matrix(g_perm, attr = "weight"))
  
  # Calcular el mayor eigenvalue (λ1) usando RSpectra (evita calcular todos los valores propios)
  lambda_1_land <- eigs(adj_land, k = 1, which = "LM")$values
  lambda_1_perm <- eigs(adj_perm, k = 1, which = "LM")$values
  
  # Almacenar los resultados en la tabla
  eigenvalues_q5_P <- rbind(eigenvalues_q5_P, 
                          data.frame(EcoR_ID = as.character(id), 
                                     Lambda_1_LAND = lambda_1_land, 
                                     Lambda_1_PERM = lambda_1_perm))
}

# Mostrar la tabla con los eigenvalues calculados para cada EcoR_ID
print(eigenvalues_q5_P)
eigenvalues_q5_P$Lambda_1_LAND <- Re(eigenvalues_q5_P$Lambda_1_LAND)
eigenvalues_q5_P$Lambda_1_PERM <- Re(eigenvalues_q5_P$Lambda_1_PERM)




print(eigenvalues_q5_P)
print(eigenvalues_q0_P)
print(eigenvalues_q5_M)
print(eigenvalues_q0_M)
write.csv(eigenvalues_q0_P, "eigenvalues_q0_P.csv", row.names = FALSE)
write.csv(eigenvalues_q0_M, "eigenvalues_q0_M.csv", row.names = FALSE)
write.csv(eigenvalues_q5_P, "eigenvalues_q5_P.csv", row.names = FALSE)
write.csv(eigenvalues_q5_M, "eigenvalues_q5_M.csv", row.names = FALSE)
write.csv(eigenvalues_q025_P, "eigenvalues_q025_P.csv", row.names = FALSE)
write.csv(eigenvalues_q025_M, "eigenvalues_q025_M.csv", row.names = FALSE)
write.csv(eigenvalues_q0, "eigenvalues_q0.csv", row.names = FALSE)
write.csv(eigenvalues_q025, "eigenvalues_q025.csv", row.names = FALSE)
write.csv(eigenvalues_q5, "eigenvalues_q5.csv", row.names = FALSE)


################################## TABLA UNICA EIGENVALUES ######
eigenvalues_q0 <- read_csv("eigenvalues_q0.csv")
eigenvalues_q025 <- read_csv("eigenvalues_q025.csv")
eigenvalues_q5 <- read_csv("eigenvalues_q5.csv")
eigenvalues_q0_M <- read_csv("eigenvalues_q0_M.csv")
eigenvalues_q0_P <- read_csv("eigenvalues_q0_P.csv")
eigenvalues_q025_M <- read_csv("eigenvalues_q025_M.csv")
eigenvalues_q025_P <- read_csv("eigenvalues_q025_P.csv")
eigenvalues_q5_M <- read_csv("eigenvalues_q5_M.csv")
eigenvalues_q5_P <- read_csv("eigenvalues_q5_P.csv")

print(eigenvalues_q5)
print(eigenvalues_q025)
print(eigenvalues_q0)
print(eigenvalues_q5_P)
print(eigenvalues_q025_P)
print(eigenvalues_q0_P)
print(eigenvalues_q5_M)
print(eigenvalues_q025_M)
print(eigenvalues_q0_M)

eigenvalues_q0_M$EcoR_ID <- 353 # sumo 353 para Malvinas, que separe de Patagonia, es originalmente 348 tb
eigenvalues_q025_M$EcoR_ID <- 353
eigenvalues_q5_M$EcoR_ID <- 353

eigenvalues_q0 <- rbind(eigenvalues_q0, eigenvalues_q0_P, eigenvalues_q0_M)
eigenvalues_q025 <- rbind(eigenvalues_q025, eigenvalues_q025_P, eigenvalues_q025_M)
eigenvalues_q5 <- rbind(eigenvalues_q5, eigenvalues_q5_P, eigenvalues_q5_M)
colnames(eigenvalues_q0)[colnames(eigenvalues_q0) == "Lambda_1_LAND"] <- "Lambda_1_LAND_q0"
colnames(eigenvalues_q0)[colnames(eigenvalues_q0) == "Lambda_1_PERM"] <- "Lambda_1_PERM_q0"
colnames(eigenvalues_q025)[colnames(eigenvalues_q025) == "Lambda_1_LAND"] <- "Lambda_1_LAND_q025"
colnames(eigenvalues_q025)[colnames(eigenvalues_q025) == "Lambda_1_PERM"] <- "Lambda_1_PERM_q025"
colnames(eigenvalues_q5)[colnames(eigenvalues_q5) == "Lambda_1_LAND"] <- "Lambda_1_LAND_q5"
colnames(eigenvalues_q5)[colnames(eigenvalues_q5) == "Lambda_1_PERM"] <- "Lambda_1_PERM_q5"
eigenvalues <- merge(eigenvalues_q0, eigenvalues_q025, by = "EcoR_ID", all = TRUE)
eigenvalues <- merge(eigenvalues, eigenvalues_q5, by = "EcoR_ID", all = TRUE)





# Crear la columna "dif" con la diferencia absoluta
eigenvalues$dif_LAND_q025 <-  eigenvalues$Lambda_1_LAND_q025 - eigenvalues$Lambda_1_LAND_q0
eigenvalues$dif_PERM_q025 <-  eigenvalues$Lambda_1_PERM_q025 - eigenvalues$Lambda_1_PERM_q0
eigenvalues$dif_LAND_q5 <- eigenvalues$Lambda_1_LAND_q5 - eigenvalues$Lambda_1_LAND_q0
eigenvalues$dif_PERM_q5 <-  eigenvalues$Lambda_1_PERM_q5 - eigenvalues$Lambda_1_PERM_q0
# Crear la columna "dif_perc" con la diferencia porcentual
eigenvalues$dif_perc_LAND_q025 <- (eigenvalues$dif_LAND_q025 / eigenvalues$Lambda_1_LAND_q0) * 100
eigenvalues$dif_perc_PERM_q025 <- (eigenvalues$dif_PERM_q025 / eigenvalues$Lambda_1_PERM_q0) * 100
eigenvalues$dif_perc_LAND_q5 <- (eigenvalues$dif_LAND_q5 / eigenvalues$Lambda_1_LAND_q0) * 100
eigenvalues$dif_perc_PERM_q5 <- (eigenvalues$dif_PERM_q5 / eigenvalues$Lambda_1_PERM_q0) * 100
# Incorpora la tabla final al archivo de resultados, sin alterar las redes.
local({
  resultados <- new.env()
  load("Resultados_redes_D50.4.RData", envir = resultados)
  resultados$eigenvalues <- eigenvalues
  temporal <- tempfile(pattern = "resultados_", tmpdir = ".", fileext = ".RData")
  save(list = ls(resultados), envir = resultados, file = temporal, compress = TRUE)
  stopifnot(file.rename(temporal, "Resultados_redes_D50.4.RData"))
})

################################## BC resumen ##########
load("Resultados_redes_D50.4.RData")
### FUNCIONES ######
ratio.in <- function(entrada, entrada_Q0, entrada_Q1, entrada_Q025, entrada_Q5) {
  out = NULL
  for (i in "grado.in.D50.4") {
    unimpacted <- entrada[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID", i)]
    Q0 <- entrada_Q0[, i]
    Q1 <- entrada_Q1[, i]
    Q025 <- entrada_Q025[, i]
    Q5 <- entrada_Q5[, i]
    P1 <- 100 - 100 * (Q1 / Q0)
    P025 <- 100 - 100 * (Q025 / Q0)
    P5 <- 100 - 100 * (Q5 / Q0)
    salida <- cbind(unimpacted[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID")], as.vector(Q0/unimpacted[[i]]), as.vector(Q1/unimpacted[[i]]), 
                    as.vector(Q025/unimpacted[[i]]), as.vector(Q5/unimpacted[[i]]), as.vector(P1), as.vector(P025), as.vector(P5))
    colnames(salida) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
  }
  out <- salida
}

ratio.out <- function(entrada, entrada_Q0, entrada_Q5, entrada_Q1, entrada_Q025) {
  out = NULL
  for (i in "grado.out.D50.4") {
    unimpacted <- entrada[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID", i)]
    Q0 <- entrada_Q0[, i]
    Q1 <- entrada_Q1[, i]
    Q025 <- entrada_Q025[, i]
    Q5 <- entrada_Q5[, i]
    P1 <- 100 - 100 * (Q1 / Q0)
    P025 <- 100 - 100 * (Q025 / Q0)
    P5 <- 100 - 100 * (Q5 / Q0)
    salida <- cbind(unimpacted[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID")], as.vector(Q0/unimpacted[[i]]), as.vector(Q1/unimpacted[[i]]), 
                    as.vector(Q025/unimpacted[[i]]), as.vector(Q5/unimpacted[[i]]), as.vector(P1), as.vector(P025), as.vector(P5))
    colnames(salida) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
  }
  out <- salida
}

ratio.bc <- function(entrada, entrada_Q0, entrada_Q5, entrada_Q1, entrada_Q025) {
  out = NULL
  for (i in "bc.D50.4") {
    unimpacted <- entrada[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID", i)]
    Q0 <- entrada_Q0[, i]
    Q1 <- entrada_Q1[, i]
    Q025 <- entrada_Q025[, i]
    Q5 <- entrada_Q5[, i]
    P0 <- 100 - 100 * (Q0 / Q0)
    P1 <- 100 - 100 * (Q1 / Q0)
    P025 <- 100 - 100 * (Q025 / Q0)
    P5 <- 100 - 100 * (Q5 / Q0)
    salida <- cbind(unimpacted[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID")], as.vector(Q0/unimpacted[[i]]), as.vector(Q1/unimpacted[[i]]), 
                    as.vector(Q025/unimpacted[[i]]), as.vector(Q5/unimpacted[[i]]), as.vector(P0), as.vector(P1), as.vector(P025), as.vector(P5))
    colnames(salida) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P0", "P1", "P025", "P5")
  }
  out <- salida
}

ratiolog.bc <- function(entrada, entrada_Q0, entrada_Q5, entrada_Q1, entrada_Q025) {
  out = NULL
  for (i in "bc.D50.4") {
    unimpacted <- entrada[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID", i)]
    Q0 <- entrada_Q0[, i]
    Q1 <- entrada_Q1[, i]
    Q025 <- entrada_Q025[, i]
    Q5 <- entrada_Q5[, i]
    log_Q0 <- log(Q0 / Q0)
    log_Q1 <- log(Q1 / Q0)
    log_Q025 <- log(Q025 / Q0)
    log_Q5 <- log(Q5 / Q0)
    log_P0 <- 100 - 100 * (Q0 / Q0)
    log_P1 <- (Q1 / Q0) / (max(Q1 / Q0))
    log_P025 <- 100 - 100 * (Q025 / Q0)
    log_P5 <- 100 - 100 * (Q5 / Q0)
    salida <- cbind(unimpacted[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID")], as.vector(log_Q0), as.vector(log_Q1), as.vector(log_Q025), 
                    as.vector(log_Q5), as.vector(log_P0), as.vector(log_P1), as.vector(log_P025), as.vector(log_P5))
    colnames(salida) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P0", "P1", "P025", "P5")
  }
  out <- salida
}

### Ratios para LAND ######
id <- 301:353 # sumo 353 para Malvinas, que separe de Patagonia, es originalmente 348 tb
id <- setdiff(id, 351) # 351 Juan Fernandez 1 sola celda
Resultados.in.land <- list()
for (i in id) {
  matrix_name_entrada <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i); entrada <- get(matrix_name_entrada)
  matrix_name_entrada_Q0 <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i); entrada_Q0 <- get(matrix_name_entrada_Q0)
  matrix_name_entrada_Q025 <- paste0("M.final.SA.LP.Q025.D50.4.LAND.ID_", i); entrada_Q025 <- get(matrix_name_entrada_Q025)
  matrix_name_entrada_Q1 <- paste0("M.final.SA.LP.Q1.D50.4.LAND.ID_", i); entrada_Q1 <- get(matrix_name_entrada_Q1)
  matrix_name_entrada_Q5 <- paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_", i); entrada_Q5 <- get(matrix_name_entrada_Q5)
  
  # Ejecutar y almacenar el resultado
  Resultados.in.land[[i]] <- ratio.in(entrada = entrada, entrada_Q0 = entrada_Q0, entrada_Q025 = entrada_Q025, entrada_Q1 = entrada_Q1, entrada_Q5 = entrada_Q5)
  assign(paste0("ratio.in.D50.4.LAND.ID_", i), Resultados.in.land[[i]])
}

Resultados.out.land <- list()
for (i in id) {
  matrix_name_entrada <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i); entrada <- get(matrix_name_entrada)
  matrix_name_entrada_Q0 <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i); entrada_Q0 <- get(matrix_name_entrada_Q0)
  matrix_name_entrada_Q025 <- paste0("M.final.SA.LP.Q025.D50.4.LAND.ID_", i); entrada_Q025 <- get(matrix_name_entrada_Q025)
  matrix_name_entrada_Q1 <- paste0("M.final.SA.LP.Q1.D50.4.LAND.ID_", i); entrada_Q1 <- get(matrix_name_entrada_Q1)
  matrix_name_entrada_Q5 <- paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_", i); entrada_Q5 <- get(matrix_name_entrada_Q5)
  
  # Ejecutar y almacenar el resultado
  Resultados.out.land[[i]] <- ratio.out(entrada = entrada, entrada_Q0 = entrada_Q0, entrada_Q025 = entrada_Q025, entrada_Q1 = entrada_Q1, entrada_Q5 = entrada_Q5)
  assign(paste0("ratio.out.D50.4.LAND.ID_", i), Resultados.out.land[[i]])
}

Resultados.bc.land <- list()
for (i in id) {
  matrix_name_entrada <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i); entrada <- get(matrix_name_entrada)
  matrix_name_entrada_Q0 <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i); entrada_Q0 <- get(matrix_name_entrada_Q0)
  matrix_name_entrada_Q025 <- paste0("M.final.SA.LP.Q025.D50.4.LAND.ID_", i); entrada_Q025 <- get(matrix_name_entrada_Q025)
  matrix_name_entrada_Q1 <- paste0("M.final.SA.LP.Q1.D50.4.LAND.ID_", i); entrada_Q1 <- get(matrix_name_entrada_Q1)
  matrix_name_entrada_Q5 <- paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_", i); entrada_Q5 <- get(matrix_name_entrada_Q5)
  
  # Ejecutar y almacenar el resultado
  Resultados.bc.land[[i]] <- ratio.bc(entrada = entrada, entrada_Q0 = entrada_Q0, entrada_Q025 = entrada_Q025, entrada_Q1 = entrada_Q1, entrada_Q5 = entrada_Q5)
  assign(paste0("ratio.bc.D50.4.LAND.ID_", i), Resultados.bc.land[[i]])
}

Resultados.logbc.land <- list()
for (i in id) {
  matrix_name_entrada <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i); entrada <- get(matrix_name_entrada)
  matrix_name_entrada_Q0 <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i); entrada_Q0 <- get(matrix_name_entrada_Q0)
  matrix_name_entrada_Q025 <- paste0("M.final.SA.LP.Q025.D50.4.LAND.ID_", i); entrada_Q025 <- get(matrix_name_entrada_Q025)
  matrix_name_entrada_Q1 <- paste0("M.final.SA.LP.Q1.D50.4.LAND.ID_", i); entrada_Q1 <- get(matrix_name_entrada_Q1)
  matrix_name_entrada_Q5 <- paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_", i); entrada_Q5 <- get(matrix_name_entrada_Q5)
  #"corregir" el bc de 0, para el ratiolog
  entrada[["bc.D50.4"]]<- ifelse(entrada[["bc.D50.4"]]==0,0.00005, entrada[["bc.D50.4"]])
  entrada_Q0[["bc.D50.4"]]<- ifelse(entrada_Q0[["bc.D50.4"]]==0,0.00005, entrada_Q0[["bc.D50.4"]])
  entrada_Q025[["bc.D50.4"]]<- ifelse(entrada_Q025[["bc.D50.4"]]==0,0.00005, entrada_Q025[["bc.D50.4"]])
  entrada_Q1[["bc.D50.4"]]<- ifelse(entrada_Q1[["bc.D50.4"]]==0,0.00005, entrada_Q1[["bc.D50.4"]])
  entrada_Q5[["bc.D50.4"]]<- ifelse(entrada_Q5[["bc.D50.4"]]==0,0.00005, entrada_Q5[["bc.D50.4"]])
  # Ejecutar y almacenar el resultado
  Resultados.logbc.land[[i]] <- ratiolog.bc(entrada = entrada, entrada_Q0 = entrada_Q0, entrada_Q025 = entrada_Q025, entrada_Q1 = entrada_Q1, entrada_Q5 = entrada_Q5)
  assign(paste0("ratiolog.bc.D50.4.LAND.ID_", i), Resultados.logbc.land[[i]])
}
### Ratios para PERM ######
id <- 301:353 # sumo 353 para Malvinas, que separe de Patagonia, es originalmente 348 tb
id <- setdiff(id, 351) # 351 Juan Fernandez 1 sola celda
Resultados.in.perm <- list()
for (i in id) {
  matrix_name_entrada <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i); entrada <- get(matrix_name_entrada)
  matrix_name_entrada_Q0 <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i); entrada_Q0 <- get(matrix_name_entrada_Q0)
  matrix_name_entrada_Q025 <- paste0("M.final.SA.LP.Q025.D50.4.PERM.ID_", i); entrada_Q025 <- get(matrix_name_entrada_Q025)
  matrix_name_entrada_Q1 <- paste0("M.final.SA.LP.Q1.D50.4.PERM.ID_", i); entrada_Q1 <- get(matrix_name_entrada_Q1)
  matrix_name_entrada_Q5 <- paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_", i); entrada_Q5 <- get(matrix_name_entrada_Q5)
  
  # Ejecutar y almacenar el resultado
  Resultados.in.perm[[i]] <- ratio.in(entrada = entrada, entrada_Q0 = entrada_Q0, entrada_Q025 = entrada_Q025, entrada_Q1 = entrada_Q1, entrada_Q5 = entrada_Q5)
  assign(paste0("ratio.in.D50.4.PERM.ID_", i), Resultados.in.perm[[i]])
}

Resultados.out.perm <- list()
for (i in id) {
  matrix_name_entrada <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i); entrada <- get(matrix_name_entrada)
  matrix_name_entrada_Q0 <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i); entrada_Q0 <- get(matrix_name_entrada_Q0)
  matrix_name_entrada_Q025 <- paste0("M.final.SA.LP.Q025.D50.4.PERM.ID_", i); entrada_Q025 <- get(matrix_name_entrada_Q025)
  matrix_name_entrada_Q1 <- paste0("M.final.SA.LP.Q1.D50.4.PERM.ID_", i); entrada_Q1 <- get(matrix_name_entrada_Q1)
  matrix_name_entrada_Q5 <- paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_", i); entrada_Q5 <- get(matrix_name_entrada_Q5)
  
  # Ejecutar y almacenar el resultado
  Resultados.out.perm[[i]] <- ratio.out(entrada = entrada, entrada_Q0 = entrada_Q0, entrada_Q025 = entrada_Q025, entrada_Q1 = entrada_Q1, entrada_Q5 = entrada_Q5)
  assign(paste0("ratio.out.D50.4.PERM.ID_", i), Resultados.out.perm[[i]])
}

Resultados.bc.perm <- list()
for (i in id) {
  matrix_name_entrada <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i); entrada <- get(matrix_name_entrada)
  matrix_name_entrada_Q0 <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i); entrada_Q0 <- get(matrix_name_entrada_Q0)
  matrix_name_entrada_Q025 <- paste0("M.final.SA.LP.Q025.D50.4.PERM.ID_", i); entrada_Q025 <- get(matrix_name_entrada_Q025)
  matrix_name_entrada_Q1 <- paste0("M.final.SA.LP.Q1.D50.4.PERM.ID_", i); entrada_Q1 <- get(matrix_name_entrada_Q1)
  matrix_name_entrada_Q5 <- paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_", i); entrada_Q5 <- get(matrix_name_entrada_Q5)
  
  # Ejecutar y almacenar el resultado
  Resultados.bc.perm[[i]] <- ratio.bc(entrada = entrada, entrada_Q0 = entrada_Q0, entrada_Q025 = entrada_Q025, entrada_Q1 = entrada_Q1, entrada_Q5 = entrada_Q5)
  assign(paste0("ratio.bc.D50.4.PERM.ID_", i), Resultados.bc.perm[[i]])
}

Resultados.logbc.perm <- list()
for (i in id) {
  matrix_name_entrada <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i); entrada <- get(matrix_name_entrada)
  matrix_name_entrada_Q0 <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i); entrada_Q0 <- get(matrix_name_entrada_Q0)
  matrix_name_entrada_Q025 <- paste0("M.final.SA.LP.Q025.D50.4.PERM.ID_", i); entrada_Q025 <- get(matrix_name_entrada_Q025)
  matrix_name_entrada_Q1 <- paste0("M.final.SA.LP.Q1.D50.4.PERM.ID_", i); entrada_Q1 <- get(matrix_name_entrada_Q1)
  matrix_name_entrada_Q5 <- paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_", i); entrada_Q5 <- get(matrix_name_entrada_Q5)
  #"corregir" el bc de 0, para el ratiolog
  entrada[["bc.D50.4"]]<- ifelse(entrada[["bc.D50.4"]]==0,0.00005, entrada[["bc.D50.4"]])
  entrada_Q0[["bc.D50.4"]]<- ifelse(entrada_Q0[["bc.D50.4"]]==0,0.00005, entrada_Q0[["bc.D50.4"]])
  entrada_Q025[["bc.D50.4"]]<- ifelse(entrada_Q025[["bc.D50.4"]]==0,0.00005, entrada_Q025[["bc.D50.4"]])
  entrada_Q1[["bc.D50.4"]]<- ifelse(entrada_Q1[["bc.D50.4"]]==0,0.00005, entrada_Q1[["bc.D50.4"]])
  entrada_Q5[["bc.D50.4"]]<- ifelse(entrada_Q5[["bc.D50.4"]]==0,0.00005, entrada_Q5[["bc.D50.4"]])
  # Ejecutar y almacenar el resultado
  Resultados.logbc.perm[[i]] <- ratiolog.bc(entrada = entrada, entrada_Q0 = entrada_Q0, entrada_Q025 = entrada_Q025, entrada_Q1 = entrada_Q1, entrada_Q5 = entrada_Q5)
  assign(paste0("ratiolog.bc.D50.4.PERM.ID_", i), Resultados.logbc.perm[[i]])
}
