################################################################################
## Reichmann et al.
## Artificial light at night reorganizes habitat connectivity 
## and erodes landscape capacity to sustain biodiversity across South America
################################################################################
# This script is part of the code accompanying the manuscript above.
#
# 05 — Ratios, maps and exports
# Compares centrality metrics across scenarios and creates and exports the
# manuscript maps.

setwd(".")
load("Resultados_redes_D50.4.RData")
load("Paisaje_SA_EcoR_M_SA_J_LP.RData")

# J.permanente se utiliza solo para las capas de mapas de J PERM.
claves_J_permanente_mapa <- paste(M.SA.J.LP$CENTROID_X, M.SA.J.LP$CENTROID_Y, sep = "|")
J_permanente_mapa <- function(B, q = 0) {
  idx <- match(paste(B$CENTROID_X, B$CENTROID_Y, sep = "|"), claves_J_permanente_mapa)
  stopifnot(!anyNA(idx))
  M.SA.J.LP$J.permanente[idx] * (1 - B$LP)^q
}

library(ggplot2)
library(scales)
library(sf)
limites_ecoregions_sf <- st_as_sf(limites_ecoregions)
id <- 301:353 ; id <- setdiff(id, 351) 

##### LIGHT POLLUTION #####
ids <-301:353 # sumo 353 para Malvinas, que separe de Patagonia, es originalmente 348 tb
ids <- setdiff(ids, 351) # 351 Juan Fernandez 1 sola celda
all_points <- data.frame()

for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  B <- get(matrix_name)
  filtered_indices <- which(B$`EF(1).T(2).P(3)` == 100)
  if (length(filtered_indices) > 0) {
    all_points <- rbind(all_points, 
                        data.frame(CENTROID_X = B$CENTROID_X[filtered_indices],
                                   CENTROID_Y = B$CENTROID_Y[filtered_indices],
                                   LP = B$LP[filtered_indices],
                                   id = i))
  }
}


ggplot() +
  geom_point(data = M.SA.J, aes(x = CENTROID_X, y = CENTROID_Y), 
             cex = 0.15, col = alpha("grey", 0.5), pch = 15) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = LP), size = 0.8, shape = 15) +
  scale_color_gradientn(colors = c("#212C52", "#6B74F7", "#DE247B", "#FF7900", "#FFDE94"), name = "LP") +
  theme_minimal() +
  labs(title = "Light Pollution", x = "CENTROID_X", y = "CENTROID_Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA) 


jpeg("LP_export.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = M.SA.J, aes(x = CENTROID_X, y = CENTROID_Y), 
             cex = 0.15, col = alpha("#212C52", 0.7), pch = 15) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = LP), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#212C52", "#6B74F7", "#DE247B", "#FF7900", "#FFDE94"), name = "LP",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "Light Pollution", x = "CENTROID_X", y = "CENTROID_Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA) 
)
dev.off()





##### LAND #####
ids <-301:353 # sumo 353 para Malvinas, que separe de Patagonia, es originalmente 348 tb
ids <- setdiff(ids, 351) # 351 Juan Fernandez 1 sola celda
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  B <- get(matrix_name)
  all_points <- rbind(all_points, 
                      data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 100],
                                 CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 100],id = i))
}


ggplot() +
  geom_point(data = M.final.SA.LP.Q0.D50.4.LAND, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "#8FBC8F", size = 0.7, shape = 15) +
  theme_minimal() +
  labs(title = "Landscape", x = "CENTROID_X", y = "CENTROID_Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA) 


jpeg("LAND.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = M.final.SA.LP.Q0.D50.4.LAND, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "#8FBC8F", size = 0.7, shape = 15) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "Landscape", x = "CENTROID_X", y = "CENTROID_Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA) 
)
dev.off()





##### LAND QUARTILES #####
ids <- 301:353  # Sumo 353 para Malvinas, que separe de Patagonia, es originalmente 348 tb
ids <- setdiff(ids, 351)  # 351 Juan Fernandez 1 sola celda

all_points <- do.call(rbind, lapply(ids, function(i) {
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name)) {
    B <- get(matrix_name)
    if (!is.null(B$LAND) && any(!is.na(B$LAND))) {
      data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y, LAND = as.numeric(as.character(B$LAND)))
    }
  }
}))

if (nrow(all_points) > 0) {
  all_points <- all_points[!is.na(all_points$LAND), ]
  if (length(all_points$LAND) > 3) {
    all_points$LAND_category <- cut(all_points$LAND, 
                                    breaks = quantile(all_points$LAND, probs = seq(0, 1, 0.25), na.rm = TRUE, include.lowest = TRUE), 
                                    labels = c("Muy Bajo", "Bajo", "Medio", "Alto"), 
                                    include.lowest = TRUE)
    
    jpeg("LAND_tonos.jpeg", width = 8000, height = 8000, res = 300)
    print(
      ggplot() +
        geom_point(data = M.final.SA.LP.Q0.D50.4.LAND, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
        geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = LAND_category), size = 0.7, shape = 15) +
        scale_color_manual(values = c("Muy Bajo" = "#D0E8C6", "Bajo" = "#A3D9A5", "Medio" = "#73C78C", "Alto" = "#8FBC8F")) +
        theme_minimal() +
        theme(plot.title = element_text(size = 60, face = "bold"), axis.title = element_text(size = 30), axis.text = element_text(size = 30), 
              legend.title = element_text(size = 30), legend.text = element_text(size = 20)) +
        labs(title = "Landscape", x = "CENTROID_X", y = "CENTROID_Y", color = "LAND Nivel") +
        geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)
    )
    dev.off()
  }
}


##### PERM #####
all_points <- data.frame()
all_background <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y, id = i))
  filtered_indices <- which(B$`EF(1).T(2).P(3)` == 3)
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[filtered_indices], CENTROID_Y = B$CENTROID_Y[filtered_indices], id = i))
}

ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "grey", size = 0.7, shape = 15, alpha = 0.5) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y), color = "#4682B4", size = 0.7, shape = 15) +
  theme_minimal() +
  labs(title = "Waterscape", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)



jpeg("PERM.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "grey", size = 0.7, shape = 15, alpha = 0.5) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y), color = "#4682B4", size = 0.7, shape = 15) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "Waterscape", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)
)
dev.off()

##### PERM CUANTILES #####
all_points <- data.frame()
all_background <- data.frame()

for (i in ids) {
  B <- get(paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i))
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y, id = i))
  filtered_indices <- which(B$`EF(1).T(2).P(3)` == 3)
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[filtered_indices], CENTROID_Y = B$CENTROID_Y[filtered_indices], MAR = as.numeric(as.character(B$MAR[filtered_indices]))))
}

if (nrow(all_points) > 0) {
  all_points <- all_points[!is.na(all_points$MAR), ]
  if (length(all_points$MAR) > 3) {  
    breaks <- unique(quantile(all_points$MAR, probs = seq(0, 1, 0.25), na.rm = TRUE))
    if (length(breaks) < 4) breaks <- seq(min(all_points$MAR, na.rm = TRUE), max(all_points$MAR, na.rm = TRUE), length.out = 5)
    all_points$MAR_category <- cut(all_points$MAR, breaks = breaks, labels = c("Muy Bajo", "Bajo", "Medio", "Alto"), include.lowest = TRUE)
    colores <- c("Muy Bajo" = "#B3DAF1", "Bajo" = "#6CA6CD", "Medio" = "#4682B4", "Alto" = "#1E4B7E")
    jpeg("PERM_tonos.jpeg", width = 8000, height = 8000, res = 300)
    print(ggplot() +
            geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "grey", size = 0.7, shape = 15, alpha = 0.5) +
            geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = MAR_category), size = 0.7, shape = 15) +
            scale_color_manual(values = colores) +
            theme_minimal() +
            theme(plot.title = element_text(size = 60, face = "bold"),
                  axis.title = element_text(size = 30),
                  axis.text = element_text(size = 30),
                  legend.title = element_text(size = 30),
                  legend.text = element_text(size = 20)) +
            labs(title = "Waterscape", x = "CENTROID_X", y = "CENTROID_Y", color = "MAR Nivel") +
            geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA))
    dev.off()
  }
}

###### MAPAS J #####
##### GGPLOT J - PERM #####
all_background <- data.frame()
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y))
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 3],
                                             CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 3],
                                             J.uk = J_permanente_mapa(B, 0)[B$`EF(1).T(2).P(3)` == 3],id = i))
}

ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.uk), size = 0.3, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=0)" , 
                        breaks = c(min(all_points$J.uk, na.rm = TRUE), max(all_points$J.uk, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =1, barheight = 5)) +
  theme_minimal() +
  labs(title = "J.q0 PERM", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)

jpeg("Jq0_PERM.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.uk), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=0)" , 
                        breaks = c(min(all_points$J.uk, na.rm = TRUE), max(all_points$J.uk, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "J.q0 PERM", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)
)
dev.off()



all_background <- data.frame()
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q025.D50.4.PERM.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y))
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 3],
                                             CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 3],
                                             J.q025 = J_permanente_mapa(B, 0.25)[B$`EF(1).T(2).P(3)` == 3],id = i))
}

jpeg("Jq025_PERM.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.q025), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=0.25)",
                        breaks = c(min(all_points$J.q025, na.rm = TRUE), max(all_points$J.q025, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "J.q025 PERM", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.005, fill = NA)
)
dev.off()


all_background <- data.frame()
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q1.D50.4.PERM.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y))
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 3],
                                             CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 3],
                                             J.q1 = J_permanente_mapa(B, 1)[B$`EF(1).T(2).P(3)` == 3],id = i))
}

jpeg("Jq1_PERM.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.q1), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=1)",
                        breaks = c(min(all_points$J.q1, na.rm = TRUE), max(all_points$J.q1, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "J.q1 PERM", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.005, fill = NA)
)
dev.off()



all_background <- data.frame()
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y))
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 3],
                                             CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 3],
                                             J.q5 = J_permanente_mapa(B, 5)[B$`EF(1).T(2).P(3)` == 3],id = i))
}

jpeg("Jq5_PERM.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.q5), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=5)",
                        breaks = c(min(all_points$J.q5, na.rm = TRUE), max(all_points$J.q5, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "J.q5 PERM", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.005, fill = NA)
  )
dev.off()


##### GGPLOT J - LAND #####
all_background <- data.frame()
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y))
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 100],
                                             CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 100],
                                             J.uk = B$J.uk[B$`EF(1).T(2).P(3)` == 100],id = i))
}

jpeg("HOY_Jq0_LAND.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.uk), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=0)",
                        breaks = c(min(all_points$J.uk, na.rm = TRUE), max(all_points$J.uk, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "J.q0 LAND", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.005, fill = NA)
)
dev.off()



all_background <- data.frame()
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q025.D50.4.LAND.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y))
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 100],
                                             CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 100],
                                             J.q025 = B$J.q025[B$`EF(1).T(2).P(3)` == 100],id = i))
}
jpeg("Jq025_LAND.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.q025), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=0.25)",
                        breaks = c(min(all_points$J.q025, na.rm = TRUE), max(all_points$J.q025, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "J.q025 LAND", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.005, fill = NA)
)
dev.off()


all_background <- data.frame()
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q1.D50.4.LAND.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y))
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 100],
                                             CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 100],
                                             J.q1 = B$J.q1[B$`EF(1).T(2).P(3)` == 100],id = i))
}

jpeg("Jq1_LAND.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.q1), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=1)",
                        breaks = c(min(all_points$J.q1, na.rm = TRUE), max(all_points$J.q1, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "J.q1 LAND", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.005, fill = NA)
)
dev.off()



all_background <- data.frame()
all_points <- data.frame()
for (i in ids) {
  matrix_name <- paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_", i)
  B <- get(matrix_name)
  all_background <- rbind(all_background, data.frame(CENTROID_X = B$CENTROID_X, CENTROID_Y = B$CENTROID_Y))
  all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 100],
                                             CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 100],
                                             J.q5 = B$J.q5[B$`EF(1).T(2).P(3)` == 100],id = i))
}

jpeg("Jq5_LAND.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "grey", size = 1, shape = 15, alpha = 0.8) +
  geom_point(data = all_points, aes(x = CENTROID_X, y = CENTROID_Y, color = J.q5), size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#FFF9C4", "#FFEE58", "#FFD54F", "#FFB74D", "#FF9800", "#F57C00", "#EF6C00", "#E65100", "red4"), name = "J (Q=5)",
                        breaks = c(min(all_points$J.q5, na.rm = TRUE), max(all_points$J.q5, na.rm = TRUE)),
                        labels = c("0", "max"),
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "J.q5 LAND", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.005, fill = NA)
)
dev.off()





###### FUNCIONES ######
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

###### Ratios para LAND ######
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
###### Ratios para PERM ######
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

############################################### GGPLOT #####
library(ggplot2)
library(dplyr)
library(sf)
limites_ecoregions_sf <- st_as_sf(limites_ecoregions)

##### GGPLOT PERM ratio IN Q5/Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()
all_background <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.in.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
    
    B_df$ValueType <- ifelse(B_df$P5 > 0, "Positive", ifelse(B_df$P5 < 0, "Negative", "Zero"))
    B_df <- B_df[is.finite(B_df$P5), ]
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "P5", "ValueType")])
    all_background <- rbind(all_background, as.data.frame(map_data[, c("CENTROID_X", "CENTROID_Y")]))
  }
}

jpeg("PERM_ratio_IN_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, ValueType == "Positive"), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = P5), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, ValueType == "Zero"), 
             aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "#4682B4") +
  scale_color_gradientn(colors = colors, name = "P5 Ratio", guide = guide_colorbar (barwidth =3, barheight = 30)) +
  theme_minimal() +  
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 20),  # Títulos de los ejes
    axis.text = element_text(size = 20),  # Texto de los ejes
    legend.title = element_text(size = 20),  # Título de la leyenda
    legend.text = element_text(size =20)  # Texto de la leyenda
  ) +
  labs(title = "PERM ratio IN Q5/Q0", x = "Centroid X", y = "Centroid Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()
# contar cedas cero
# Crear dataframe vacío para almacenar los resultados
all_points <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.in.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    
    # Corregimos el nombre de columnas
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P5", "P025", "P5_2")
    
    # Agregamos la columna EcoR_ID correctamente
    B_df$EcoR_ID <- i
    
    # Clasificar en ValueType
    B_df$ValueType <- ifelse(B_df$P5 > 0, "Positive", ifelse(B_df$P5 < 0, "Negative", "Zero"))
    
    # Filtrar valores finitos
    B_df <- B_df[is.finite(B_df$P5), ]
    
    # Agregar a all_points
    all_points <- rbind(all_points, B_df[, c("EcoR_ID", "CENTROID_X", "CENTROID_Y", "P5", "ValueType")])
  }
}

# Ahora sí podemos contar correctamente
count_results <- data.frame()

for (i in id) {
  # Filtrar los datos de all_points para la EcoR_ID actual
  subset_data <- subset(all_points, EcoR_ID == i)
  
  # Contar cuántos son "Zero"
  zero_count <- sum(subset_data$ValueType == "Zero", na.rm = TRUE)
  total_count <- nrow(subset_data)
  
  # Agregar resultados al dataframe
  count_results <- rbind(count_results, data.frame(EcoR_ID = i, Zero_Count = zero_count, Total_Count = total_count))
}

# Calcular suma total
sum(count_results$Zero_Count)
sum(count_results$Total_Count)
(sum(count_results$Zero_Count) / sum(count_results$Total_Count)) * 100

# Mostrar resultados
print(count_results)

##### GGPLOT PERM ratio IN Q025/Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()
all_background <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.in.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
    
    B_df$ValueType <- ifelse(B_df$P025 > 0, "Positive", ifelse(B_df$P025 < 0, "Negative", "Zero"))
    B_df <- B_df[is.finite(B_df$P025), ]
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "P025", "ValueType")])
    all_background <- rbind(all_background, as.data.frame(map_data[, c("CENTROID_X", "CENTROID_Y")]))
  }
}

jpeg("PERM_ratio_IN_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, ValueType == "Positive"), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = P025), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, ValueType == "Zero"), 
             aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "black") +
  scale_color_gradientn(colors = colors, name = "P025 Ratio", guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "PERM ratio IN Q025/Q0", x = "Centroid X", y = "Centroid Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()

##### GGPLOT LAND ratio IN Q5/Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()
all_background <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.in.D50.4.LAND.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
    
    B_df$ValueType <- ifelse(B_df$P5 > 0, "Positive", ifelse(B_df$P5 < 0, "Negative", "Zero"))
    B_df <- B_df[is.finite(B_df$P5), ]
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "P5", "ValueType")])
    all_background <- rbind(all_background, as.data.frame(map_data[, c("CENTROID_X", "CENTROID_Y")]))
  }
}

jpeg("LAND_ratio_IN_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, ValueType == "Positive"), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = P5), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, ValueType == "Zero"), 
             aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "#8FBC8F") +
  scale_color_gradientn(colors = colors, name = "P5 Ratio",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND ratio IN Q5/Q0", x = "Centroid X", y = "Centroid Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()

# Crear dataframe vacío para almacenar los resultados
all_points <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.in.D50.4.LAND.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    
    # Corregimos el nombre de columnas
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P5", "P025", "P5_2")
    
    # Agregamos la columna EcoR_ID correctamente
    B_df$EcoR_ID <- i
    
    # Clasificar en ValueType
    B_df$ValueType <- ifelse(B_df$P5 > 0, "Positive", ifelse(B_df$P5 < 0, "Negative", "Zero"))
    
    # Filtrar valores finitos
    B_df <- B_df[is.finite(B_df$P5), ]
    
    # Agregar a all_points
    all_points <- rbind(all_points, B_df[, c("EcoR_ID", "CENTROID_X", "CENTROID_Y", "P5", "ValueType")])
  }
}

# Ahora sí podemos contar correctamente
count_results <- data.frame()

for (i in id) {
  # Filtrar los datos de all_points para la EcoR_ID actual
  subset_data <- subset(all_points, EcoR_ID == i)
  
  # Contar cuántos son "Zero"
  zero_count <- sum(subset_data$ValueType == "Zero", na.rm = TRUE)
  total_count <- nrow(subset_data)
  
  # Agregar resultados al dataframe
  count_results <- rbind(count_results, data.frame(EcoR_ID = i, Zero_Count = zero_count, Total_Count = total_count))
}

# Calcular suma total
sum(count_results$Zero_Count)
sum(count_results$Total_Count)
(sum(count_results$Zero_Count) / sum(count_results$Total_Count)) * 100

# Mostrar resultados
print(count_results)

##### GGPLOT LAND ratio IN Q025/Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()
all_background <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.in.D50.4.LAND.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
    
    B_df$ValueType <- ifelse(B_df$P025 > 0, "Positive", ifelse(B_df$P025 < 0, "Negative", "Zero"))
    B_df <- B_df[is.finite(B_df$P025), ]
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "P025", "ValueType")])
    all_background <- rbind(all_background, as.data.frame(map_data[, c("CENTROID_X", "CENTROID_Y")]))
  }
}

jpeg("LAND_ratio_IN_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, ValueType == "Positive"), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = P025), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, ValueType == "Zero"), 
             aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "#8FBC8F") +
  scale_color_gradientn(colors = colors, name = "P025 Ratio",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND ratio IN Q025/Q0", x = "Centroid X", y = "Centroid Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()

##### GGPLOT PERM ratio OUT Q5/Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()
all_background <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.out.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P5", "P025", "P5")
    
    B_df$ValueType <- ifelse(B_df$P5 > 0, "Positive", ifelse(B_df$P5 < 0, "Negative", "Zero"))
    B_df <- B_df[is.finite(B_df$P5), ]
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "P5", "ValueType")])
    all_background <- rbind(all_background, as.data.frame(map_data[, c("CENTROID_X", "CENTROID_Y")]))
  }
}

jpeg("00_PERM_ratio_OUT_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, ValueType == "Positive"), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = P5), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, ValueType == "Zero"), 
             aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "#758283") +
  scale_color_gradientn(colors = colors, name = "P5 Ratio",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "PERM ratio OUT Q5/Q0", x = "Centroid X", y = "Centroid Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


# contar cedas cero
# Crear dataframe vacío para almacenar los resultados
all_points <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.out.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    
    # Corregimos el nombre de columnas
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P5", "P025", "P5_2")
    
    # Agregamos la columna EcoR_ID correctamente
    B_df$EcoR_ID <- i
    
    # Clasificar en ValueType
    B_df$ValueType <- ifelse(B_df$P5 > 0, "Positive", ifelse(B_df$P5 < 0, "Negative", "Zero"))
    
    # Filtrar valores finitos
    B_df <- B_df[is.finite(B_df$P5), ]
    
    # Agregar a all_points
    all_points <- rbind(all_points, B_df[, c("EcoR_ID", "CENTROID_X", "CENTROID_Y", "P5", "ValueType")])
  }
}

# Ahora sí podemos contar correctamente
count_results <- data.frame()

for (i in id) {
  # Filtrar los datos de all_points para la EcoR_ID actual
  subset_data <- subset(all_points, EcoR_ID == i)
  
  # Contar cuántos son "Zero"
  zero_count <- sum(subset_data$ValueType == "Zero", na.rm = TRUE)
  total_count <- nrow(subset_data)
  
  # Agregar resultados al dataframe
  count_results <- rbind(count_results, data.frame(EcoR_ID = i, Zero_Count = zero_count, Total_Count = total_count))
}

# Calcular suma total
sum(count_results$Zero_Count)
sum(count_results$Total_Count)
(sum(count_results$Zero_Count) / sum(count_results$Total_Count)) * 100

# Mostrar resultados
print(count_results)


##### GGPLOT PERM ratio OUT Q025/Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()
all_background <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.out.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
    
    B_df$ValueType <- ifelse(B_df$P025 > 0, "Positive", ifelse(B_df$P025 < 0, "Negative", "Zero"))
    B_df <- B_df[is.finite(B_df$P025), ]
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "P025", "ValueType")])
    all_background <- rbind(all_background, as.data.frame(map_data[, c("CENTROID_X", "CENTROID_Y")]))
  }
}

jpeg("PERM_ratio_OUT_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, ValueType == "Positive"), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = P025), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, ValueType == "Zero"), 
             aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "#758283") +
  scale_color_gradientn(colors = colors, name = "P025 Ratio",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "PERM ratio OUT Q025/Q0", x = "Centroid X", y = "Centroid Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


##### GGPLOT LAND ratio OUT Q5/Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()
all_background <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.out.D50.4.LAND.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
    
    B_df$ValueType <- ifelse(B_df$P5 > 0, "Positive", ifelse(B_df$P5 < 0, "Negative", "Zero"))
    B_df <- B_df[is.finite(B_df$P5), ]
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "P5", "ValueType")])
    all_background <- rbind(all_background, as.data.frame(map_data[, c("CENTROID_X", "CENTROID_Y")]))
  }
}

jpeg("LAND_ratio_OUT_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, ValueType == "Positive"), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = P5), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, ValueType == "Zero"), 
             aes(x = CENTROID_X, y = CENTROID_Y), size = 0.1, shape = 15, color = "grey") +
  scale_color_gradientn(colors = colors, name = "P5 Ratio",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND ratio OUT Q5/Q0", x = "Centroid X", y = "Centroid Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()

# contar cedas cero
# Crear dataframe vacío para almacenar los resultados
all_points <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.out.D50.4.LAND.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    
    # Corregimos el nombre de columnas
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P5", "P025", "P5_2")
    
    # Agregamos la columna EcoR_ID correctamente
    B_df$EcoR_ID <- i
    
    # Clasificar en ValueType
    B_df$ValueType <- ifelse(B_df$P5 > 0, "Positive", ifelse(B_df$P5 < 0, "Negative", "Zero"))
    
    # Filtrar valores finitos
    B_df <- B_df[is.finite(B_df$P5), ]
    
    # Agregar a all_points
    all_points <- rbind(all_points, B_df[, c("EcoR_ID", "CENTROID_X", "CENTROID_Y", "P5", "ValueType")])
  }
}

# Ahora sí podemos contar correctamente
count_results <- data.frame()

for (i in id) {
  # Filtrar los datos de all_points para la EcoR_ID actual
  subset_data <- subset(all_points, EcoR_ID == i)
  
  # Contar cuántos son "Zero"
  zero_count <- sum(subset_data$ValueType == "Zero", na.rm = TRUE)
  total_count <- nrow(subset_data)
  
  # Agregar resultados al dataframe
  count_results <- rbind(count_results, data.frame(EcoR_ID = i, Zero_Count = zero_count, Total_Count = total_count))
}

# Calcular suma total
sum(count_results$Zero_Count)
sum(count_results$Total_Count)
(sum(count_results$Zero_Count) / sum(count_results$Total_Count)) * 100

# Mostrar resultados
print(count_results)

##### GGPLOT LAND ratio OUT Q025/Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()
all_background <- data.frame()

for (i in id) {
  matrix_name <- paste0("ratio.out.D50.4.LAND.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    B_df <- as.data.frame(B)
    colnames(B_df) <- c("CENTROID_X", "CENTROID_Y", "EcoR_ID", "Q0", "Q1", "Q025", "Q5", "P1", "P025", "P5")
    
    B_df$ValueType <- ifelse(B_df$P025 > 0, "Positive", ifelse(B_df$P025 < 0, "Negative", "Zero"))
    B_df <- B_df[is.finite(B_df$P025), ]
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "P025", "ValueType")])
    all_background <- rbind(all_background, as.data.frame(map_data[, c("CENTROID_X", "CENTROID_Y")]))
  }
}

jpeg("LAND_ratio_OUT_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, ValueType == "Positive"), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = P025), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, ValueType == "Zero"), 
             aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "gray") +
  scale_color_gradientn(colors = colors, name = "P025 Ratio",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND ratio OUT Q025/Q0", x = "Centroid X", y = "Centroid Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.7, fill = NA)
)
dev.off()

########### BC LAND Q5/Q0 ########## #####
id <- 301:353 ; id <- setdiff(id, 351)  
all_points <- all_background <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric(), Value = numeric(), Type = character())

for (i in id) {
  matrix_name <- paste0("ratiolog.bc.D50.4.LAND.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    pos <- which(B$Q5 > 0 & B$Q5 <= log(10))
    mayores <- which(B$Q5 > log(10))
    neg <- which(B$Q5 < 0 & B$Q5 >= log(0.10))
    menores <- which(B$Q5 < log(0.10))
    zeros <- which(abs(B$Q5) < 1e-10)
    if (length(pos) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[pos], CENTROID_Y = B$CENTROID_Y[pos], Value = B$Q5[pos], Type = "Positive"))
    if (length(mayores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[mayores], CENTROID_Y = B$CENTROID_Y[mayores], Value = NA, Type = "Mayores"))
    if (length(neg) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[neg], CENTROID_Y = B$CENTROID_Y[neg], Value = B$Q5[neg], Type = "Negative"))
    if (length(menores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[menores], CENTROID_Y = B$CENTROID_Y[menores], Value = NA, Type = "Menores"))
    if (length(zeros) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[zeros], CENTROID_Y = B$CENTROID_Y[zeros], Value = 0, Type = "Zero"))
    all_background <- rbind(all_background, data.frame(CENTROID_X = map_data$CENTROID_X, CENTROID_Y = map_data$CENTROID_Y))
  }
}


jpeg("LAND_logBC_pos_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, Type %in% c("Positive", "Mayores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("gray","#C6DBEF","#A5C9E0", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(0, log(10)), name = "Q5 LogBC",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND LogratioBC Positives Q5/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()

jpeg("LAND_logBC_neg_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, Type %in% c("Negative", "Menores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","#FCD2C0","gray"), limits = c(log(0.10), 0), name = "Q5 LogBC",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND LogratioBC Negatives Q5/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()

jpeg("LAND_logBC_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Positive"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Negative"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","gray","#A5C9E0", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(log(0.10), log(10)), name = "Q5 LogBC",
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.7, shape = 15) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND LogratioBC Q5/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()

# Crear dataframe para almacenar los resultados
count_results <- data.frame(EcoR_ID = numeric(), Zero_Count = numeric(), Total_Count = numeric())

for (i in id) {
  # Filtrar los datos de all_points para la EcoR_ID actual
  subset_data <- subset(all_points, Type %in% c("Positive", "Negative", "Zero") & CENTROID_X %in% get(paste0("ratiolog.bc.D50.4.LAND.ID_", i))$CENTROID_X)
  
  # Contar cuántos son "Zero"
  zero_count <- sum(subset_data$Type == "Zero", na.rm = TRUE)
  
  # Contar el total de datos para esa EcoR_ID (sin incluir "Mayores" y "Menores")
  total_count <- nrow(subset_data)
  
  # Agregar resultados al dataframe
  count_results <- rbind(count_results, data.frame(EcoR_ID = i, Zero_Count = zero_count, Total_Count = total_count))
}

# Calcular porcentaje de "Zero" sobre el total
count_results$Percent_Zero <- (count_results$Zero_Count / count_results$Total_Count) * 100

# Manejar el caso donde Total_Count sea 0 (evitar NaN o Inf)
count_results$Percent_Zero[is.na(count_results$Percent_Zero) | is.infinite(count_results$Percent_Zero)] <- 0

# Calcular suma total
total_zeros <- sum(count_results$Zero_Count)
total_points <- sum(count_results$Total_Count)
percent_total_zero <- (total_zeros / total_points) * 100

# Agregar fila de totales al dataframe
count_results <- rbind(count_results, data.frame(EcoR_ID = "Total", Zero_Count = total_zeros, Total_Count = total_points, Percent_Zero = percent_total_zero))

# Mostrar resultados
print(count_results)



########### BC LAND Q025/Q0 ########## #####
id <- 301:353 ; id <- setdiff(id, 351) 
all_points <- all_background <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric(), Value = numeric(), Type = character())

for (i in id) {
  matrix_name <- paste0("ratiolog.bc.D50.4.LAND.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    pos <- which(B$Q025 > 0 & B$Q025 <= log(10))
    mayores <- which(B$Q025 > log(10))
    neg <- which(B$Q025 < 0 & B$Q025 >= log(0.10))
    menores <- which(B$Q025 < log(0.10))
    zeros <- which(abs(B$Q025) < 1e-10)
    if (length(pos) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[pos], CENTROID_Y = B$CENTROID_Y[pos], Value = B$Q025[pos], Type = "Positive"))
    if (length(mayores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[mayores], CENTROID_Y = B$CENTROID_Y[mayores], Value = NA, Type = "Mayores"))
    if (length(neg) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[neg], CENTROID_Y = B$CENTROID_Y[neg], Value = B$Q025[neg], Type = "Negative"))
    if (length(menores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[menores], CENTROID_Y = B$CENTROID_Y[menores], Value = NA, Type = "Menores"))
    if (length(zeros) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[zeros], CENTROID_Y = B$CENTROID_Y[zeros], Value = 0, Type = "Zero"))
    all_background <- rbind(all_background, data.frame(CENTROID_X = map_data$CENTROID_X, CENTROID_Y = map_data$CENTROID_Y))
  }
}

jpeg("LAND_logBC_pos_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, Type %in% c("Positive", "Mayores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("gray","#C6DBEF","#A5C9E0", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(0, log(10)), name = "Q025 LogBC",
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND LogratioBC Positives Q025/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


jpeg("LAND_logBC_neg_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = subset(all_points, Type %in% c("Negative", "Menores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","#FCD2C0","gray"), limits = c(log(0.10), 0), name = "Q025 LogBC",guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND LogratioBC Negatives Q025/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


jpeg("LAND_logBC_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Positive"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Negative"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","gray","#A5C9E0", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(log(0.10), log(10)), name = "Q025 LogBC",
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.7, shape = 15) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "LAND LogratioBC Q025/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()

########### BC PERM Q5/Q0 ########## #####
id <- 301:353 ; id <- setdiff(id, 351) 
#id <- c(345,348) 
all_points <- all_background <- perm_cells <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric(), Value = numeric(), Type = character())

for (i in id) {
  matrix_name <- paste0("ratiolog.bc.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    pos <- which(B$Q5 > 0 & B$Q5 <= log(10))
    mayores <- which(B$Q5 > log(10))
    neg <- which(B$Q5 < 0 & B$Q5 >= log(0.10))
    menores <- which(B$Q5 < log(0.10))
    zeros <- which(abs(B$Q5) < 1e-10)
    if (length(pos) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[pos], CENTROID_Y = B$CENTROID_Y[pos], Value = B$Q5[pos], Type = "Positive"))
    if (length(mayores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[mayores], CENTROID_Y = B$CENTROID_Y[mayores], Value = NA, Type = "Mayores"))
    if (length(neg) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[neg], CENTROID_Y = B$CENTROID_Y[neg], Value = B$Q5[neg], Type = "Negative"))
    if (length(menores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[menores], CENTROID_Y = B$CENTROID_Y[menores], Value = NA, Type = "Menores"))
    if (length(zeros) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[zeros], CENTROID_Y = B$CENTROID_Y[zeros], Value = 0, Type = "Zero"))
    all_background <- rbind(all_background, data.frame(CENTROID_X = map_data$CENTROID_X, CENTROID_Y = map_data$CENTROID_Y))
    perm_cells <- rbind(perm_cells, data.frame(CENTROID_X = map_data$CENTROID_X[map_data$`EF(1).T(2).P(3).` == 3], CENTROID_Y = map_data$CENTROID_Y[map_data$`EF(1).T(2).P(3).` == 3]))
      }
}

jpeg("PERM_logBC_pos_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type %in% c("Positive", "Mayores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#758283","#C6DBEF","#A5C9E0", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(0, log(10)), name = "Q5 LogBC",
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "PERM LogratioBC Positives Q5/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


jpeg("PERM_logBC_neg_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type %in% c("Negative", "Menores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","#FCD2C0","#758283"), limits = c(log(0.10), 0), name = "Q5 LogBC",
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "PERM LogratioBC Negatives Q5/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


jpeg("PERM_logBC_Q5_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.7, shape = 15) +
  geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Positive"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Negative"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","#758283","#A5C9E0", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(log(0.10), log(10)), name = "Q5 LogBC",
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.7, shape = 15) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "PERM LogratioBC Q5/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


# Crear dataframe para almacenar los resultados
count_results <- data.frame(EcoR_ID = numeric(), Zero_Count = numeric(), Total_Count = numeric())

for (i in id) {
  # Filtrar los datos de all_points para la EcoR_ID actual
  subset_data <- subset(all_points, Type %in% c("Positive", "Negative", "Zero") & CENTROID_X %in% get(paste0("ratiolog.bc.D50.4.PERM.ID_", i))$CENTROID_X)
  
  # Contar cuántos son "Zero"
  zero_count <- sum(subset_data$Type == "Zero", na.rm = TRUE)
  
  # Contar el total de datos para esa EcoR_ID (sin incluir "Mayores" y "Menores")
  total_count <- nrow(subset_data)
  
  # Agregar resultados al dataframe
  count_results <- rbind(count_results, data.frame(EcoR_ID = i, Zero_Count = zero_count, Total_Count = total_count))
}

# Calcular porcentaje de "Zero" sobre el total
count_results$Percent_Zero <- (count_results$Zero_Count / count_results$Total_Count) * 100

# Manejar el caso donde Total_Count sea 0 (evitar NaN o Inf)
count_results$Percent_Zero[is.na(count_results$Percent_Zero) | is.infinite(count_results$Percent_Zero)] <- 0

# Calcular suma total
sum(count_results$Zero_Count)
sum(count_results$Total_Count)
(sum(count_results$Zero_Count) / sum(count_results$Total_Count)) * 100




########### BC PERM Q025/Q0 ########## #####
id <- 301:353 ; id <- setdiff(id, 351) 
#id <- c(345,348) 
all_points <- all_background <- perm_cells <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric(), Value = numeric(), Type = character())

for (i in id) {
  matrix_name <- paste0("ratiolog.bc.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  B <- get(matrix_name)
  map_data <- get(map_name)
  pos <- which(B$Q025 > 0 & B$Q025 <= log(10))
  mayores <- which(B$Q025 > log(10))
  neg <- which(B$Q025 < 0 & B$Q025 >= log(0.10))
  menores <- which(B$Q025 < log(0.10))
  zeros <- which(abs(B$Q025) < 1e-10)
  if (length(pos) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[pos], CENTROID_Y = B$CENTROID_Y[pos], Value = B$Q025[pos], Type = "Positive"))
  if (length(mayores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[mayores], CENTROID_Y = B$CENTROID_Y[mayores], Value = NA, Type = "Mayores"))
  if (length(neg) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[neg], CENTROID_Y = B$CENTROID_Y[neg], Value = B$Q025[neg], Type = "Negative"))
  if (length(menores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[menores], CENTROID_Y = B$CENTROID_Y[menores], Value = NA, Type = "Menores"))
  if (length(zeros) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[zeros], CENTROID_Y = B$CENTROID_Y[zeros], Value = 0, Type = "Zero"))
  all_background <- rbind(all_background, data.frame(CENTROID_X = map_data$CENTROID_X, CENTROID_Y = map_data$CENTROID_Y))
  perm_cells <- rbind(perm_cells, data.frame(CENTROID_X = map_data$CENTROID_X[map_data$`EF(1).T(2).P(3).` == 3], CENTROID_Y = map_data$CENTROID_Y[map_data$`EF(1).T(2).P(3).` == 3]))
}

jpeg("PERM_logBC_pos_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
  ggplot() +
    geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
    geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.7, shape = 15) +
    geom_point(data = subset(all_points, Type %in% c("Positive", "Mayores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
    geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.7, shape = 15) +
    scale_color_gradientn(colors = c("#758283","#C6DBEF","#A5C9E0", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(0, log(10)), name = "Q025 LogBC",
                          guide = guide_colorbar (barwidth =5, barheight = 50)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    labs(title = "PERM LogratioBC Positives Q025/Q0", x = "Centroid X", y = "Centroid Y") +
    geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.7, fill = NA)
)
dev.off()


jpeg("PERM_logBC_neg_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type %in% c("Negative", "Menores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","#FCD2C0","#758283"), limits = c(log(0.10), 0), name = "Q025 LogBC",
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "PERM LogratioBC Negatives Q025/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


jpeg("PERM_logBC_Q025_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.7, shape = 15, alpha = 0.8) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.7, shape = 15) +
  geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Positive"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Negative"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.7, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.7, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","#758283","#A5C9E0", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(log(0.10), log(10)), name = "Q025 LogBC",
                        guide = guide_colorbar (barwidth =5, barheight = 50)) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.7, shape = 15) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "PERM LogratioBC Q025/Q0", x = "Centroid X", y = "Centroid Y")+
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


########### BC PERM Q025/Q0 ########## #####
id <- setdiff(301:353, 351)  # Excluir ID 351
all_points <- all_background <- perm_cells <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric(), Value = numeric(), Type = character())

for (i in id) {
  matrix_name <- paste0("ratiolog.bc.D50.4.PERM.ID_", i)
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(matrix_name) && exists(map_name)) {
    B <- get(matrix_name)
    map_data <- get(map_name)
    pos <- which(B$Q025 > 0 & B$Q025 <= log(10))
    mayores <- which(B$Q025 > log(10))
    neg <- which(B$Q025 < 0 & B$Q025 >= log(0.10))
    menores <- which(B$Q025 < log(0.10))
    zeros <- which(abs(B$Q025) < 1e-10)
    if (length(pos) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[pos], CENTROID_Y = B$CENTROID_Y[pos], Value = B$Q025[pos], Type = "Positive"))
    if (length(mayores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[mayores], CENTROID_Y = B$CENTROID_Y[mayores], Value = NA, Type = "Mayores"))
    if (length(neg) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[neg], CENTROID_Y = B$CENTROID_Y[neg], Value = B$Q025[neg], Type = "Negative"))
    if (length(menores) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[menores], CENTROID_Y = B$CENTROID_Y[menores], Value = NA, Type = "Menores"))
    if (length(zeros) > 0) all_points <- rbind(all_points, data.frame(CENTROID_X = B$CENTROID_X[zeros], CENTROID_Y = B$CENTROID_Y[zeros], Value = 0, Type = "Zero"))
    all_background <- rbind(all_background, data.frame(CENTROID_X = map_data$CENTROID_X, CENTROID_Y = map_data$CENTROID_Y))
    perm_cells <- rbind(perm_cells, data.frame(CENTROID_X = map_data$CENTROID_X[map_data$`EF(1).T(2).P(3).` == 3],CENTROID_Y = map_data$CENTROID_Y[map_data$`EF(1).T(2).P(3).` == 3]))
  }
}

ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.3, shape = 15, alpha = 0.8) +
  geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.3, shape = 15) +
  geom_point(data = subset(all_points, Type %in% c("Positive", "Mayores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.3, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.2, shape = 15) +
  scale_color_gradientn(colors = c("white","#E3F2FD","#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(0, log(10)), name = "Q025 LogBC") +
  theme_minimal() +
  labs(title = "PERM LogratioBC Positives Q025/Q0", x = "Centroid X", y = "Centroid Y")

ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.3, shape = 15, alpha = 0.8) +
  geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.3, shape = 15) +
  geom_point(data = subset(all_points, Type %in% c("Negative", "Menores")), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.3, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.2, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","#FFF5F0","white"), limits = c(log(0.10), 0), name = "Q025 LogBC") +
  theme_minimal() +
  labs(title = "PERM LogratioBC Negatives Q025/Q0", x = "Centroid X", y = "Centroid Y")

ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), color = "gray", size = 0.3, shape = 15, alpha = 0.8) +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.2, shape = 15) +
  geom_point(data = perm_cells, aes(x = CENTROID_X, y = CENTROID_Y), color = "#758283", size = 0.3, shape = 15) +
  geom_point(data = subset(all_points, Type == "Positive"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.3, shape = 15) +
  geom_point(data = subset(all_points, Type == "Negative"), aes(x = CENTROID_X, y = CENTROID_Y, color = Value), size = 0.3, shape = 15) +
  geom_point(data = subset(all_points, Type == "Mayores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#08306B", size = 0.2, shape = 15) +
  geom_point(data = subset(all_points, Type == "Menores"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#A50F15", size = 0.2, shape = 15) +
  scale_color_gradientn(colors = c("#F8766D", "#FC9272", "#FDB294", "#FBB5A0","white","#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD"), limits = c(log(0.10), log(10)), name = "Q025 LogBC") +
  #geom_point(data = subset(all_points, Type == "Zero"), aes(x = CENTROID_X, y = CENTROID_Y), color = "#8FBC8F", size = 0.2, shape = 15) +
  theme_minimal() +
  labs(title = "PERM LogratioBC Q025/Q0", x = "Centroid X", y = "Centroid Y")

##### BC_resumen LAND #######
BC_resumen <- data.frame(EcoR_ID = numeric(), Q0_sum = numeric(), Q025_sum = numeric(), Q1_sum = numeric(), Q5_sum = numeric())
ids <- setdiff(301:353, 351)
for (i in ids) {
  B <- get(paste0("ratiolog.bc.D50.4.LAND.ID_", i))
  BC_resumen <- rbind(BC_resumen, data.frame(
    EcoR_ID = i,
    Q0_sum = sum(B$Q0, na.rm = TRUE),
    Q025_sum = sum(B$Q025, na.rm = TRUE),
    Q1_sum = sum(B$Q1, na.rm = TRUE),
    Q5_sum = sum(B$Q5, na.rm = TRUE)
  ))
}

# Generar puntos con datos
all_points <- data.frame()
for (i in ids) {
  B <- get(paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i))
  all_points <- rbind(all_points, data.frame(
    CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 100],
    CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 100],
    id = i
  ))
}

# Calcular magnitud absoluta y rango global (sin incluir SA)
BC_resumen$Q025_abs <- abs(BC_resumen$Q025_sum)
BC_resumen$Q5_abs <- abs(BC_resumen$Q5_sum)
BC_resumen_noSA <- BC_resumen[BC_resumen$EcoR_ID != "SA", ]
global_range <- range(c(BC_resumen_noSA$Q025_abs, BC_resumen_noSA$Q5_abs), na.rm = TRUE)

jpeg("LAND_logBC_resumen_Q025.jpeg", width = 8000, height = 8000, res = 300)
ggplot() +
  geom_point(data = M.final.SA.LP.Q0.D50.4.LAND, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.8, shape = 15, alpha = 0.8) +
  geom_point(data = merge(all_points, BC_resumen[, c("EcoR_ID", "Q025_abs")], by.x = "id", by.y = "EcoR_ID", all.x = TRUE), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = Q025_abs), 
             size = 0.8, shape = 15) +
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = global_range, guide = guide_colorbar (barwidth =5, barheight = 50),
    name = "|Loss| (Abs)", 
    oob = scales::squish) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "Magnitude of Q025 Loss in Ecoregions (PERM)", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)
dev.off()

jpeg("LAND_logBC_resumen_Q5.jpeg", width = 8000, height = 8000, res = 300)

ggplot() +
  geom_point(data = M.final.SA.LP.Q0.D50.4.LAND, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.8, shape = 15, alpha = 0.8) +
  geom_point(data = merge(all_points, BC_resumen[, c("EcoR_ID", "Q5_abs")], by.x = "id", by.y = "EcoR_ID", all.x = TRUE), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = Q5_abs), 
             size = 0.8, shape = 15) +
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = global_range, guide = guide_colorbar (barwidth =5, barheight = 50),
    name = "|Loss| (Abs)", 
    oob = scales::squish) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "Magnitude of Q5 Loss in Ecoregions (PERM)", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)

dev.off()


##### BC_resumen PERM #######
BC_resumen <- data.frame(EcoR_ID = numeric(), Q0_sum = numeric(), Q025_sum = numeric(), Q1_sum = numeric(), Q5_sum = numeric())
ids <- setdiff(301:353, 351)
for (i in ids) {
  B <- get(paste0("ratiolog.bc.D50.4.PERM.ID_", i))
  BC_resumen <- rbind(BC_resumen, data.frame(
    EcoR_ID = i,
    Q0_sum = sum(B$Q0, na.rm = TRUE),
    Q025_sum = sum(B$Q025, na.rm = TRUE),
    Q1_sum = sum(B$Q1, na.rm = TRUE),
    Q5_sum = sum(B$Q5, na.rm = TRUE)
  ))
}

BC_resumen$Q025_abs <- abs(BC_resumen$Q025_sum)
BC_resumen$Q5_abs <- abs(BC_resumen$Q5_sum)

id <- setdiff(301:353, 351)
all_points <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric(), id = numeric(), Value = numeric())
all_background <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric())

for (i in id) {
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(map_name)) {
    map_data <- get(map_name)
    if (!is.null(map_data$CENTROID_X) && !is.null(map_data$CENTROID_Y)) {
      all_background <- rbind(all_background, data.frame(CENTROID_X = map_data$CENTROID_X, CENTROID_Y = map_data$CENTROID_Y))
    }
  }
}

for (i in ids) {
  B <- get(paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i))
  if (!is.null(B$CENTROID_X) && !is.null(B$CENTROID_Y) && any(B$`EF(1).T(2).P(3)` == 3)) {
    all_points <- rbind(all_points, data.frame(
      CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 3],
      CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 3],
      id = i
    ))
  }
}

jpeg("PERM_logBC_resumen_Q025.jpeg", width = 8000, height = 8000, res = 300)

ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.8, shape = 15, alpha = 0.8) +
  geom_point(data = merge(all_points, BC_resumen[, c("EcoR_ID", "Q025_abs")], by.x = "id", by.y = "EcoR_ID", all.x = TRUE), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = Q025_abs), 
             size = 0.8, shape = 15) +
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = global_range, guide = guide_colorbar (barwidth =5, barheight = 50),
    name = "|Loss| (Abs)", 
    oob = scales::squish) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "Magnitude of Q025 Loss in Ecoregions (PERM)", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)

dev.off()

jpeg("PERM_logBC_resumen_Q5.jpeg", width = 8000, height = 8000, res = 300)

ggplot() +
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.8, shape = 15, alpha = 0.8) +
  geom_point(data = merge(all_points, BC_resumen[, c("EcoR_ID", "Q5_abs")], by.x = "id", by.y = "EcoR_ID", all.x = TRUE), 
             aes(x = CENTROID_X, y = CENTROID_Y, color = Q5_abs), 
             size = 0.8, shape = 15) +
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = global_range, guide = guide_colorbar (barwidth =5, barheight = 50),
    name = "|Loss| (Abs)", 
    oob = scales::squish) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20) # Texto de la leyenda
  ) +
  labs(title = "Magnitude of Q5 Loss in Ecoregions (PERM)", x = "CENTROID_X", y = "CENTROID_Y") +
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)

dev.off()



##### LAMBDA Q5 resumen LAND #####
library(readr)
stopifnot(exists("eigenvalues")) # Cargado desde Resultados_redes_D50.4.RData
ids <- setdiff(301:353, 351)

# Definir el rango global de los valores de dif_perc_LAND
global_range <- range(eigenvalues$dif_perc_LAND_q5, na.rm = TRUE)

# Generar puntos con datos
all_points <- data.frame()
for (i in ids) {
  B <- get(paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i))
  all_points <- rbind(all_points, data.frame(
    CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 100],
    CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 100],
    id = i
  ))
}

# Fusionar los puntos con eigenvalues para agregar dif_perc_LAND
merged_data <- merge( all_points, transform(eigenvalues[, c("EcoR_ID", "dif_perc_LAND_q5")],
            abs_dif_perc_LAND_q5 = abs(dif_perc_LAND_q5)),by.x = "id", by.y = "EcoR_ID", all.x = TRUE)


# Guardar la imagen como JPEG con el nuevo nombre
jpeg("paper_LAND_lambda_q5.jpeg", width = 8000, height = 8000, res = 300)

ggplot() +
  # Puntos de fondo en gris
  geom_point(data = M.final.SA.LP.Q0.D50.4.LAND, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.8, shape = 15, alpha = 0.8) +
  # Puntos coloreados según la magnitud de dif_perc_LAND
  geom_point(data = merged_data, aes(x = CENTROID_X, y = CENTROID_Y, color = abs_dif_perc_LAND_q5), 
             size = 0.8, shape = 15) +
  # Escala de colores con el rango corregido
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = c(0,35), guide = guide_colorbar(barwidth = 2, barheight = 25),
    name = "Percentage Loss", 
    oob = scales::squish) +
  # Estilo minimalista
  theme_minimal() +
  theme(
    plot.title = element_text(size = 40, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.position= c(0.8, 0.18), # x e y entre 0 y 1 (relativo al gráfico)
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20)  # Texto de la leyenda
  ) +
  # Etiquetas
  labs(title = "Proportional change in ecoregion carrying capacity (LAND - Q5)", x = "", y = "") +
  # Límites de ecoregiones
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.5, fill = NA)
dev.off()


##### LAMBDA Q5 resumen PERM #####
ids <- setdiff(301:353, 351)

# Definir el rango global de los valores de dif_perc_PERM
global_range <- range(eigenvalues$dif_perc_PERM_q5, na.rm = TRUE)

# Generar puntos con datos
all_points <- data.frame()
for (i in ids) {
  B <- get(paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i))
  if (!is.null(B$CENTROID_X) && !is.null(B$CENTROID_Y) && any(B$`EF(1).T(2).P(3)` == 3)) {
    all_points <- rbind(all_points, data.frame(
      CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 3],
      CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 3],
      id = i
    ))
  }
}

all_background <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric())
for (i in ids) {
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(map_name)) {
    map_data <- get(map_name)
    if (!is.null(map_data$CENTROID_X) && !is.null(map_data$CENTROID_Y)) {
      all_background <- rbind(all_background, data.frame(CENTROID_X = map_data$CENTROID_X, CENTROID_Y = map_data$CENTROID_Y))
    }
  }
}

# Fusionar los puntos con eigenvalues para agregar dif_perc_PERM

merged_data <- merge( all_points, transform(eigenvalues[, c("EcoR_ID", "dif_perc_PERM_q5")],
                                            abs_dif_perc_PERM_q5 = abs(dif_perc_PERM_q5)),by.x = "id", by.y = "EcoR_ID", all.x = TRUE)

jpeg("paper_PERM_lambda_q5.jpeg", width = 8000, height = 8000, res = 300)

ggplot() +
  # Puntos de fondo en gris para TODAS las EcoR_ID
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.8, shape = 15, alpha = 0.8) +
  # Puntos coloreados según la magnitud de dif_perc_PERM
  geom_point(data = merged_data, aes(x = CENTROID_X, y = CENTROID_Y, color = abs_dif_perc_PERM_q5), 
             size = 0.8, shape = 15) +
  # Escala de colores con el rango corregido
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = c(0,75), guide = guide_colorbar(barwidth = 2, barheight = 25),
    name = "Percentage Loss", 
    oob = scales::squish) +
  # Estilo minimalista
  theme_minimal() +
  theme(
    plot.title = element_text(size = 40, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.position= c(0.8, 0.18), # x e y entre 0 y 1 (relativo al gráfico)
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20)  # Texto de la leyenda
  ) +
  # Etiquetas
  labs(title = "Proportional change in ecoregion carrying capacity” (PERM - Q5)", x = "", y = "") +
  # Límites de ecoregiones
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 1, fill = NA)

dev.off()


##### ECOR_ COMPLETA #####
ids <- setdiff(301:353, 351)

# Preparamos el dataframe grande con TODAS las celdas, por EcoR_ID
all_regions <- data.frame()

for (i in ids) {
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  
  if (exists(map_name)) {
    B <- get(map_name)
    
    if (!is.null(B$CENTROID_X) && !is.null(B$CENTROID_Y)) {
      temp <- data.frame(
        CENTROID_X = B$CENTROID_X,
        CENTROID_Y = B$CENTROID_Y,
        id = i
      )
      all_regions <- rbind(all_regions, temp)
    }
  }
}

# Tomamos solo EcoR_ID + dif_perc_PERM_q5, y su valor absoluto
ev <- transform(
  eigenvalues[, c("EcoR_ID", "dif_perc_PERM_q5")],
  abs_dif_perc_PERM_q5 = abs(dif_perc_PERM_q5)
)

# Fusionamos con todas las celdas de cada EcoR_ID
merged_regions <- merge(
  all_regions,
  ev,
  by.x = "id",
  by.y = "EcoR_ID",
  all.x = TRUE
)

# Abrimos archivo
jpeg("paper_full_PERM_lambda_q5_ECOREGIONES.jpeg", width = 8000, height = 8000, res = 300)

ggplot() +
  # Fondo gris para todas las celdas sin clasificación
  geom_point(
    data = all_regions, 
    aes(x = CENTROID_X, y = CENTROID_Y),
    color = "gray80", size = 0.7, shape = 15
  ) +
  
  # Ahora pintamos TODAS las celdas de cada EcoR_ID según su pérdida (%)
  geom_point(
    data = merged_regions,
    aes(x = CENTROID_X, y = CENTROID_Y, color = abs_dif_perc_PERM_q5),
    size = 0.8, shape = 15
  ) +
  
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = c(0, 75),
    guide = guide_colorbar(barwidth = 2, barheight = 25),
    name = "Percentage Loss",
    oob = scales::squish
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 40, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.position = c(0.8, 0.18),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 20)
  ) +
  
  labs(
    title = "Proportional change in ecoregion carrying capacity (PERM - Q5)",
    x = "", 
    y = ""
  ) +
  
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 1, fill = NA)

dev.off()


##### LAMBDA Q025 resumen LAND #####
ids <- setdiff(301:353, 351)

# Definir el rango global de los valores de dif_perc_LAND
global_range <- range(eigenvalues$dif_perc_LAND_q025, na.rm = TRUE)

# Generar puntos con datos
all_points <- data.frame()
for (i in ids) {
  B <- get(paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i))
  all_points <- rbind(all_points, data.frame(
    CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 100],
    CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 100],
    id = i
  ))
}

# Fusionar los puntos con eigenvalues para agregar dif_perc_LAND
merged_data <- merge( all_points, transform(eigenvalues[, c("EcoR_ID", "dif_perc_LAND_q025")],
                                            abs_dif_perc_LAND_q025 = abs(dif_perc_LAND_q025)),by.x = "id", by.y = "EcoR_ID", all.x = TRUE)


# Guardar la imagen como JPEG con el nuevo nombre
jpeg("HOYLAND_lambda_q025.jpeg", width = 8000, height = 8000, res = 300)

ggplot() +
  # Puntos de fondo en gris
  geom_point(data = M.final.SA.LP.Q0.D50.4.LAND, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.8, shape = 15, alpha = 0.8) +
  # Puntos coloreados según la magnitud de dif_perc_LAND
  geom_point(data = merged_data, aes(x = CENTROID_X, y = CENTROID_Y, color = abs_dif_perc_LAND_q025), 
             size = 0.8, shape = 15) +
  # Escala de colores con el rango corregido
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = c(0,35), guide = guide_colorbar(barwidth = 2, barheight = 25),
    name = "Percentage Loss", 
    oob = scales::squish) +
  # Estilo minimalista
  theme_minimal() +
  theme(
    plot.title = element_text(size = 40, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.position= c(0.8, 0.18), # x e y entre 0 y 1 (relativo al gráfico)
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20)  # Texto de la leyenda
  ) +
  # Etiquetas
  labs(title = "Proportional change in ecoregion carrying capacity” (LAND - Q025)", x = "CENTROID_X", y = "CENTROID_Y") +
  # Límites de ecoregiones
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)

dev.off()


##### LAMBDA Q025 resumen PERM #####
ids <- setdiff(301:353, 351)

# Definir el rango global de los valores de dif_perc_PERM
global_range <- range(eigenvalues$dif_perc_PERM_q025, na.rm = TRUE)

# Generar puntos con datos
all_points <- data.frame()
for (i in ids) {
  B <- get(paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i))
  if (!is.null(B$CENTROID_X) && !is.null(B$CENTROID_Y) && any(B$`EF(1).T(2).P(3)` == 3)) {
    all_points <- rbind(all_points, data.frame(
      CENTROID_X = B$CENTROID_X[B$`EF(1).T(2).P(3)` == 3],
      CENTROID_Y = B$CENTROID_Y[B$`EF(1).T(2).P(3)` == 3],
      id = i
    ))
  }
}

all_background <- data.frame(CENTROID_X = numeric(), CENTROID_Y = numeric())
for (i in ids) {
  map_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  if (exists(map_name)) {
    map_data <- get(map_name)
    if (!is.null(map_data$CENTROID_X) && !is.null(map_data$CENTROID_Y)) {
      all_background <- rbind(all_background, data.frame(CENTROID_X = map_data$CENTROID_X, CENTROID_Y = map_data$CENTROID_Y))
    }
  }
}

# Fusionar los puntos con eigenvalues para agregar dif_perc_PERM
merged_data <- merge( all_points, transform(eigenvalues[, c("EcoR_ID", "dif_perc_PERM_q025")],
                                            abs_dif_perc_PERM_q025 = abs(dif_perc_PERM_q025)),by.x = "id", by.y = "EcoR_ID", all.x = TRUE)


jpeg("PERM_lambda_q025.jpeg", width = 8000, height = 8000, res = 300)

ggplot() +
  # Puntos de fondo en gris para TODAS las EcoR_ID
  geom_point(data = all_background, aes(x = CENTROID_X, y = CENTROID_Y), 
             color = "gray", size = 0.8, shape = 15, alpha = 0.8) +
  # Puntos coloreados según la magnitud de dif_perc_PERM
  geom_point(data = merged_data, aes(x = CENTROID_X, y = CENTROID_Y, color = abs_dif_perc_PERM_q025), 
             size = 0.8, shape = 15) +
  # Escala de colores con el rango corregido
  scale_color_gradientn(
    colors = c("#FFF5F0", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
    limits = c(0,75), guide = guide_colorbar(barwidth = 2, barheight = 25),
    name = "Percentage Loss", 
    oob = scales::squish) +
  # Estilo minimalista
  theme_minimal() +
  theme(
    plot.title = element_text(size = 40, face = "bold"),  # Título del gráfico
    axis.title = element_text(size = 30),  # Títulos de los ejes
    axis.text = element_text(size = 30),  # Texto de los ejes
    legend.position= c(0.8, 0.18), # x e y entre 0 y 1 (relativo al gráfico)
    legend.title = element_text(size = 30),  # Título de la leyenda
    legend.text = element_text(size = 20)  # Texto de la leyenda
  ) +
  # Etiquetas
  labs(title = "Proportional change in ecoregion carrying capacity” (PERM - Q025)", x = "CENTROID_X", y = "CENTROID_Y") +
  # Límites de ecoregiones
  geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.05, fill = NA)
dev.off()




###### POST DEFENSA Q DIRECTOS #####
##### GGPLOT LAND grado OUT Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()

for (i in id) {
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name)) {
    B <- get(matrix_name)
    B_df <- as.data.frame(B)
    
    B_df$grado_out <- B_df$grado.out.D50.4
    B_df <- B_df[B_df$LAND > 0 & B_df$J.uk > 0, ]
    B_df <- B_df[is.finite(B_df$grado_out), ]
    B_df <- B_df[!duplicated(B_df[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID")]), ]
    
    B_df$ValueType <- ifelse(B_df$grado_out > 0, "Positive", ifelse(B_df$grado_out < 0, "Negative", "Zero"))
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "grado_out", "ValueType")])
  }
}

jpeg("LAND_grado_OUT_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
  ggplot() +
    geom_point(data = subset(all_points, ValueType == "Zero"), 
               aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "grey80") +
    geom_point(data = subset(all_points, ValueType == "Positive"), 
               aes(x = CENTROID_X, y = CENTROID_Y, color = grado_out), size = 0.7, shape = 15) +
    scale_color_gradientn(colors = colors, name = "Grado out", guide = guide_colorbar(barwidth = 5, barheight = 50)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    labs(title = "LAND grado OUT Q0", x = "Centroid X", y = "Centroid Y") +
    geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()



##### GGPLOT LAND grado OUT Q5 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()

for (i in id) {
  matrix_name <- paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_", i)
  if (exists(matrix_name)) {
    B <- get(matrix_name)
    B_df <- as.data.frame(B)
    
    B_df$grado_out <- B_df$grado.out.D50.4
    B_df <- B_df[B_df$LAND > 0 & B_df$J.uk > 0, ]
    B_df <- B_df[is.finite(B_df$grado_out), ]
    B_df <- B_df[!duplicated(B_df[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID")]), ]
    
    B_df$ValueType <- ifelse(B_df$grado_out > 0, "Positive", ifelse(B_df$grado_out < 0, "Negative", "Zero"))
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "grado_out", "ValueType")])
  }
}

jpeg("LAND_grado_OUT_Q5.jpeg", width = 8000, height = 8000, res = 300)
print(
  ggplot() +
    geom_point(data = subset(all_points, ValueType == "Zero"), 
               aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "grey80") +
    geom_point(data = subset(all_points, ValueType == "Positive"), 
               aes(x = CENTROID_X, y = CENTROID_Y, color = grado_out), size = 0.7, shape = 15) +
    scale_color_gradientn(colors = colors, name = "Grado out", guide = guide_colorbar(barwidth = 5, barheight = 50)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    labs(title = "LAND grado OUT Q5", x = "Centroid X", y = "Centroid Y") +
    geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()




##### GGPLOT LAND grado IN Q0 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()

for (i in id) {
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.LAND.ID_", i)
  if (exists(matrix_name)) {
    B <- get(matrix_name)
    B_df <- as.data.frame(B)
    
    B_df$grado_in <- B_df$grado.in.D50.4
    B_df <- B_df[B_df$LAND > 0 & B_df$J.uk > 0, ]
    B_df <- B_df[is.finite(B_df$grado_in), ]
    B_df <- B_df[!duplicated(B_df[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID")]), ]
    
    B_df$ValueType <- ifelse(B_df$grado_in > 0, "Positive", ifelse(B_df$grado_in < 0, "Negative", "Zero"))
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "grado_in", "ValueType")])
  }
}

jpeg("LAND_grado_IN_Q0.jpeg", width = 8000, height = 8000, res = 300)
print(
  ggplot() +
    geom_point(data = subset(all_points, ValueType == "Zero"), 
               aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "grey80") +
    geom_point(data = subset(all_points, ValueType == "Positive"), 
               aes(x = CENTROID_X, y = CENTROID_Y, color = grado_in), size = 0.7, shape = 15) +
    scale_color_gradientn(colors = colors, name = "Grado in", guide = guide_colorbar(barwidth = 5, barheight = 50)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    labs(title = "LAND grado IN Q0", x = "Centroid X", y = "Centroid Y") +
    geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()



##### GGPLOT LAND grado IN Q5 #####
id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)
all_points <- data.frame()

for (i in id) {
  matrix_name <- paste0("M.final.SA.LP.Q5.D50.4.LAND.ID_", i)
  if (exists(matrix_name)) {
    B <- get(matrix_name)
    B_df <- as.data.frame(B)
    
    B_df$grado_in <- B_df$grado.in.D50.4
    B_df <- B_df[B_df$LAND > 0 & B_df$J.uk > 0, ]
    B_df <- B_df[is.finite(B_df$grado_in), ]
    B_df <- B_df[!duplicated(B_df[, c("CENTROID_X", "CENTROID_Y", "EcoR_ID")]), ]
    
    B_df$ValueType <- ifelse(B_df$grado_in > 0, "Positive", ifelse(B_df$grado_in < 0, "Negative", "Zero"))
    
    all_points <- rbind(all_points, B_df[, c("CENTROID_X", "CENTROID_Y", "grado_in", "ValueType")])
  }
}

jpeg("LAND_grado_IN_Q5.jpeg", width = 8000, height = 8000, res = 300)
print(
  ggplot() +
    geom_point(data = subset(all_points, ValueType == "Zero"), 
               aes(x = CENTROID_X, y = CENTROID_Y), size = 0.7, shape = 15, color = "grey80") +
    geom_point(data = subset(all_points, ValueType == "Positive"), 
               aes(x = CENTROID_X, y = CENTROID_Y, color = grado_in), size = 0.7, shape = 15) +
    scale_color_gradientn(colors = colors, name = "Grado in", guide = guide_colorbar(barwidth = 5, barheight = 50)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    labs(title = "LAND grado IN Q5", x = "Centroid X", y = "Centroid Y") +
    geom_sf(data = limites_ecoregions_sf, color = "black", size = 0.1, fill = NA)
)
dev.off()


##### GGPLOT PERM grado IN Q0 #####

id <- 301:353 ; id <- setdiff(id, 351) 

colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)

all_background <- data.frame()
water_background <- data.frame()
all_points <- data.frame()

for (i in id) {
  
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  map_name    <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  
  if (exists(matrix_name) && exists(map_name)) {
    
    B <- get(matrix_name)
    map_data <- get(map_name)
    
    B_df <- as.data.frame(B)
    map_data <- as.data.frame(map_data)
    
    B_df$grado_in <- B_df$grado.in.D50.4
    
    # fondo general: todas las celdas de la ecorregión
    all_background <- rbind(
      all_background,
      map_data[, c("CENTROID_X", "CENTROID_Y")]
    )
    
    # fondo de agua dulce: usar J.uk > 0, no MAR > 0
    # porque en PERM J.uk es el que está indicando celdas con hábitat efectivo
    water_background <- rbind(
      water_background,
      map_data[map_data$J.uk > 0, c("CENTROID_X", "CENTROID_Y")]
    )
    
    # capa de color: solo celdas con agua y grado positivo
    B_df <- B_df[
      B_df$J.uk > 0 &
        is.finite(B_df$grado_in) &
        B_df$grado_in > 0,
    ]
    
    all_points <- rbind(
      all_points,
      B_df[, c("CENTROID_X", "CENTROID_Y", "grado_in")]
    )
  }
}

jpeg("00_PERM_grado_IN_Q0.jpeg", width = 8000, height = 8000, res = 300)

print(
  ggplot() +
    
    # fondo continental / ecorregiones en gris claro
    geom_point(
      data = all_background,
      aes(x = CENTROID_X, y = CENTROID_Y),
      color = "grey90",
      size = 0.7,
      shape = 15,
      alpha = 0.8
    ) +
    
    # todas las celdas con agua dulce en gris oscuro
    geom_point(
      data = water_background,
      aes(x = CENTROID_X, y = CENTROID_Y),
      color = "#758283",
      size = 0.7,
      shape = 15
    ) +
    
    # arriba: grado in positivo en gradiente
    geom_point(
      data = all_points,
      aes(x = CENTROID_X, y = CENTROID_Y, color = grado_in),
      size = 0.7,
      shape = 15
    ) +
    
    scale_color_gradientn(
      colors = colors,
      name = "Grado in",
      guide = guide_colorbar(barwidth = 5, barheight = 50)
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    
    labs(
      title = "PERM grado IN Q0",
      x = "Centroid X",
      y = "Centroid Y"
    ) +
    
    geom_sf(
      data = limites_ecoregions_sf,
      color = "black",
      size = 0.1,
      fill = NA
    )
)

dev.off()
##### GGPLOT PERM grado IN Q5 #####

id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)

all_background <- data.frame()
water_background <- data.frame()
all_points <- data.frame()

for (i in id) {
  
  matrix_name <- paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_", i)
  map_name    <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  
  if (exists(matrix_name) && exists(map_name)) {
    
    B_df <- as.data.frame(get(matrix_name))
    map_data <- as.data.frame(get(map_name))
    
    B_df$grado_in <- B_df$grado.in.D50.4
    
    all_background <- rbind(
      all_background,
      map_data[, c("CENTROID_X", "CENTROID_Y")]
    )
    
    water_background <- rbind(
      water_background,
      map_data[map_data$J.uk > 0, c("CENTROID_X", "CENTROID_Y")]
    )
    
    B_df <- B_df[
      B_df$J.uk > 0 &
        is.finite(B_df$grado_in) &
        B_df$grado_in > 0,
    ]
    
    all_points <- rbind(
      all_points,
      B_df[, c("CENTROID_X", "CENTROID_Y", "grado_in")]
    )
  }
}

jpeg("00_PERM_grado_IN_Q5.jpeg", width = 8000, height = 8000, res = 300)

print(
  ggplot() +
    geom_point(
      data = all_background,
      aes(x = CENTROID_X, y = CENTROID_Y),
      color = "grey90",
      size = 0.7,
      shape = 15,
      alpha = 0.8
    ) +
    geom_point(
      data = water_background,
      aes(x = CENTROID_X, y = CENTROID_Y),
      color = "#758283",
      size = 0.7,
      shape = 15
    ) +
    geom_point(
      data = all_points,
      aes(x = CENTROID_X, y = CENTROID_Y, color = grado_in),
      size = 0.7,
      shape = 15
    ) +
    scale_color_gradientn(
      colors = colors,
      name = "Grado in",
      guide = guide_colorbar(barwidth = 5, barheight = 50)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    labs(
      title = "PERM grado IN Q5",
      x = "Centroid X",
      y = "Centroid Y"
    ) +
    geom_sf(
      data = limites_ecoregions_sf,
      color = "black",
      size = 0.1,
      fill = NA
    )
)

dev.off()



##### GGPLOT PERM grado OUT Q0 #####

id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)

all_background <- data.frame()
water_background <- data.frame()
all_points <- data.frame()

for (i in id) {
  
  matrix_name <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  map_name    <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  
  if (exists(matrix_name) && exists(map_name)) {
    
    B_df <- as.data.frame(get(matrix_name))
    map_data <- as.data.frame(get(map_name))
    
    B_df$grado_out <- B_df$grado.out.D50.4
    
    all_background <- rbind(
      all_background,
      map_data[, c("CENTROID_X", "CENTROID_Y")]
    )
    
    water_background <- rbind(
      water_background,
      map_data[map_data$J.uk > 0, c("CENTROID_X", "CENTROID_Y")]
    )
    
    B_df <- B_df[
      B_df$J.uk > 0 &
        is.finite(B_df$grado_out) &
        B_df$grado_out > 0,
    ]
    
    all_points <- rbind(
      all_points,
      B_df[, c("CENTROID_X", "CENTROID_Y", "grado_out")]
    )
  }
}

jpeg("00_PERM_grado_OUT_Q0.jpeg", width = 8000, height = 8000, res = 300)

print(
  ggplot() +
    geom_point(
      data = all_background,
      aes(x = CENTROID_X, y = CENTROID_Y),
      color = "grey90",
      size = 0.7,
      shape = 15,
      alpha = 0.8
    ) +
    geom_point(
      data = water_background,
      aes(x = CENTROID_X, y = CENTROID_Y),
      color = "#758283",
      size = 0.7,
      shape = 15
    ) +
    geom_point(
      data = all_points,
      aes(x = CENTROID_X, y = CENTROID_Y, color = grado_out),
      size = 0.7,
      shape = 15
    ) +
    scale_color_gradientn(
      colors = colors,
      name = "Grado out",
      guide = guide_colorbar(barwidth = 5, barheight = 50)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    labs(
      title = "PERM grado OUT Q0",
      x = "Centroid X",
      y = "Centroid Y"
    ) +
    geom_sf(
      data = limites_ecoregions_sf,
      color = "black",
      size = 0.1,
      fill = NA
    )
)

dev.off()



##### GGPLOT PERM grado OUT Q5 #####

id <- 301:353 ; id <- setdiff(id, 351) 
colors <- colorRampPalette(c("#fdf7c3", "#fde98a", "#fdd55f", "#fc9c46", "#f06642", "#c71f1f"))(100)

all_background <- data.frame()
water_background <- data.frame()
all_points <- data.frame()

for (i in id) {
  
  matrix_name <- paste0("M.final.SA.LP.Q5.D50.4.PERM.ID_", i)
  map_name    <- paste0("M.final.SA.LP.Q0.D50.4.PERM.ID_", i)
  
  if (exists(matrix_name) && exists(map_name)) {
    
    B_df <- as.data.frame(get(matrix_name))
    map_data <- as.data.frame(get(map_name))
    
    B_df$grado_out <- B_df$grado.out.D50.4
    
    all_background <- rbind(
      all_background,
      map_data[, c("CENTROID_X", "CENTROID_Y")]
    )
    
    water_background <- rbind(
      water_background,
      map_data[map_data$J.uk > 0, c("CENTROID_X", "CENTROID_Y")]
    )
    
    B_df <- B_df[
      B_df$J.uk > 0 &
        is.finite(B_df$grado_out) &
        B_df$grado_out > 0,
    ]
    
    all_points <- rbind(
      all_points,
      B_df[, c("CENTROID_X", "CENTROID_Y", "grado_out")]
    )
  }
}

jpeg("00_PERM_grado_OUT_Q5.jpeg", width = 8000, height = 8000, res = 300)

print(
  ggplot() +
    geom_point(
      data = all_background,
      aes(x = CENTROID_X, y = CENTROID_Y),
      color = "grey90",
      size = 0.7,
      shape = 15,
      alpha = 0.8
    ) +
    geom_point(
      data = water_background,
      aes(x = CENTROID_X, y = CENTROID_Y),
      color = "#758283",
      size = 0.7,
      shape = 15
    ) +
    geom_point(
      data = all_points,
      aes(x = CENTROID_X, y = CENTROID_Y, color = grado_out),
      size = 0.7,
      shape = 15
    ) +
    scale_color_gradientn(
      colors = colors,
      name = "Grado out",
      guide = guide_colorbar(barwidth = 5, barheight = 50)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 60, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 20)
    ) +
    labs(
      title = "PERM grado OUT Q5",
      x = "Centroid X",
      y = "Centroid Y"
    ) +
    geom_sf(
      data = limites_ecoregions_sf,
      color = "black",
      size = 0.1,
      fill = NA
    )
)

dev.off()