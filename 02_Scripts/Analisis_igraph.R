# Análisis de red con igraph

#Siguiendo la metodología de Álvarez (2024).

## Medidas de centralidad

# Degree centrality
# para determinar cuales son los primeros cinco hubs 
hubs <- sort(degree(fg), decreasing = TRUE)[1:5] # toma los cinco hubs principales 
hubs
#para visualizar los hubs 
plot(fg)

#Closeness centrality
cc <- sort(closeness(fg), decreasing = TRUE)[1:5]
cc

#betweenness centrality 
bc <- betweenness(fg, normalized = TRUE)[1:5]
bc

# detección de comunidades 
cfg <- cluster_walktrap(fg, steps = 4) #determina los "pasos" para que ueda determinar los cluster

plot(cfg, fg, layout = layout_with_fr(fg), #el layout = layout_with_fr(fg) es para ubicar los nodos y que no se vayan a encimas uno sobre otro
     vertex.size = 12, vertex.label.cex = 0.5, # tamaño del nodo y cambia el tamaño de las letras 
     main = paste0("random Walk (steps=4)\nQ = ", round(modularity(fg, membership(cfg)), 4), #\n sirve como "salto de línea" (como presionar Enter). Hace que el título no quede largo hacia los lados, sino que se divida en dos renglones.
                   # round(modularity(fg, membership(cfg)), 4) mide que tan buena fue la division de los cluster y redondea los numeros a 4 dígitos
                   " | k = ", length(cfg))) # numero de cluster por la red 



