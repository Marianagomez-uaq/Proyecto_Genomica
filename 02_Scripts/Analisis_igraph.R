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

#betweness centrality 
bc <- betweenness(fg, normalized = TRUE)[1:5]
bc

# detección de comunidades 
cfg <- cluster_walktrap(fg, steps = 4)

plot(cfg, fg, layout = lay,
     vertex.size = 12, vertex.label.cex = 0.5,
     main = paste0("Walktrap (steps=4)\nQ = ", round(modularity(fg, membership(cfg)), 4),
                   " | k = ", length(cfg)))

plot_dendrogram(cfg, mode = "hclust", main = "Dendrograma")



