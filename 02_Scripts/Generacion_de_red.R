# Generación de red a partir de datos filtrados

## instalar los paquetes que vamos a utilizar

install.packages("ggraph")
install.packages("corrr")
install.packages("tidygraph")

#cargar las librerias 
library(igraph)
library(corrr)
library(ggraph)
library(ggplot2)
library(tidygraph)

tabla_abundancia <- as.data.frame(t(otu_table(ps_filtered))) ## como ps_filtered esta en formato de phyloseq, lo pase a una matriz para que lo tome como un data.frame y poder hacer la red de correlación
#la (t) es para transponer las muestras en las columnas porque queremos una matriz de los ASV y los ID para la correlación 

# 1. Calcular la matriz de correlación
# 2. Enfocar (mantener solo las relaciones)
# 3. Convertir a un formato de red (nodos y enlaces)
red_psf <- tabla_abundancia %>% # %>% es un pipe, sirve para encadenar funciones 
  correlate(method = "pearson") %>% # Calcula que tanto se parecen los perfiles de abundancia entre cada par de bacterias. El método de Pearson busca relaciones lineales
  stretch(na.rm = TRUE) %>% # elimina los valores de NA en caso de que tenga 
  filter(abs(r) > 0.3) # Filtramos para mostrar solo correlaciones fuertes, correlación entre media y fuerte 

fg <- as_tbl_graph(red_psf, directed = FALSE) # hace una red dirigida
# as_tbl_graph es de la libreria tidygraph, 

ggraph(fg, layout = "graphopt") + # ggraph se utiliza en lugar de ggplot, ya que esta es una libreria que se utiliza para graficar redes 
  #fg es la base de datos 
  #layout = "graphopt" sirve para ayuda a ordenar la posición de los nodos basado en que tan fuerte es la interacción entre los nodos 
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) + # geom_edge_link sirve para colocar las conexiones entre nodos 
  #edge_alpha = abs(r) ayuda a que los nodos que tengan correlaciones cercanas a cero casi no se noten (sean más transparentes)
  #edge_width = abs(r) determina dependiendo de la fuerza de la correlación si la línea es más gruesa o más delgada   
  #color = (r) este determina el color de la línea, si el valor de la correlacion es positivo o negativo 
  geom_node_point(size = 4, color = "cyan3") + # Añadir los nodos, determina el tamaño y el color 
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold") + # Añade llos nombres de los nodos
  #label = name extrae los nombres 
  #repel = TRUE evita que se encimen las etiquqetas y los nodos 
  #fontface = "bold" pone el texto en negrita 
  scale_edge_color_gradient2(low = "pink", mid = "purple4", high = "red4", midpoint = 0) +  # Determina la escala de colores (Azul = Positiva, Rojo = Negativa)
  # midpoint = 0 hace que la transicion de los cambios de color sea a partir de 0 
  theme_graph() + # elimina el fondo y lo deja blanco 
  labs(title = "Red de Correlación", # pone el título de la imagen 
       edge_color = "Coef. Correlación (r)") # pone el nombre de la escala de colores
