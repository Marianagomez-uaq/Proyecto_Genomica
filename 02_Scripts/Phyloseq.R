# Phyloseq bacteria

# install.packages("BiocManager")
# BiocManager::install(version = "3.22") # Para actualizar la versión de Bioconductor
BiocManager::install("phyloseq") # Phyloseq no es de CRAN, es de Bioconductor

library(phyloseq)

library(ggplot2)


# Crear archivo phyloseq
seqtab_b <- seqtab_bacteria
metadata_b <- read.csv("01_RawData/sample_data_bacteria.csv") # Lee el .csv
orden <- rownames(seqtab_b) # Esto se lo pedimos a la IA, crea un objeto llamado orden con los nombres de las secuencias
metadata_b$ID_sequencing <- as.character(metadata_b$ID_sequencing) # Convierte los ID de numéricos a caracter
metadata_ordenado <- metadata_b[match(orden, metadata_b$ID_sequencing), ] # Ordena los metadatos basados en el ID_sequencing, en el orden que aparecen en seqtab_bacteria

metadata_ordenado$ID_sequencing

ID <- rownames (seqtab_b) # genera un objeto con los ID
Continent <- metadata_ordenado$Continent # genera un objeto con las localizaciones de cada lectura

samdf_b <- data.frame (ID = ID, Continent = Continent) #genera un dataframe con los ID y el Continente al que pertenecen
samdf_b

rownames(samdf_b) <- rownames(seqtab_b) # le pusimos los nombres a las filas

###########################################################################################
# Distribución

# DADA2 produce: seqtab.nochim, taxa, samdf. Las siguientes líneas juntan todo esto en un objeto phyloseq (ps)
ps <- phyloseq( 
  otu_table(seqtab.nochim_b, taxa_are_rows = FALSE),
  tax_table(taxa_bacteria),
  sample_data(samdf_b)
)
ps

saveRDS(ps, file="03_Results/ps.RDS") # guarda el objeto ps en un archivo para poder subirlo al repositorio

sample_sums_vec <- sample_sums(ps) # esto nos permite identificar como se distribuyen las lecturas en cada muestra
distribucion <- summary(sample_sums_vec)

# min 0: indica que hubo al menos una muestra con 0 lecturas
# 1st Qu. 0: indica que hasta el primer cuartil (25%), ninguna muestra tenía lecturas
# Median (2nd Qu) 3: indica que hasta el segundo cuartil (50%), las muestras tienen máximo 3 lecturas
# 3rd Qu 13: indica que hasta el tercer cuartil (75%), las muestras tienen máximo 13 lecturas.
# Max 79: indica que hubo al menos una muestra con 79 lecturas, y ninguna muestra tuvo más
# Mean 8: en promedio las muestras tuvieron 8 lecturas, lo cual es muy malo, porque inicialmente se tenían más de 100 millones de lecturas. 

sum (sample_sums_vec) # Nosotros nos quedamos con 1907 lecturas. Y el autor original se quedó con 6,278,448 lecturas en este paso.
ntaxa(ps) # 731 ASV (el autor se quedó con 60,538 de bacterias)

df_reads <- data.frame( # genera un dataframe que junta: 
  Sample     = names(sample_sums_vec), # nombre de la muestra
  TotalReads = sample_sums_vec, # Total de lecturas que tuvo
  SampleType = sample_data(ps)$Continent # En qué continente se encontró
)
df_reads


png("03_Results/Distribucion_Bacteria.png") # gráfica para visualizar las distribuciones
ggplot(df_reads, aes(x = reorder(Sample, TotalReads), # acomoda las muestras según el total de lecturas
                     y = TotalReads, fill = Continent)) + # el eje y es el total de lecturas, y cada barra se colorea según el continente
  geom_bar(stat = "identity") + # identity deja los datos como son, no hace ninguna modificación
  coord_flip() + # hace que las barras sean horizontales y no verticales
  # scale_y_continuous(labels = comma) + # Esto no es necesario porque nuestras lecturas no llegan a los miles
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Total de Lecturas por Muestra", # título del gráfico
       x = "Muestra", y = "Numero de lecturas", fill = "Continente") + # nombre de los ejes y de la leyenda 
  theme_bw(base_size = 11) + # tamaño del texto
  theme(legend.position = "bottom") # indica donde se va a posicionar la leyenda
dev.off()

###################################################################################################
# Filtrado
taxa_are_rows(ps)
summary (taxa_sums(ps))

ps_filtered <- prune_taxa(taxa_sums(ps) > 1, ps) # Borra los ASVs con 1 lectura (para hacer un filtro de ejemplo)

# prevalence_threshold <- 0.25 * nsamples(ps_filtered) # calcula el 25% del total de muestras
# este es un filtro demasiado estricto considerando que tenemos muy pocas muestras, deja con cero muestras
# este filtro solo permite que se queden los que están presentes en más de 59 muestras, y eso no ocurre en ningún caso

prevalence_threshold <- 2 # Esto también remueve la mayoría (692), por lo que mejor no haremos filtrado por prevalencia 

taxa_prevalence      <- colSums(otu_table(ps_filtered) > 0) # verifica si el ASV está en la muestra (... > 0) y cuenta cuantas veces aparece cada ASV (rowsums)
ps_filtered <- prune_taxa(taxa_prevalence >= prevalence_threshold, ps_filtered) # selecciona solo los ASVs presentes en exactamente o más del 25% de las muestras, corta los que no cumplan la condición

cat("ASVs originales:", ntaxa(ps), "\n")
cat("ASVs tras filtrado:", ntaxa(ps_filtered), "\n")
cat("ASVs removidos:", ntaxa(ps) - ntaxa(ps_filtered), "\n")

saveRDS (ps_filtered, file="03_Results/ps_filtered.RDS") # después del filtrado nos quedamos con 692 ASV, mientras que el autor original se quedó con 4913

