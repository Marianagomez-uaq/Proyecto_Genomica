# Asignación taxonómica de hongos con DADA2

# install.packages("BiocManager")
# BiocManager::install("dada2")

library(dada2)

# Siguiendo la metodología de Álvarez (2024).

path_h <- "C:/Users/anabe/Documents/mariana/genomica/datos_proyecto/Hongos" # Esta línea funciona únicamente en mi computadora, para correrla en otra máquina, reemplazar la dirección por la local.
# Se usa un path fuera del proyecto debido al peso de los archivos, los cuales están disponibles para descarga en 01_RawData/datos_secuenciacion
list.files (path_h) # Para asegurarnos que los archivos que necesitaremos sí están en la carpeta


fnFs_h <- sort(list.files(path_h, pattern="_R1.fastq", full.names = TRUE)) # Genera una lista ordenada de las direcciones completas de los archivos forward (R1) con extensión .fastq (incluye los .fastq y los fastq.gz)
fnRs_h <- sort(list.files(path_h, pattern="_R2.fastq", full.names = TRUE)) # Reverse (R2)


sample.names_h <- sapply(strsplit(basename(fnFs_h), "_"), `[`, 1) # Toma los nombres de las muestras (como R1 y R2 son de las mismas muestras, solo es necesario tomar en cuenta uno)
# basename quita la ruta y solo deja el nombre del archivo
# "_" indica que se separe el nombre en donde está el guión bajo, y pone las dos partes en un vector
# '[', 1 indica que tome el primer elemento del vector generado (el nombre)


pdf("03_Results/Quality_Forward_Hongo.pdf",width=13,height = 8) # genera un pdf con cuadros de  13x8  de la calidad Phred de las lecturas de los primeros 10 archivos fastq
plotQualityProfile(fnFs_h[1:10])
dev.off()
pdf("03_Results/Quality_Reverse_Hongo.pdf",width=13,height = 8)
plotQualityProfile(fnRs_h[1:10])
dev.off()


filtFs_h <- file.path(path_h, "filtered", paste0(sample.names_h, "_F_filt.fastq.gz"))
# Este código no filtra todavía, simplemente genera una ruta a un subdirectorio del path (que no se ha creado), y le agrega a cada nombre la terminación "_F_filt.fastq.gz".
filtRs_h <- file.path(path_h, "filtered", paste0(sample.names_h, "_R_filt.fastq.gz"))
names(filtFs_h) <- sample.names_h # indica cual es el nombre de la muestra de cada ruta
names(filtRs_h) <- sample.names_h


# Debido a que el marcador fue ITS, no es recomendable usar el argumento truncLen, pero encontré esta información después de haber corrido todo, entonces mantuve el argumento. No se recomienda debido a la variabilidad de la longitud de esta secuencia (DADA2 Pipeline Tutorial (1.16), n.d.).
out_h <- filterAndTrim(fnFs_h, filtFs_h, fnRs_h, filtRs_h, truncLen=c(235,185), #Se seleccionaron los números de truncLen basado en el phred score de los gráficos de calidad
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, # maxN (máximo de bases ambiguas), maxEE (máximo error esperado), truncQ (elimina lecturas con calidad menor a 2), rm.phix (elimina bases parecidas al genoma del fago phiX (control illumina))
                       compress=TRUE, multithread=F)


errF_h <- learnErrors(filtFs, multithread=F) # Estima la probabilidad de que haya transiciones de nucleótidos por la secuenciación. Permite diferenciar entre variación biológica y errores de secuenciación (Benjjneb, 2025a).
saveRDS (errF_h,file="03_Results/errF_hongo.RDS") # Guarda el archivo en un RDS, es como un checkpoint que podemos utilizar sin tener que correr todo de nuevo al clonar el repositorio.

errR_h <- learnErrors(filtRs, multithread=F)
saveRDS (errR_h,file="03_Results/errR_hongo.RDS")

png("03_Results/Errores_Forward_Hongo.png") # genera imagen .png para visualizar si el modelo de error que generó DADA2 se ajusta a los datos de secuenciación (Benjjneb, 2025c).
plotErrors(errF_h, nominalQ=TRUE) # nominalQ son los scores de calidad
dev.off()

png("03_Results/Errores_Reverse_Hongo.png")
plotErrors(errR_h, nominalQ=TRUE)
dev.off()


##################### pendiente de correr:



dadaFs_h <- dada(filtFs_h, err=errF_h, multithread=F) # La función dada hace la "inferencia de muestras", que usa las secuencias filtradas y el modelo de error para generar ASV (amplicon sequence variants), que serán procesadas para poder hacer la asignación taxonómica (Benjjneb, 2025b)
dadaRs_h <- dada(filtRs_h, err=errR_h, multithread=F)


mergers_h <- mergePairs(dadaFs_h, filtFs_h, dadaRs_h, filtRs_h, verbose=TRUE) 
saveRDS(mergers_h, file="03_Results/mergers_hongo.RDS")


seqtab_h <- makeSequenceTable(mergers_h) 
dim(seqtab_h)
saveRDS(seqtab_h, file="03_Results/seqtab_hongo.RDS")


seqtab.nochim_h <- removeBimeraDenovo(seqtab_h, method="consensus", multithread=F, verbose=TRUE)
sum(seqtab.nochim_h)/sum(seqtab_h) # 0.944527


getN <- function(x) sum(getUniques(x))
track_h <- cbind(out_h, sapply(dadaFs_h, getN), sapply(dadaRs_h, getN), sapply(mergers_h, getN), rowSums(seqtab.nochim_h))


colnames(track_h) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_h) <- sample.names_h
head(track_h)


# Asignación utilizando base de datos SILVA (Quast et al., 2013; Yilmaz et al., 2014; Callahan, 2024).
taxa_h <- assignTaxonomy(seqtab.nochim_h, "C:/Users/anabe/Documents/mariana/genomica/datos_proyecto/Silva/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=F) # llega a nivel de género
saveRDS(taxa_h, file="03_Results/taxa_hongo.RDS")


taxa_h <- addSpecies(taxa_h, "C:/Users/anabe/Documents/mariana/genomica/datos_proyecto/Silva/silva_v138.2_assignSpecies.fa.gz") # llega a nivel de especie
saveRDS(taxa_h, file="03_Results/taxa_hongo.RDS")


taxa.print_h <- taxa_h # Se hace una copia del archivo taxa, para eliminar los nombres de las secuencias, solo para visualizar
rownames(taxa.print_h) <- NULL # quita rownames
saveRDS(taxa.print_h, file="03_Results/taxa.print_hongo.RDS")