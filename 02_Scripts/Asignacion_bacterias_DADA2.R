# Asignación taxonómica de bacterias con DADA2

# install.packages("BiocManager")
# BiocManager::install("dada2")

library(dada2)

# Siguiendo la metodología de Álvarez (2024).

path_b <- "C:/Users/anabe/Documents/mariana/genomica/datos_proyecto/Bacterias" # Esta línea funciona únicamente en mi computadora, para correrla en otra máquina, reemplazar la dirección por la local.
# Se usa un path fuera del proyecto debido al peso de los archivos, los cuales están disponibles para descarga en 01_RawData/datos_secuenciacion
list.files (path_b) # Para asegurarnos que los archivos que necesitaremos sí están en la carpeta


fnFs_b <- sort(list.files(path_b, pattern="_R1.fastq", full.names = TRUE)) # Genera una lista ordenada de las direcciones completas de los archivos forward (R1) con extensión .fastq (incluye los .fastq y los fastq.gz)
fnRs_b <- sort(list.files(path_b, pattern="_R2.fastq", full.names = TRUE)) # Reverse (R2)


sample.names_b <- sapply(strsplit(basename(fnFs_b), "_"), `[`, 1) # Toma los nombres de las muestras (como R1 y R2 son de las mismas muestras, solo es necesario tomar en cuenta uno)
# basename quita la ruta y solo deja el nombre del archivo
# "_" indica que se separe el nombre en donde está el guión bajo, y pone las dos partes en un vector
# '[', 1 indica que tome el primer elemento del vector generado (el nombre)


pdf("03_Results/Quality_Forward_Bacteria.pdf",width=13,height = 8) # genera un pdf con cuadros de  13x8  de la calidad Phred de las lecturas de los primeros 10 archivos fastq
plotQualityProfile(fnFs_b[1:10])
dev.off()
pdf("03_Results/Quality_Reverse_Bacteria.pdf",width=13,height = 8)
plotQualityProfile(fnRs_b[1:10])
dev.off()


filtFs_b <- file.path(path_b, "filtered", paste0(sample.names_b, "_F_filt.fastq.gz"))
# Este código no filtra todavía, simplemente genera una ruta a un subdirectorio del path (que no se ha creado), y le agrega a cada nombre la terminación "_F_filt.fastq.gz".
filtRs_b <- file.path(path_b, "filtered", paste0(sample.names_b, "_R_filt.fastq.gz"))
names(filtFs_b) <- sample.names_b # indica cual es el nombre de la muestra de cada ruta
names(filtRs_b) <- sample.names_b


out_b <- filterAndTrim(fnFs_b, filtFs_b, fnRs_b, filtRs_b, truncLen=c(230,180), #Se seleccionaron los números de truncLen basado en el phred score de los gráficos de calidad
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, # maxN (máximo de bases ambiguas), maxEE (máximo error esperado), truncQ (elimina lecturas con calidad menor a 2), rm.phix (elimina bases parecidas al genoma del fago phiX (control illumina))
                     compress=TRUE, multithread=F)


errF_b <- learnErrors(filtFs, multithread=F) # Estima la probabilidad de que haya transiciones de nucleótidos por la secuenciación. Permite diferenciar entre variación biológica y errores de secuenciación (Benjjneb, 2025a).
saveRDS (errF_b,file="03_Results/errF_bacteria.RDS") # Guarda el archivo en un RDS, es como un checkpoint que podemos utilizar sin tener que correr todo de nuevo al clonar el repositorio.
errR_b <- learnErrors(filtRs, multithread=F)
saveRDS (errR_b,file="03_Results/errR_bacteria.RDS")


png("03_Results/Errores_Forward_Bacteria.png") # genera imagen .png para visualizar si el modelo de error que generó DADA2 se ajusta a los datos de secuenciación (Benjjneb, 2025c).
plotErrors(errF_b, nominalQ=TRUE) # nominalQ son los scores de calidad
dev.off()
png("03_Results/Errores_Reverse_Bacteria.png")
plotErrors(errR_b, nominalQ=TRUE)
dev.off()


dadaFs_b <- dada(filtFs_b, err=errF_b, multithread=F) # La función dada hace la "inferencia de muestras", que usa las secuencias filtradas y el modelo de error para generar ASV (amplicon sequence variants), que serán procesadas para poder hacer la asignación taxonómica (Benjjneb, 2025cb
dadaRs_b <- dada(filtRs_b, err=errR_b, multithread=F) # borrar el otro y dejar este después de correrlo


mergers_b <- mergePairs(dadaFs_b, filtFs_b, dadaRs_b, filtRs_b, verbose=TRUE) 
saveRDS(mergers_b, file="03_Results/mergers_bacteria.RDS")


seqtab_b <- makeSequenceTable(mergers_b) 
dim(seqtab_b)
saveRDS(seqtab_b, file="03_Results/seqtab_bacteria.RDS")


seqtab.nochim_b <- removeBimeraDenovo(seqtab_b, method="consensus", multithread=F, verbose=TRUE)
sum(seqtab.nochim_b)/sum(seqtab_b) # 0.944527


getN <- function(x) sum(getUniques(x))
track_b <- cbind(out_b, sapply(dadaFs_b, getN), sapply(dadaRs_b, getN), sapply(mergers_b, getN), rowSums(seqtab.nochim_b))


colnames(track_b) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_b) <- sample.names_b
head(track_b)


# Asignación utilizando base de datos SILVA (Quast et al., 2013; Yilmaz et al., 2014; Callahan, 2024).
taxa_b <- assignTaxonomy(seqtab.nochim_b, "C:/Users/anabe/Documents/mariana/genomica/datos_proyecto/Silva/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=F) # llega a nivel de género
saveRDS(taxa_b, file="03_Results/taxa_bacteria.RDS")
# como los archivos silva son muy pesados, también los trabajé localmente, los links están disponibles para descarga en 01_RawData/datos_secuenciacion

taxa_b <- addSpecies(taxa_b, "C:/Users/anabe/Documents/mariana/genomica/datos_proyecto/Silva/silva_v138.2_assignSpecies.fa.gz") # llega a nivel de especie
saveRDS(taxa_b, file="03_Results/taxa_bacteria.RDS")


taxa.print_b <- taxa_b # Se hace una copia del archivo taxa, para eliminar los nombres de las secuencias, solo para visualizar
rownames(taxa.print_b) <- NULL # quita rownames
saveRDS(taxa.print_b, file="03_Results/taxa.print_bacteria.RDS")



#############PENDIENTE DE CORRER

library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
save(ps,file="03_Results/ps.RDS")


jpeg("03_Results/Alfa_Diversity.jpeg",width=610,height = 367)
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
dev.off()


# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
jpeg("03_Results/NMDS.jpeg")
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
dev.off()


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
jpeg("03_Results/BarPlot.jpeg")
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
dev.off()