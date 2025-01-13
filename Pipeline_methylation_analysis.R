################################################################################
#   LIBRERÍAS #

# install.packages("BiocManager")
# BiocManager::install(version = "3.20")
# BiocManager::install("remotes")
# BiocManager::install("minfi")
# BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
# install.packages("RColorBrewer")
# BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
# BiocManager::install("matrixStats")
# BiocManager::install("limma")
# BiocManager::install("DMRcate")
# BiocManager::install("pheatmap")
# BiocManager::install("DMRcatedata")
# install.packages("ggplot2")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("ReactomePA")
# BiocManager::install("msigdbr")


library("limma") 
library("matrixStats")
library("minfi")
library("RColorBrewer")
library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
library("pheatmap")
library("DMRcate")
library("IlluminaHumanMethylationEPICv2manifest")
library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")
library("dplyr")
library("ReactomePA")
library("msigdbr")


dataDirectory <- "ruta_directorio_trabajo"
setwd(dataDirectory)




################################################################################
#   INTRODUCCIÓN DE LOS DATOS CRUDOS EN R   #

#Leer la SampleSheet con los metadatos
targets <- read.metharray.sheet(dataDirectory, pattern = "SampleSheet.csv")
targets

#Crea el RGChannelSet que contiene los datos de intensidad sin procesar
rgSet <- read.metharray.exp(targets=targets)

#La anotación del array es la del EPICv2, hay que incluirla 
annotation <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

#Comprobación de los datos 
rgSet
pData(rgSet) 
annotation(rgSet)
dim(rgSet)

#Cambiamos los nombres de las muestras para que aparezca también si son controles o muestras de pacientes junto con el ID
targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep = ".")
sampleNames(rgSet) <- targets$ID
pData(rgSet)


#Convertimos el RGChannelSet en un MethylSet para tener los valores de beta previos a la normalización
MSet <- preprocessRaw(rgSet)
MSet
dim(MSet)
phenoData <- pData(MSet)



################################################################################
#   CONTROL DE CALIDAD   #

#1#

#QCReport: interpretación en la vignette de minfi. 

qcReport(
  rgSet, 
  sampNames = targets$ID, 
  sampGroups = targets$Sample_Group,
  pdf = "qcReport.pdf"
)

#2#

#Control de los valores de p de detección

options(matrixStats.useNames.NA = "deprecated")
detP <- detectionP(rgSet)
head(detP)


################################################################################
#   NORMALIZACIÓN   #

# Normalización con preprocessFunnorm, el output es un GenomicRatioSet
GenRset <- preprocessFunnorm(rgSet)
dim(GenRset)

#Obtenemos los valores de beta y m después de la normalización
beta_after <- getBeta(GenRset)
m_after <- getM(GenRset)


################################################################################
#   FILTRADO DE SONDAS   #

#1. Filtramos las sondas que fallan en una o mas muestras 
detP_filter <- detP[match(featureNames(GenRset), rownames(detP)),]

keep_probes <- rowSums(detP_filter < 0.01) == ncol(GenRset) 
table(keep_probes)

GenRset_filter <- GenRset[keep_probes,]
GenRset_filter
dim(GenRset_filter)

#2. Filtramos las sondas de los cromosomas X e Y
keep_X_Y <- !(featureNames(GenRset_filter) %in% annotation$Name[annotation$chr %in% 
                                                                  c("chrX","chrY")])
table(keep_X_Y)
GenRset_filter_1_2 <- GenRset_filter[keep_X_Y,]
GenRset_filter_1_2
dim(GenRset_filter_1_2)

m_normalized_filtered_1_2 <- getM(GenRset_filter_1_2)
beta_normalized_filtered_1_2 <- getBeta(GenRset_filter_1_2)

#3. Filtramos las sondas que hibridan en regiones con SNPs y las que tienen reacciones cruzadas

m_norm_filt_DMRcate <- rmSNPandCH(m_normalized_filtered_1_2, rmcrosshyb = TRUE)
beta_norm_filt_DMRcate <- rmSNPandCH(beta_normalized_filtered_1_2, rmcrosshyb = TRUE)

dim(m_norm_filt_DMRcate)
dim(beta_norm_filt_DMRcate)

################################################################################

#Comparación de los valores de beta sin procesar, después de la normalización y después del filtrado:

#png("beta_comparison_all.png", width = 1200, height = 600)
#par(mfrow=c(1,3))
densityPlot(MSet, sampGroups = phenoData$Sample_Group, main = "Beta Values Before Normalization")
densityPlot(beta_after, sampGroups = phenoData$Sample_Group, main = "Beta Values After Normalization")
densityPlot(beta_norm_filt_DMRcate, sampGroups = phenoData$Sample_Group, main = "Beta Values After Filtering")


################################################################################
#   ANÁLISIS DE ESCALAMIENTO MULTIDIMENSIONAL #


#MDS en función de los grupos: KBG o Control 
#png("MDS_KBG_normal_filter_normalized_control_all.png", width = 800, height = 600)
pal <- c("grey", "red")
plotMDS(m_norm_filt_DMRcate, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], pch=16, cex=2)
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)
par(cex.axis = 1.5)  
par(cex.lab = 1.8)

#MDS dependiendo del array
#colors <- colorRampPalette(brewer.pal(8, "Set3"))(65)
#png("MDS_Array_filter_normalized_control_all.png", width = 800, height = 600)

#plotMDS(m_norm_filt_DMRcate,gene.selection="common",  
#        col=colors[factor(targets$Slide)], pch=16)
#legend("top", legend=levels(factor(targets$Slide)), text.col=pal,
#       bg="white", cex=0.7)

#MDS en función de los grupos: KBG o Control (separados por edad)
pal2 <- c("grey", "red", "grey", "grey", "grey") #cambiar el color en función del grupo de edad

plotMDS(m_norm_filt_DMRcate, gene.selection="common", 
        col=pal2[factor(targets$Pool_ID)], pch=16, cex=2)
legend("topleft", legend=levels(factor(targets$Pool_ID)), text.col=pal2,
       bg="white", cex=0.7)
par(cex.axis = 1.5)  
par(cex.lab = 1.8)


################################################################################
#   SONDAS DIFERENCIALMENTE METILADAS (DMPs) #

#Los valores de m son mejores para el análisis estadístico
dim(m_norm_filt_DMRcate)

#Para la matriz de diseño cogemos el sample_Group para incluir el factor de la edad ponemos el Pool_ID (contiene las categorías de edad) 

Sample_Group <- factor(targets$Sample_Group, levels = c("CONTROL", "KBG"))
Edad <- factor(targets$Pool_ID)


design_matrix_1 <- model.matrix(~0 + Sample_Group + Edad, data = targets)
colnames(design_matrix_1) <- c(levels(Sample_Group), levels(Edad) [-1])

#Creamos el modelo lineal

fit_1 <- lmFit(m_norm_filt_DMRcate, design_matrix_1)

contMatrix_1 <- makeContrasts(CONTROL - KBG, levels=design_matrix_1)
contMatrix_1


#Cálculo de los p-values de la expresión diferencial y ranking de genes 
fit_contMatrix1 <- contrasts.fit(fit_1, contMatrix_1)


fit_bayes_1 <- eBayes(fit_contMatrix1)
DMPs_1 <- topTable(fit_bayes_1, num=Inf, coef=1)

head(DMPs_1)

#Añadimos la anotación de las sondas a las DMPs
annotation_m <- annotation[match(rownames(m_norm_filt_DMRcate), annotation$Name), c(1:4,12:19,24:ncol(annotation))]
DMPs_info_1 <- topTable(fit_bayes_1, num=Inf, coef=1, genelist=annotation_m)
head(DMPs_info_1)
write.table(DMPs_info_1, file="DMPs_info_1.csv", sep=",", row.names = FALSE)


#Calculamos el EpiScore
DMPs_info_1$logP <- -log10(DMPs_info_1$adj.P.Val)

DMPs_info_1$EpiScore <- abs(DMPs_info_1$logFC) * DMPs_info_1$logP

#Filtramos las DMPs en base al valor p ajustado y al logFC
DMPs_filtered_1 <- DMPs_info_1[DMPs_info_1$adj.P.Val < 0.05 & abs(DMPs_info_1$logFC) > 1, ]

#Ordenamos las DMPs filtradas en base al EpiScore
DMPs_filtered_1 <- DMPs_filtered_1[order(-DMPs_filtered_1$EpiScore), ]

write.csv(DMPs_filtered_1, file = "DMPs_filtered_1.csv", row.names = FALSE)

#### VOLCANO PLOT

DMPs_1$logP <- -log10(DMPs_1$adj.P.Val)
logFC_threshold <- 1 
p_adj_value_threshold <- 0.05

DMPs_1$Significancia <- with(DMPs_1, ifelse( adj.P.Val < p_adj_value_threshold & abs(logFC) > logFC_threshold, "Significativo", "No significativo") )



#png("volcano_plot.png", width = 1000, height = 600)

ggplot(DMPs_1, aes(x=logFC, y = logP, color = Significancia)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Significativo" = "red", "No Significativo" = "grey")) +
  theme_minimal() +
  labs(title = "Sondas Diferencialmente Metiladas",  x = "Log2 Fold-Change", y = "Log 10 P-value ajustado", color = "Significancia")+
  theme(
    axis.title = element_text(size = 14),   
    axis.text = element_text(size = 12),    
    plot.title = element_text(size = 16, hjust = 0.5) 
  )


#### BOXPLOT DE LA SONDA CON MAYOR EPISCORE

#Visualización de las sondas diferencialmente metiladas (se hace con beta 
#porque es mejor para visualizar)

head(DMPs_filtered_1)
targets$Sample_Group <- as.factor(targets$Sample_Group)

#png("top_cpgs.png", width = 800, height = 600)
plotCpg(beta_norm_filt_DMRcate, cpg="cg10951305_TC11", pheno=targets$Sample_Group, ylab="Beta values")


#### HEATMAP DE LAS 100 SONDAS CON MAYOR EPISCORE
# Subset de las 100 DMPs
DMPs_top100 <- head(DMPs_filtered_1, 100)

names_top_100_DMPs <- rownames(head(DMPs_top100, 100))

beta_matrix_significant <- beta_norm_filt_DMRcate[names_top_100_DMPs, ]

# Seleccionar los KBGs y 30 controles para la gráfica
control_cols <- grep("^CONTROL", colnames(beta_norm_filt_DMRcate), value = TRUE)
kbg_cols <- grep("^KBG", colnames(beta_norm_filt_DMRcate), value = TRUE)
selected_cols <- c(control_cols[1:30], kbg_cols)

beta_subset <- beta_matrix_significant[, selected_cols]

#png("pheatmap_peq.png", width = 1200, height = 600)
pheatmap(beta_subset, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, 
         main = "Heatmap de top 100 DMPs")


################################################################################
#   REGIONES DIFERENCIALMENTE METILADAS (DMRs)   #

#Partimos de:
design_1 <- model.matrix(~ Sample_Group + Edad, data = targets)


#Anotamos con la información cromosómica de las sondas y colapsamos las sondas duplicadas en su media

myannotation <- cpg.annotate(datatype = "array", object = m_norm_filt_DMRcate,
                             what = "M", arraytype = "EPICv2", epicv2Remap = F, 
                             analysis.type = "differential",
                             design = design_1, coef = 2, epicv2Filter = "mean", fdr=0.001)

#Identificación de las DMRs
dmrcoutput <- dmrcate(myannotation) 

#Obtenemos la información de las DMRs y la guardamos en formato csv
results.ranges <- extractRanges(dmrcoutput = dmrcoutput, genome = "hg38")
length(results.ranges)

results_df <- as.data.frame(results.ranges)
write.csv(results_df, file = "results_ranges.csv", row.names = FALSE)


#### gráfico DMRs por cromosoma, coloreado por tamaño de las DMRs
# Crear categorías de tamaño de región
results_df$region_size <- cut(results_df$width, breaks = c(0, 500, 1000, 5000, Inf),
                                     labels = c("0-500", "501-1000", "1001-5000", ">5000"))

# Barplot apilado
#png("DMR_result1_control_todas_edades_2.png", width = 1200, height = 600)
ggplot(results_df, aes(x = seqnames, fill = region_size)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  labs(title = "Distribución de Tamaño de Regiones por Cromosoma",
       x = "Cromosoma",
       y = "Número de Regiones",
       fill = "Tamaño de Región (nt)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### gráfico de DMRs normalizado en función del tamaño del cromosoma
# Crear el DataFrame vacío
DMRs_graph <- data.frame(
  Chromosome = character(),
  Total_DMR_nt = numeric(),
  Chromosome_Length = numeric(),
  Normalized_DMR_nt = numeric(),
  stringsAsFactors = FALSE
)

# Calcular el total de nt por cromosoma en las DMRs
dmr_nt_per_chromosome <- aggregate(width ~ seqnames, data = as.data.frame(results_df), sum)

# Agregar la longitud total de cada cromosoma
lengths_by_chromosome <- c(chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, 
                           chr4 = 190214555, chr5 = 181538259, chr6 = 170805979, 
                           chr7 = 159345973, chr8 = 145138636, chr9 = 138394717, 
                           chr10 = 133797422, chr11 = 135086622, chr12 = 133275309, 
                           chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, 
                           chr16 = 90338345, chr17 = 83257441, chr18 = 80373285, 
                           chr19 = 58617616, chr20 = 64444167, chr21 = 46709983, 
                           chr22 = 50818468)

# Convertir a un DataFrame para combinar con el resultado anterior
lengths_df <- data.frame(
  Chromosome = names(lengths_by_chromosome),
  Chromosome_Length = lengths_by_chromosome,
  stringsAsFactors = FALSE
)

# Renombrar columnas para hacer el merge
colnames(dmr_nt_per_chromosome) <- c("Chromosome", "Total_DMR_nt")

# Combinar los datos de las DMRs con las longitudes totales de los cromosomas
DMRs_graph <- merge(dmr_nt_per_chromosome, lengths_df, by = "Chromosome", all = TRUE)

# Calcular la normalización
DMRs_graph$Normalized_DMR_nt <- DMRs_graph$Total_DMR_nt / DMRs_graph$Chromosome_Length

# Guardar el DataFrame en un archivo CSV
write.csv(DMRs_graph, "DMRs_graph.csv", row.names = FALSE)

#gráfico de barras
#png("DMR_chr_2_control_.png", width = 1000, height = 600)
ggplot(DMRs_graph, aes(x = Chromosome, y = Normalized_DMR_nt)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Proporción Normalizada de nt en DMRs por Cromosoma",
    x = "Cromosoma",
    y = "Proporción Normalizada de nt en DMRs"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.title = element_text(hjust = 0.5, size = 16) 
  )
theme(
  axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
  axis.text.y = element_text(size = 12), 
  axis.title.x = element_text(size = 14, margin = margin(t = 10)), 
  axis.title.y = element_text(size = 14, margin = margin(r = 10)), 
  plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = 15)) 
)

### Pie chart 

# Calcular la distribución de las categorías de tamaño
size_distribution <- results_df %>%
  group_by(region_size) %>%
  summarise(count = n()) %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) # Calcular porcentajes

# Crear el gráfico
#png("pie_chart.png", width = 800, height = 600)
ggplot(size_distribution, aes(x = "", y = count, fill = region_size)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(
    title = "Distribución de Tamaño de Regiones (DMRs)",
    fill = "Tamaño de Región (nt)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = 10))
  ) +
  scale_fill_brewer(palette = "Set3")

################################################################################

#   ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL DE LOS GENES DE LAS DMRs   #

#Filtramos las DMRs en base a la diferencia de expresión media 

DMRs_filtered <- results.ranges[abs(results.ranges$meandiff) > 0.1, ]
length(DMRs_filtered)

#guardamos como archivo csv las DMRs filtradas
results_df_filter <- as.data.frame(DMRs_filtered)
write.csv(results_df_filter, file = "results_ranges_filtered.csv", row.names = FALSE)

#listado de genes de las regiones diferencialmente metiladas
DMGenes <- results_df_filter$overlapping.genes

DMGenes_list <- unlist(strsplit(as.character(DMGenes), split=","))

DMGenes_unique <- unique(DMGenes_list)

DMGenes_unique <- na.omit(DMGenes_unique) #elimina NA
length(DMGenes_unique)

DMGenes_unique <- trimws(DMGenes_unique)  # Eliminar espacios en blanco
DMGenes_unique <- toupper(DMGenes_unique)  # Convertir a mayúsculas si aplica
DMGenes_unique <- unique(DMGenes_unique)  # Eliminar duplicados

#Conversión de nombre de genes a Entrez IDs
entrez_genes <- bitr(DMGenes_unique, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
print(entrez_genes)

unmapped_genes <- setdiff(DMGenes_unique, entrez_genes$SYMBOL)

############Análisis de enriquecimiento funcional

###GO

#BP
go_enrichment_BP <- enrichGO(gene = entrez_genes$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "BP", # Ontologías disponibles: "BP", "MF", "CC"
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)

#png("go_enrichment_BP.png", width = 800, height = 600)
dotplot(go_enrichment_BP, showCategory = 10) + ggtitle("GO Enrichment (BP)")

write.csv(as.data.frame(go_enrichment_BP), "GO_enrichment_BP_results.csv", row.names = FALSE)


#CC
go_enrichment_CC <- enrichGO(gene = entrez_genes$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "CC", # Ontologías disponibles: "BP", "MF", "CC"
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)

#png("go_enrichment_CC.png", width = 800, height = 600)
dotplot(go_enrichment_CC, showCategory = 10) + ggtitle("GO Enrichment (CC)")

write.csv(as.data.frame(go_enrichment_CC), "GO_enrichment_CC_results.csv", row.names = FALSE)


#MF
go_enrichment_MF <- enrichGO(gene = entrez_genes$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "MF", # Ontologías disponibles: "BP", "MF", "CC"
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)

#png("go_enrichment_MF.png", width = 800, height = 600)
dotplot(go_enrichment_MF, showCategory = 10) + ggtitle("GO Enrichment (MF)")

write.csv(as.data.frame(go_enrichment_MF), "GO_enrichment_MF_results.csv", row.names = FALSE)



###KEGG

kegg_enrichment <- enrichKEGG(gene = entrez_genes$ENTREZID,
                              organism = 'hsa', # Para humanos
                              pvalueCutoff = 0.05, 
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2)

#png("kegg_enrichment_BP.png", width = 800, height = 600)
dotplot(kegg_enrichment, showCategory = 10) + ggtitle("KEGG Pathways Enrichment")

write.csv(as.data.frame(kegg_enrichment), "KEGG_enrichment_results.csv", row.names = FALSE)


###HPOs

#Descargar el archivo: phenotype to genes que hay en: https://hpo.jax.org/data/annotations
file_path<-"ruta_al _archivo"

phenotype_to_genes <- readr::read_delim(file_path, delim = "\t", col_names = TRUE)

#agrupar cada HPO con todos sus genes
unique_genes_per_hpo <- phenotype_to_genes %>%
  group_by(hpo_id, hpo_name) %>%
  summarize(unique_genes = list(unique(gene_symbol)), .groups = 'drop')

# Run enricher
HPO_enrichment <- enricher(
  gene = entrez_genes_2$SYMBOL,      
  TERM2GENE = phenotype_to_genes[,c(1,4)],        
  TERM2NAME = unique_genes_per_hpo[,c(1:2)],       
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1)

enrich_HPO_df <- as.data.frame(HPO_enrichment)
dotplot(HPO_enrichment,title = "Enriquecimiento Funcional de los HPOs")


### PANELES
load("ruta_archivo/paneles_genes_para_gsea.Rda")

unique_panels<-unique(paneles_genes_para_gsea$FileName)
unique_panels_df<-data.frame(panel_ID=c(1:length(unique_panels)),panel_name=unique_panels)

paneles_enrichment <- enricher(
  gene = entrez_genes_2$SYMBOL,      
  TERM2GENE = paneles_genes_para_gsea,        # Gene sets formatted for GSEA: paneles con genes (repetidos paneles, 1 por gen -> rollo 2 filas que tienen panel_DHR:ABCA4, panel_DHR=USH2A)
  TERM2NAME = unique_panels_df,pvalueCutoff = 1,qvalueCutoff = 1)       # HPO_ID and HPO_name: 2-column datafram)


enrich_df <- as.data.frame(paneles_enrichment)
dotplot(paneles_enrichment, title = "Enriquecimiento Funcional de los paneles")


###DO
do_results <- DOSE::enrichDO(gene = entrez_genes_2$ENTREZID,
                             organism = "hsa",
                             pAdjustMethod = "BH")


enrich_df_do <- as.data.frame(do_results)
dotplot(do_results,title = "Análisis de enriquecimiento funcional con Disease Ontology (DO)")


###REACTOME

reactome_results <- ReactomePA::enrichPathway(
  gene = entrez_genes_2$ENTREZID,     
  organism = "human",               
  pAdjustMethod = "BH",               
  pvalueCutoff = 1
)

enrich_df_reactome <- as.data.frame(reactome_results)


###MSigDB


msig_data <- msigdbr(species = "Homo sapiens", category = "H") # Hallmark gene sets

msig_results <- enricher(
  gene = entrez_genes_2$SYMBOL,
  TERM2GENE = msig_data[, c("gs_name", "gene_symbol")],
  pAdjustMethod = "BH"
)

enrich_df_mSigDB <- as.data.frame(msig_results)


