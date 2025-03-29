library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(dplyr)

## PHYLOSEQ

data("dietswap", package = "microbiome")
ps <- dietswap
ps

# ¿Cuántas muestras y taxones contiene el objeto?
# ¿Qué variables están disponibles en los metadatos de las muestras?

ps_prune <- prune_taxa(taxa_sums(ps) > 1, ps)

View(sample_data(ps))
View(otu_table(ps))

# CURVAS DE RAREFACCION

otu.rare <- otu_table(ps_prune)
otu.rare <- as.data.frame(t(otu.rare))
sample_names <- rownames(otu.rare)

otu.rarecurve = rarecurve(otu.rare[1:20], step = 10000, label=TRUE, cex.lab=0.1)

otu.rarecurve = rarecurve(otu.rare[1:20], step = 10000,label=FALSE)

# DIVERSIDAD ALFA

plot_richness(ps_prune,x="bmi_group")

# FILTRADO

# DIVERSIDAD BETA

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="bmi_group", title="Bray NMDS")

# GRAFICA RANK ABUNDANCE

par(mar = c(10, 4, 4, 2) + 0.1)
barplot(sort(taxa_sums(ps.prop), TRUE)[1:40]/nsamples(ps.prop), las=2, cex.names = 0.63,
        main = "Curvas de Rango-Abundancia (40 Taxa más abundantes)",
        ylab = "Abundancia Relativa")

# Gráficas apiladas de abundancia por taxón

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="nationality", fill="Family") + facet_wrap(~timepoint, scales="free_x")

## GLOBAL PATTERNS

data("GlobalPatterns")
gp <- GlobalPatterns
gp 
View(sample_data(gp))
View(otu_table(gp))
View(tax_table(gp))

# Preprocesamiento
# 1. Filtrar taxa con menos de 5 lecturas en al menos 20% de las muestras
tax <- c()
for (i in 1:(dim(otu_table(gp))[1])){
  x <- 0
  for (j in 1:(dim(otu_table(gp))[2])){
    if (otu_table(gp)[i,j]<5){
      x <- x + 1
    }
  }
  if ((x/(dim(otu_table(gp))[2]))>=0.20){
    tax <- c(tax, rownames(otu_table(gp))[i])
  }
}

gp.taxfilt <- prune_taxa(tax, gp)
taxa_names(gp.taxfilt)

gp.taxfilt
View(otu_table(gp.taxfilt))

# 2. Aglomerar a nivel de Familia

gp.taxfam <- tax_glom(gp.taxfilt, taxrank = "Family")
gp.taxfam
View(sample_data(gp.taxfam))

# Transformar a abundancias relativas (%)

gp.taxfam.rel <- transform_sample_counts(gp.taxfam, function(otu) otu/sum(otu))
gp.taxfam.rel

# Subset para incluir solo muestras de: Soil, Feces, Skin

gp.filt.no.rel <- subset_samples(gp.taxfam, sample_data(gp)$SampleType == "Soil" | sample_data(gp)$SampleType == "Feces" | sample_data(gp)$SampleType == "Skin")
gp.filt <- subset_samples(gp.taxfam.rel, sample_data(gp)$SampleType == "Soil" | sample_data(gp)$SampleType == "Feces" | sample_data(gp)$SampleType == "Skin")
gp.filt
View(otu_table(gp.filt))
View(sample_data(gp.filt))
taxa_names(gp.filt)

# Diversidad alfa

head(otu_table(gp.taxfam))
head(otu_table(gp.taxfam.rel))

alfa.rich <- estimate_richness(gp.filt_prune, measures = c("Observed", "Shannon", "Simpson"))
alfa.rich$SampleType <- sample_data(gp.filt_prune)$SampleType
alfa.rich

gp.filt_prune <- prune_taxa(taxa_sums(gp.filt.no.rel) > 1, gp.filt.no.rel)

plot_richness(gp.filt_prune,x = "SampleType",measures = c("Observed","Shannon","Simpson"))

plot_richness(gp.filt_prune, x = "SampleType", color = "SampleType", measures = c("Observed","Shannon","Simpson")) + 
  geom_boxplot() + 
  theme_bw() + 
  ggtitle("Indices de diversidad alfa entre tipo de muestra") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Obs.kruskal <- kruskal.test(Observed ~ SampleType, alfa.rich)
Obs.kruskal
Shan.kruskal <- kruskal.test(Shannon ~ SampleType, alfa.rich)
Shan.kruskal
Simp.kruskal <- kruskal.test(Simpson ~ SampleType, alfa.rich)
Simp.kruskal

# Curvas de Rango-Abundancia

tax_table(gp.filt)[,5]
familia <- as.character(tax_table(gp.filt)[, 5])
duplicados <- familia[duplicated(familia)]
print(duplicados)
familia.unicos <- make.unique(familia)
taxa_names(gp.filt) <- familia.unicos
otu_table(gp.filt)

par(mar = c(10, 4, 4, 2) + 0.1)
barplot(sort(taxa_sums(gp.filt), TRUE)[1:40]/nsamples(gp.filt), las=2, cex.names = 0.63)

barplot(log10(sort(taxa_sums(gp.filt), TRUE)[1:40] + 1), las = 2, cex.names = 0.63, 
        names.arg = taxa_names(gp.filt)[1:40],  
        main = "Curvas de Rango-Abundancia (40 Taxa más abundantes)",
        ylab = "Log10(Abundancia Relativa)")

gp.soil <- subset_samples(gp.filt, sample_data(gp.filt)$SampleType == "Soil")
gp.feces <- subset_samples(gp.filt, sample_data(gp.filt)$SampleType == "Feces")
gp.skin <- subset_samples(gp.filt, sample_data(gp.filt)$SampleType == "Skin")

par(mfrow = c(1, 3), mar = c(10, 4, 4, 2) + 0.1)

barplot(log10(sort(taxa_sums(gp.soil), TRUE)[1:40] + 1), las = 2, cex.names = 0.63, 
        names.arg = taxa_names(gp.soil)[1:40],  
        main = "Curvas de Rango-Abundancia de Soil (40 Taxa más abundantes)",
        ylab = "Log10(Abundancia Relativa)")
barplot(log10(sort(taxa_sums(gp.feces), TRUE)[1:40] + 1), las = 2, cex.names = 0.63, 
        names.arg = taxa_names(gp.feces)[1:40],  
        main = "Curvas de Rango-Abundancia de feces (40 Taxa más abundantes)",
        ylab = "Log10(Abundancia Relativa)")
barplot(log10(sort(taxa_sums(gp.skin), TRUE)[1:40] + 1), las = 2, cex.names = 0.63, 
        names.arg = taxa_names(gp.skin)[1:40],  
        main = "Curvas de Rango-Abundancia de skin (40 Taxa más abundantes)",
        ylab = "Log10(Abundancia Relativa)")


# Perfil taxonómico



# Diversidad Beta