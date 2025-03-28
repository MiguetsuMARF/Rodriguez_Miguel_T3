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

dim(otu_table(gp))
otu_table(gp)[1,1]>1
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
# Preprocesamiento

gp.less5 <- prune_taxa(taxa_sums(gp) < 5, gp)
gp.less5
gp.filt <- subset_samples(gp, sample_data(gp)$SampleType == "Soil" | sample_data(gp)$SampleType == "Feces" | sample_data(gp)$SampleType == "Skin")
gp.filt <- transform_sample_counts(gp.filt, function(otu) otu/sum(otu))
View(otu_table(gp.filt))
View(sample_data(gp.filt))
# Diversidad alfa

# Curvas de Rango-Abundancia

# Perfil taxonómico

# Diversidad Beta