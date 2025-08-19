# Microbiome16S
Recommend making a environment just for Dada2.

```r
# Load required packages
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(decontam)
library(RColorBrewer)
library(vegan)

BiocManager::install("decontam") 
library(decontam)  # Decontamination package

# Set path to your data
Make sure - data is trimmed 
path <- "/Users/yourpath"  ## CHANGE to your fastq directory
list.files(path)

# Identify forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_trimmed.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)  # Extracts "CP_S104" from "CP_S104_R1.fastq"

# Define control sample name **without** R1/R2
control_sample_name <- "CP_S104"

# Quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,275),
                     maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy
*Implement with your own taxa assignments
taxa <- assignTaxonomy(seqtab.nochim, "/Users/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1")
rownames(taxa) <- getSequences(seqtab.nochim)
head(taxa)

## Define control sample
control_sample <- "CP"

## Create sample data frame
samdf <- data.frame(Sample = rownames(seqtab.nochim), SampleID = rownames(seqtab.nochim))
rownames(samdf) <- rownames(seqtab.nochim)
samdf <- data.frame(SampleID = rownames(seqtab.nochim))
rownames(samdf) <- rownames(seqtab.nochim)

## Add missing samples with placeholder metadata
#missing_samples <- c("201O", "213O", "222O")
#missing_samples <- c("222U")
missing_samples <- c("201N", "202N", "204N", "212N", "222N")
missing_rows <- data.frame(SampleID = missing_samples)
rownames(missing_rows) <- missing_samples

## Combine original and missing samples
samdf <- rbind(samdf, missing_rows)

## Add grouping info to samdf for ordination coloring
samdf$Group <- substr(rownames(samdf), nchar(rownames(samdf)), nchar(rownames(samdf)))
samdf$SampleID <- rownames(samdf)  # ensure SampleID is aligned with row names

## Create phyloseq object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(samdf)
)

## Order samples numerically
sample_names(ps) <- as.character(sample_names(ps))
sample_names(ps) <- sort(sample_names(ps))

## Identify contaminants using prevalence method
valid_samples <- sample_names(ps)[sample_names(ps) != control_sample & sample_names(ps) %in% rownames(seqtab.nochim)]

# Extract OTU table from ps
otu_mat <- as(otu_table(ps), "matrix")

# Get the control sample counts
control_counts <- otu_mat[control_sample, , drop=FALSE]

# Subtract control counts from all other samples
otu_mat_adj <- sweep(otu_mat, 2, control_counts, FUN="-", check.margin=FALSE)

# Replace negative values with zeros
otu_mat_adj[otu_mat_adj < 0] <- 0

# Update the phyloseq object with the adjusted OTU table
otu_table(ps) <- otu_table(otu_mat_adj, taxa_are_rows=FALSE)

## Remove control sample from phyloseq object
ps <- prune_samples(sample_names(ps) != control_sample, ps)

## Aggregate data at Genus level
ps_genus <- tax_glom(ps, taxrank = "Genus")

## Transform to relative abundance
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) if(sum(x) > 0) (x / sum(x)) * 100 else x)

## Summarize genus abundance
plot_data <- psmelt(ps_genus_rel)

## Standardize naming convention for Escherichia-Shigella
plot_data$Genus[plot_data$Genus == "Escherichia-Shigella"] <- "Escherichia_Shigella"

## Get top 20 genera by overall abundance
overall_abundance <- aggregate(Abundance ~ Genus, data = plot_data, sum)
top20_genera <- head(overall_abundance[order(-overall_abundance$Abundance), "Genus"], 20)

## Classify anything outside top 20 as 'Other'
plot_data$Genus[!plot_data$Genus %in% top20_genera] <- "Other"

## Add missing samples with zero abundance
missing_samples_in_plot <- setdiff(missing_samples, unique(plot_data$Sample))
if (length(missing_samples_in_plot) > 0) {
  columns_needed <- colnames(plot_data)
  empty_rows <- as.data.frame(matrix(NA, nrow = length(missing_samples_in_plot), ncol = length(columns_needed)))
  colnames(empty_rows) <- columns_needed
  empty_rows$Sample <- missing_samples_in_plot
  empty_rows$Abundance <- 0
  empty_rows$Genus <- "NA"
  plot_data <- rbind(plot_data, empty_rows)
}

## Order genera correctly
unique_genera <- unique(plot_data$Genus)
plot_data$Genus <- factor(plot_data$Genus, levels = c("Staphylococcus", setdiff(unique_genera, c("Staphylococcus", "Other", "NA")), "Other", "NA"))

## Generate the color palette
custom_palette_path <- "/home/bioi488/16S_Sandra/Urobiome.colors.Package.R" #Can load any other palette you want
if (file.exists(custom_palette_path)) {
  source(custom_palette_path)
  palette_colors <- urobiome.colors

  ## Add Ligilactobacillus, Azospirillum, Dolosigranulum, Sphingomonas, Phyllobacterium, Lautropia if not present
  ##These were just not present in packet
  if (!"Ligilactobacillus" %in% names(palette_colors)) {
    palette_colors["Ligilactobacillus"] <- rgb(204, 0, 204, maxColorValue = 255)
  }
  if (!"Azospirillum" %in% names(palette_colors)) {
    palette_colors["Azospirillum"] <- rgb(255, 20, 147, maxColorValue = 255)
  }
  if (!"Dolosigranulum" %in% names(palette_colors)) {
    palette_colors["Dolosigranulum"] <- rgb(128, 0, 0, maxColorValue = 255)
  }
  if (!"Sphingomonas" %in% names(palette_colors)) {
    palette_colors["Sphingomonas"] <- rgb(252, 250, 167, maxColorValue = 255)
  }
  if (!"Phyllobacterium" %in% names(palette_colors)) {
    palette_colors["Phyllobacterium"] <- rgb(255, 255, 102, maxColorValue = 255)
  }
  if (!"Lautropia" %in% names(palette_colors)) {
    palette_colors["Lautropia"] <- rgb(255, 0, 0, maxColorValue = 255)
  }
} else {
  stop("urobiome.colors file not found at", custom_palette_path)
}

## Graph 1: Genus Abundance Plot
plot1 <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(palette_colors, "NA" = "#D3D3D3"), na.translate = TRUE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Genus")
```

# Alpha diversity 
```r
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)

#Calculate alpha diversity measures
alpha_df <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

#Add metadata (e.g., anatomical site)
metadata_df <- data.frame(sample_data(ps)) %>%
  mutate(SampleID = rownames(.))

#Merge with diversity data
alpha_df <- alpha_df %>%
  mutate(SampleID = rownames(.)) %>%
  left_join(metadata_df, by = "SampleID")

#Reshape to long format (so Shannon and Simpson are stacked in one column)
alpha_long <- alpha_df %>%
  pivot_longer(cols = c("Shannon", "Simpson"),
               names_to = "Index",
               values_to = "Diversity")

#Graph 3: Alpha diversity
p_alpha <- ggplot(alpha_long, aes(x = Index, y = Diversity, fill = Index)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1) +
  theme_minimal() +
  labs(title = "Alpha Diversity (Single Site)",
       x = "Diversity Index",
       y = "Diversity Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
```

# Beta diversity 
```r
### Load required libraries
library(vegan)
library(ggplot2)
library(reshape2)

### Load CSVs
nasal  <- read.csv("NasalMarch27th_Abundance_Per_Sample.csv")
oral   <- read.csv("OralMarch27th_Abundance_Per_Sample.csv")
urine  <- read.csv("UrineagainMarch27th_Abundance_Per_Sample.csv")

### Add Source labels
nasal$Source <- "Nasal"
oral$Source  <- "Oral"
urine$Source <- "Urine"

### Combine into one long dataframe
all_abund <- rbind(nasal, oral, urine)

### Convert to wide format
abund_wide <- dcast(all_abund, Sample + Source ~ Genus, value.var = "Abundance", fill = 0)

### Separate sample metadata
sample_info <- abund_wide[, c("Sample", "Source")]
rownames(sample_info) <- sample_info$Sample

### Create abundance matrix (samples x taxa)
abund_matrix <- abund_wide[, !(names(abund_wide) %in% c("Sample", "Source"))]
rownames(abund_matrix) <- sample_info$Sample

### Remove samples with all-zero abundances
keep <- rowSums(abund_matrix) > 0
abund_matrix <- abund_matrix[keep, ]
sample_info  <- sample_info[rownames(sample_info) %in% rownames(abund_matrix), ]

### Bray-Curtis distance
bray_dist <- vegdist(abund_matrix, method = "bray")

### PCoA
pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)
var_explained <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 1)
pcoa_df <- data.frame(Sample = rownames(pcoa_result$points),
                      Axis.1 = pcoa_result$points[, 1],
                      Axis.2 = pcoa_result$points[, 2],
                      stringsAsFactors = FALSE)
pcoa_df <- merge(pcoa_df, sample_info, by = "Sample")

### Plot PCoA
p <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Source)) +
  geom_point(size = 3) +
  stat_ellipse(type = "norm", level = 0.95) +
  theme_minimal(base_size = 14) +
  labs(title = "PCoA of Microbial Communities (Bray-Curtis)",
       x = paste0("PCoA Axis 1 (", var_explained[1], "% variance)"),
       y = paste0("PCoA Axis 2 (", var_explained[2], "% variance)")) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

### Graph 2
ggsave("PCoA_RelAbund_Combined.png", p, width = 8, height = 6)

### PERMANOVA (overall)
adonis_out <- adonis2(bray_dist ~ Source, data = sample_info)

### Pairwise PERMANOVAs
nasal_urine_idx <- sample_info$Source %in% c("Urine", "Nasal")
oral_urine_idx  <- sample_info$Source %in% c("Urine", "Oral")
oral_nasal_idx  <- sample_info$Source %in% c("Oral", "Nasal")

pw_urine_nasal <- adonis2(as.dist(as.matrix(bray_dist)[nasal_urine_idx, nasal_urine_idx]) ~ Source,
                          data = sample_info[nasal_urine_idx, ])
pw_urine_oral <- adonis2(as.dist(as.matrix(bray_dist)[oral_urine_idx, oral_urine_idx]) ~ Source,
                         data = sample_info[oral_urine_idx, ])
pw_oral_nasal <- adonis2(as.dist(as.matrix(bray_dist)[oral_nasal_idx, oral_nasal_idx]) ~ Source,
                         data = sample_info[oral_nasal_idx, ])

### Adjust p-values for multiple testing
raw_pvals <- c(pw_urine_nasal$`Pr(>F)`[1],
               pw_urine_oral$`Pr(>F)`[1],
               pw_oral_nasal$`Pr(>F)`[1])

fdr_adj   <- p.adjust(raw_pvals, method = "fdr")
bonf_adj  <- p.adjust(raw_pvals, method = "bonferroni")

# Combine results
stat_results <- data.frame(
  Comparison = c("Urine vs Nasal", "Urine vs Oral", "Oral vs Nasal"),
  Raw_P = raw_pvals,
  FDR_Adjusted_P = fdr_adj,
  Bonferroni_Adjusted_P = bonf_adj
)

### Save statistics
write.csv(stat_results, "PCoA_Pairwise_PERMANOVA_Stats.csv", row.names = FALSE)
```


