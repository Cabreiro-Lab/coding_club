# ==============================================================================
# R SCRIPT FOR THE ANALYSIS OF PROTEOMICS DATA
# ==============================================================================
#
# Author: Jules
#
# Description:
# This script is designed to analyze proteomics data to identify differentially
# expressed proteins. It includes the following steps:
#
# 1.  **Data Loading**: Load the proteomics data from the CSV file.
# 2.  **Data Cleaning**: Preprocess and clean the data for analysis.
# 3.  **Statistical Analysis**: Perform t-tests to find significant changes.
# 4.  **Data Visualization**: Create boxplots, volcano plots, PCA plots, and
#     heatmaps to visualize the results.
# 5.  **Enrichment Analysis**: Prepare data for functional enrichment analysis.
#
# ==============================================================================
#
# LIBRARIES
#
# ==============================================================================

# a function to install and load packages
# do not use it if you are in a high-performance cluster.
# in that case, you should ask your administrator to install them for you
ipak = function(pkg){
  new.pkg = pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(c('tidyverse', 'readxl', 'here', 'glue', 'rstatix', 'cowplot',
       'ComplexHeatmap', 'ggrepel', 'broom'))


# ==============================================================================
#
# LOAD DATA
#
# ==============================================================================

# Set a general theme for all the plots
theme_set(theme_cowplot(14))

# Load the proteomics data
prot = read_csv(here('data', 'data_coding_club.csv')) %>%
  # Rename columns for easier access
  rename(
    gene = Gene_names,
    sample = Sample,
    protein_id = Protein_IDs,
    protein_name = Protein_names,
    intensity = Intensity,
    kegg = KEGG_name
  ) %>%
  # Remove unnecessary column
  select(-Majority_protein_IDs) %>%
  # Create a replicate number for each protein in each sample
  group_by(gene, protein_id, sample, protein_name) %>%
  mutate(replicate = 1:n(), .before = protein_id) %>%
  ungroup()

# ==============================================================================
#
# DATA CLEANING
#
# ==============================================================================

# Remove proteins with only a single instance across all samples
single_instance_proteins = prot %>%
  group_by(gene, protein_id, sample) %>%
  count() %>%
  filter(n == 1) %>%
  select(-n)

prot_cleaned = prot %>%
  anti_join(single_instance_proteins, by = c("gene", "protein_id", "sample"))

# Remove proteins that are only present in one sample group
single_group_proteins = prot_cleaned %>%
  distinct(gene, protein_id, sample) %>%
  group_by(gene, protein_id) %>%
  count() %>%
  filter(n == 1) %>%
  select(-n)

prot_cleaned = prot_cleaned %>%
  anti_join(single_group_proteins, by = c("gene", "protein_id"))

# ==============================================================================
#
# STATISTICAL ANALYSIS
#
# ==============================================================================

# Perform t-tests for each protein, comparing each sample to the control
stats = prot_cleaned %>%
  group_by(gene, protein_id) %>%
  t_test(intensity ~ sample, ref.group = 'Control', detailed = TRUE) %>%
  # Adjust p-values for multiple comparisons using FDR
  adjust_pvalue(method = 'fdr') %>%
  add_significance('p.adj')

# ==============================================================================
#
# DATA VISUALIZATION
#
# ==============================================================================

# --- Boxplots ---

# Function to create a boxplot for a given set of genes
plot_protein_boxplot = function(genes) {
  prot_cleaned %>%
    filter(gene %in% genes) %>%
    ggplot(aes(x = sample, y = intensity, fill = sample)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    facet_wrap(~ gene, scales = 'free_y') +
    labs(
      title = 'Protein Intensity Comparison',
      x = 'Sample',
      y = 'Intensity'
    )
}

# Example: Create a boxplot for the 'lon' gene
plot_protein_boxplot(c('lon'))


# --- Volcano Plots ---

# Function to create a volcano plot for a given comparison
plot_volcano = function(group1_filter, group2_filter) {
  log2F_thres = 0.5
  logpval_thres = -log10(0.05)

  volcano_data = stats %>%
    filter(group1 == group1_filter & group2 == group2_filter) %>%
    mutate(logpval = -log10(p.adj)) %>%
    mutate(
      color_point = ifelse(logpval > logpval_thres & abs(estimate) > log2F_thres,
                             'significant', 'not_significant'),
      gene_label = ifelse(logpval > logpval_thres & abs(estimate) > log2F_thres,
                            gene, '')
    )

  ggplot(volcano_data, aes(x = estimate, y = logpval)) +
    geom_hline(yintercept = logpval_thres, linetype = 'dashed', color = 'grey50') +
    geom_vline(xintercept = c(-log2F_thres, log2F_thres), linetype = 'dashed', color = 'grey50') +
    geom_point(aes(color = color_point), show.legend = FALSE) +
    scale_color_manual(values = c('significant' = 'red', 'not_significant' = 'grey70')) +
    geom_text_repel(aes(label = gene_label), max.overlaps = 15) +
    labs(
      x = 'log2(Fold Change)',
      y = '-log10(Adjusted p-value)',
      title = paste(group1_filter, 'vs', group2_filter)
    ) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Example: Create a volcano plot for 'clpX' vs 'Control'
plot_volcano('clpX', 'Control')


# --- PCA Plot ---

# Prepare data for PCA
prot_pca_ready = prot_cleaned %>%
  select(gene, sample, replicate, intensity) %>%
  pivot_wider(names_from = gene, values_from = intensity, values_fn = mean) %>%
  ungroup() %>%
  # Remove columns with any NA values
  select_if(~ !any(is.na(.)))

# Perform PCA
pca_fit = prot_pca_ready %>%
  select(where(is.numeric)) %>%
  prcomp(scale = TRUE)

# Plot PCA results
pca_plot = pca_fit %>%
  augment(prot_pca_ready) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color = sample)) +
  geom_point(size = 2) +
  stat_ellipse() +
  labs(
    title = 'PCA of Proteomics Data',
    x = 'PC1',
    y = 'PC2',
    color = 'Sample'
  )

pca_plot

# --- Heatmap ---

# Prepare data for heatmap (e.g., for TCA cycle proteins)
tca_proteins = prot %>%
  filter(str_detect(kegg, 'Citrate cycle')) %>%
  group_by(gene, sample) %>%
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) %>%
  group_by(gene) %>%
  # Scale the intensity values (z-score)
  mutate(scaled_intensity = scale(mean_intensity)[, 1]) %>%
  select(-mean_intensity) %>%
  pivot_wider(names_from = sample, values_from = scaled_intensity) %>%
  column_to_rownames('gene') %>%
  as.matrix()

# Create the heatmap
Heatmap(
  tca_proteins,
  name = 'Z-score',
  column_title = 'Samples',
  row_title = 'TCA Cycle Proteins',
  row_km = 3, # cluster rows into 3 groups
  column_km = 2 # cluster columns into 2 groups
)

# ==============================================================================
#
# ENRICHMENT ANALYSIS PREPARATION
#
# ==============================================================================

# Extract significantly up- and down-regulated genes for a comparison
# Example: 'DM' vs 'Control'
sig_genes_dm = stats %>%
  filter(group1 == 'DM' & group2 == 'Control') %>%
  filter(p.adj < 0.05, abs(estimate) > 1)

# Get up-regulated genes
up_genes_dm = sig_genes_dm %>%
  filter(estimate > 0) %>%
  select(gene)

# Get down-regulated genes
down_genes_dm = sig_genes_dm %>%
  filter(estimate < 0) %>%
  select(gene)

# Write the gene lists to text files
write_delim(up_genes_dm, here('data', 'DM_UP.txt'), col_names = FALSE)
write_delim(down_genes_dm, here('data', 'down_genes_DM.txt'), col_names = FALSE)

# Create an Excel file with both lists of genes
list_of_datasets = list(
  "DM_UP" = up_genes_dm,
  "DM_DOWN" = down_genes_dm
)
write.xlsx(list_of_datasets, file = here('data', 'DM_enrich.xlsx'))

# ==============================================================================
#
# END OF SCRIPT
#
# ==============================================================================
