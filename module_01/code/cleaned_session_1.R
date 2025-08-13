# ==============================================================================
# R SCRIPT FOR THE ANALYSIS OF BACTERIAL GROWTH DATA
# ==============================================================================
#
# Author: Jules
#
# Description:
# This script is designed to analyze bacterial growth data from a high-throughput
# screening experiment. It covers the following steps:
#
# 1.  **Data Loading**: Load the summary data from the experiment.
# 2.  **Data Cleaning**: Clean and preprocess the data for analysis.
# 3.  **Data Exploration**: Perform basic exploration of the data.
# 4.  **Data Visualization**: Generate plots to visualize the results,
#     including boxplots and growth curves.
# 5.  **Statistical Analysis**: Perform statistical tests to identify
#     significant differences between conditions.
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

ipak(c('tidyverse', 'here', 'cowplot', 'viridis', 'openxlsx',
       'rstatix', 'broom', 'gtools'))

# ==============================================================================
#
# LOAD DATA
#
# ==============================================================================

# Set a general theme for all the plots
theme_set(theme_cowplot(14))

# Load the summary data containing the Area Under the Curve (AUC) values
auc = read_csv(here('data', 'Output_595', 'Summary.csv'))  %>%
  # Convert Metformin_mM to a factor with specified levels
  mutate(Metformin_mM = factor(Metformin_mM,
                               levels = c(0, 50, 100, 200)))

# Load the time-series data for growth curve analysis
time_data = read_csv(here('data', 'Output_595', 'Timeseries.csv')) %>%
  # Filter for the relevant data type (595nm_f)
  filter(Data == '595nm_f') %>%
  # Remove rows with missing Strain information
  drop_na(Strain)  %>%
  # Pivot the data from wide to long format
  pivot_longer(cols = matches('\\d'), names_to = 'Time_s', values_to = 'OD') %>%
  # Remove rows with missing OD values
  drop_na(OD) %>%
  # Convert Time_s to numeric and calculate Time_h
  mutate(
    Time_s = as.numeric(Time_s),
    Time_h = Time_s / 3600
  ) %>%
  # Select and rename relevant columns
  select(-File, -Data, -Reader, -Pattern) %>%
  # Convert character columns to factors
  mutate_at(c('Well', 'Strain', 'Metformin_mM'), as.factor)

# ==============================================================================
#
# DATA CLEANING
#
# ==============================================================================

# Clean the AUC data
auc_cleaned = auc %>%
  # Remove wells with no strain information
  drop_na(Strain) %>%
  # Select relevant columns
  select(Well, Replicate, Metformin_mM, PG, Strain, `595nm_f_AUC`) %>%
  # Rename the AUC column for easier access
  rename(AUC = `595nm_f_AUC`)

# ==============================================================================
#
# DATA EXPLORATION AND MANIPULATION
#
# ==============================================================================

# --- Basic Data Exploration ---

# Count the number of unique strains
n_strains = auc_cleaned %>%
  distinct(Strain) %>%
  count()

# Count the number of unique plates (PG)
n_plates = auc_cleaned %>%
  distinct(PG) %>%
  count()

# Count the number of replicates
n_replicates = auc_cleaned %>%
  distinct(Replicate)

# --- Data Filtering Examples ---

# Filter for a specific strain
strain_nt12060 = auc_cleaned %>%
  filter(Strain == 'NT12060') %>%
  arrange(Metformin_mM, Replicate)

# Filter for multiple strains
strains_to_filter = c('NT12001', 'NT12002', 'NT12003')
multiple_strains = auc_cleaned %>%
  filter(Strain %in% strains_to_filter) %>%
  arrange(Strain, Metformin_mM, Replicate)

# Filter using logical operators
complex_filter = auc_cleaned %>%
  filter(Strain == 'NT12060' & (Replicate == 1 | Replicate == 2))

# --- Summary Statistics ---

# Calculate summary statistics for a single strain
summary_stats_nt12001 = auc_cleaned %>%
  filter(Strain == 'NT12001') %>%
  group_by(Metformin_mM) %>%
  summarise(
    AUC_mean = mean(AUC),
    AUC_sd = sd(AUC),
    AUC_median = median(AUC),
    n = n()
  )

# Calculate summary statistics for all strains
summary_stats_all = auc_cleaned %>%
  group_by(Strain, Metformin_mM) %>%
  summarise(
    AUC_mean = mean(AUC),
    AUC_sd = sd(AUC),
    n = n()
  ) %>%
  ungroup()

# Export summary statistics to a CSV file
write_csv(summary_stats_all, here('table', 'summary_stats.csv'))

# ==============================================================================
#
# DATA VISUALIZATION
#
# ==============================================================================

# --- Boxplots ---

# Create a boxplot for a single strain
plot_boxplot_nt12060 = auc_cleaned %>%
  filter(Strain == 'NT12060') %>%
  ggplot(aes(x = Metformin_mM, y = AUC, fill = Metformin_mM)) +
  geom_boxplot(color = 'black') +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(
    title = 'Growth of Strain NT12060 with Metformin',
    x = 'Metformin (mM)',
    y = 'Area Under the Curve (AUC)'
  ) +
  scale_fill_viridis_d(name = 'Metformin (mM)') +
  theme(legend.position = 'top')

# Save the plot
ggsave(here('plots', 'boxplot_NT12060.png'), plot_boxplot_nt12060,
       width = 8, height = 6)
ggsave(here('plots', 'boxplot_NT12060.pdf'), plot_boxplot_nt12060,
       width = 8, height = 6)


# --- Growth Curves ---

# Summarize time-series data for plotting
time_summary = time_data %>%
  group_by(Strain, PG, Metformin_mM, Well, Time_h) %>%
  summarise(Mean_OD = mean(OD),
            SD_OD = sd(OD)) %>%
  ungroup()

# Function to plot growth curve for a single strain
plot_single_growth = function(strain_id, data = time_summary) {
  strain_data = data %>%
    filter(Strain == strain_id)

  ggplot(strain_data, aes(x = Time_h, y = Mean_OD,
                          color = Metformin_mM, fill = Metformin_mM)) +
    geom_ribbon(aes(ymin = Mean_OD - SD_OD, ymax = Mean_OD + SD_OD),
                color = NA, alpha = 0.3) +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 24, by = 6)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    labs(
      title = paste('Growth Curve for Strain', strain_id),
      x = 'Time (hours)',
      y = 'Optical Density (OD)'
    ) +
    theme(legend.position = 'top')
}

# Example: Plot growth curve for strain NT12001
plot_growth_nt12001 = plot_single_growth('NT12001')
plot_growth_nt12001

# Plot growth curves for an entire plate
plot_plate_growth = function(plate_id, data = time_summary) {
  plate_data = data %>%
    filter(PG == plate_id)

  ggplot(plate_data, aes(x = Time_h, y = Mean_OD, color = Metformin_mM)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ Well, ncol = 12, nrow = 8) +
    labs(
      title = paste('Growth Curves for Plate', plate_id),
      x = 'Time (hours)',
      y = 'Optical Density (OD)'
    ) +
    theme(legend.position = 'top')
}

# Example: Plot growth curves for plate 1
plot_growth_pg1 = plot_plate_growth(1)
plot_growth_pg1


# ==============================================================================
#
# STATISTICAL ANALYSIS
#
# ==============================================================================

# --- T-tests for all strains ---

# Perform pairwise t-tests for each strain against the control (0 mM Metformin)
stats_all_strains = auc_cleaned %>%
  group_by(Strain) %>%
  t_test(AUC ~ Metformin_mM, ref.group = '0') %>%
  adjust_pvalue(method = 'fdr') %>%
  add_significance('p.adj')

# Export the statistical results to a CSV file
write_csv(stats_all_strains, here('table', 'strains_stats.csv'))


# --- ANOVA and Tukey's HSD for a single strain ---

# Example with strain NT12061
stats_nt12061 = auc_cleaned %>%
  filter(Strain == 'NT12061')

# Perform ANOVA
anova_model = aov(AUC ~ Metformin_mM, data = stats_nt12061)
summary(anova_model)

# Perform Tukey's HSD test for pairwise comparisons
tukey_results = TukeyHSD(anova_model)
print(tukey_results)

# Tidy the Tukey HSD results into a data frame
tukey_df = tidy(tukey_results)
print(tukey_df)

# ==============================================================================
#
# END OF SCRIPT
#
# ==============================================================================
