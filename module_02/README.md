# Module 02: Proteomics Data Analysis

This module covers the analysis of proteomics data, from data cleaning and statistical analysis to visualization and functional enrichment analysis. These are the first and essential analyses or representations one would like to do when dealing with such datasets. One would use the PCA to see if conditions are separating well, if replicates behave as they should, or if there is a potential batch effect among other things. Volcano plots depend on the statistical test used, and they are a good way to visualise proteins that may be driving the differences between conditions. Finally, by doing heatmaps you will be able to focus on specific comparisons that may be of interest. 

## Session Content

In this session, we will work with a proteomics dataset to identify differentially expressed proteins between different conditions. The analysis pipeline includes:

1.  **Data Cleaning:** Loading the raw data, renaming columns, and handling missing values.
2.  **Statistical Analysis:** Performing t-tests to identify statistically significant changes in protein abundance.
3.  **Data Visualization:** Creating various plots to visualize the results, including:
    *   **Boxplots:** To compare protein intensity across different samples.
    *   **Volcano plots:** To visualize differentially expressed proteins.
    *   **PCA plots:** To explore the relationships between samples.
    *   **Heatmaps:** To visualize the expression patterns of protein sets.
4.  **Enrichment Analysis:** Identifying enriched biological processes and pathways within the sets of differentially expressed proteins.

## Core Concepts

### Statistical Analysis

-   **T-test**: A t-test is used to determine if there is a significant difference between the means of two groups. In our case, we use it to compare the protein intensity between a control group and a treatment group. The output gives us a p-value, which indicates the probability of observing the data if there were no real difference. A small p-value (typically < 0.05) suggests that the difference is statistically significant.

-   **False Discovery Rate (FDR) Correction**: When performing many t-tests at once (one for each protein), the chance of getting a false positive (a significant p-value just by chance) increases. The FDR correction (using methods like Benjamini-Hochberg) adjusts the p-values to control for this, reducing the number of false positives. We typically use the adjusted p-value (or q-value) for significance testing.

Usually these tests are done by using the `lm` function in R or functions that depend on it, which corresponds to the linear model function. 

### Data Visualization Techniques

-   **Volcano Plot**: This plot is a powerful way to visualize the results of differential expression analysis. It plots statistical significance (-log10 of the p-value) on the y-axis versus the magnitude of change (log2 fold change) on the x-axis. This allows us to easily identify proteins that are both statistically significant (high on the y-axis) and have a large magnitude of change (far from zero on the x-axis).

-   **Principal Component Analysis (PCA)**: PCA is a technique used to reduce the dimensionality of large datasets. It transforms the data into a new set of variables (principal components) that capture the most variance in the data. By plotting the first two principal components, we can get a 2D overview of our data and see how samples cluster together based on their overall protein expression profiles. This is useful for quality control and identifying patterns in the data.

-   **Heatmap**: A heatmap is a graphical representation of data where the individual values contained in a matrix are represented as colors. In proteomics, we often use heatmaps to visualize the expression levels of a set of proteins across different samples. This can reveal patterns of up- and down-regulation and help to identify clusters of proteins that behave similarly.

### Functional Enrichment Analysis

After identifying a list of significantly up- or down-regulated proteins, we want to understand what biological roles these proteins have. Functional enrichment analysis is a method used to determine if our list of proteins is significantly enriched with proteins from specific biological pathways, molecular functions, or cellular components (e.g., "glycolysis", "DNA repair"). This helps us to interpret the biological meaning of the changes we observe. The `string_api_MULTI.py` script performs this analysis by querying the STRING database. For a more in-depth information about what this script does, please visit my [GitHub repo](https://github.com/dmartimarti/STRINGDB_analyser). 

## Repository Contents

### Code

-   `stable_script.R`: The main R script for this module. It performs data loading, cleaning, statistical analysis, and generates all the plots.
-   `live_session.R`: An R script used for live coding during the session.
-   `string_api_MULTI.py`: A Python script that interacts with the STRING database API to perform functional enrichment analysis. It takes lists of up- and down-regulated genes/proteins and returns enrichment results, network diagrams, and radar charts.

### Data

-   `data_coding_club.csv`: The primary dataset for this module, containing protein intensity data.
-   `DM_UP.txt` and `down_genes_DM.txt`: Text files containing lists of up- and down-regulated genes.
-   `DM_enrich.xlsx`: An Excel file with the enrichment results.
-   `enrich_out_DM/`: A directory containing the output from the enrichment analysis.

### Exploration

This directory contains exploratory plots and results, such as heatmaps and enrichment analysis outputs.
