# Session 01 - Bacterial Growth Assay

## Data description

In this session, we will analyze a dataset consisting of an array of E. coli strains that were grown under different metformin concentrations (0, 50, 100, and 200 mM). The original data was read using the Biospa machine in the lab, and the tables were summarized using an in-house Python script.

Given that we have several hundred strains at four conditions, each with triplicates, we will perform basic analyses on the data. For example, we will check the number of strains, conditions, and replicates per condition using R commands. After this, we will learn how to clean the data by filtering and transforming variables into a usable format. Finally, we will conclude the session by learning the basic steps of creating simple boxplots using ggplot2, and how to modify certain aspects of the plots to suit our needs.

The data looks like this once cleaned.

| **Well** | **Replicate** | **Metformin_mM** | **PG** | **Strain** | **AUC** |
|:---------|--------------:|-----------------:|-------:|:-----------|--------:|
| B2       |             1 |              200 |     10 | NT12534    |   0.488 |
| B3       |             1 |              200 |     10 | NT12540    |   0.004 |
| B4       |             1 |              200 |     10 | NT12546    |   0.961 |
| B5       |             1 |              200 |     10 | NT12557    |   0.012 |
| B6       |             1 |              200 |     10 | NT12010    |   3.510 |
| B7       |             1 |              200 |     10 | NT12016    |   0.257 |

During the course of the first session, we will learn:

-   What are vectors in R and how to index them

-   What are dataframes in R and how to index them

-   To use the basic functions from `dplyr` such as `filter` or `select`

-   How to make filter multiple items in the data using logical operators (`AND` and `OR`)

-   How to generate very simple dot plots and boxplots

-   A very brief introduction on how to modify the plot (e.g., labels, colours, font size...)

In the end you should be able to generate a plot similar to this:
[Boxplot example](https://github.com/Cabreiro-Lab/coding_club/blob/main/session_01/plots/boxplot_NT12060.pdf)
![My plot](plots/boxplot_NT12060.png)
