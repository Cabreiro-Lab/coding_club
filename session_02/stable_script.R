library(tidyverse)
library(readxl)
library(here)
library(glue)
library(rstatix)
library(cowplot)
library(ComplexHeatmap)

# set up cowplot theme for the project

theme_set(theme_cowplot(14))


# load data ---------------------------------------------------------------

prot = read_csv('data/data_coding_club.csv')

prot = prot %>% rename(gene = Gene_names, 
                       sample = Sample,
                       protein_id = Protein_IDs,
                       protein_name = Protein_names,
                       intensity = Intensity) %>% 
  select(-Majority_protein_IDs)



prot = prot %>% 
  group_by(gene, protein_id, sample, protein_name) %>% 
  mutate(replicate = 1:n(), .before = protein_id) %>% 
  ungroup %>% 
  drop_na(intensity)

# how many groups
unique(prot$sample)

# how many reps
unique(prot$replicate)

prot


# clean the data ---------------------------------------

prot %>%
  group_by(gene, protein_id, sample) %>% 
  count() %>% 
  arrange(n)

# 1. first round of data cleaning: genes and groups with only 1 instance
single_instance = prot %>%
  group_by(gene, protein_id, sample) %>% 
  count() %>% 
  filter(n == 1) %>% 
  select(-n)

prot_clean = prot %>% anti_join(single_instance)

# 2. second round: genes present in only one group
single_groups = prot_clean %>%
  distinct(gene, protein_id, sample) %>% 
  group_by(gene, protein_id) %>% 
  count() %>% 
  filter(n == 1) %>% 
  select(-n)

prot_clean = prot_clean %>% 
  anti_join(single_groups)







# stats -------------------------------------------------------------------

stats = prot_clean %>% 
  group_by(gene, protein_id) %>% 
  t_test(intensity ~ sample, ref.group = 'Control',  detailed = TRUE)

stats %>% 
  adjust_pvalue(method = 'fdr') %>% 
  add_significance('p.adj') %>% view


# boxplots ---------------------------------------


prot_clean %>% 
  filter(gene == 'lon') %>% 
  ggplot(aes(x = sample, y = intensity, fill = sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  facet_wrap(~gene)

box_prot = function(genes){
  prot_clean %>% 
    filter(gene %in% genes) %>% 
    ggplot(aes(x = sample, y = intensity, fill = sample)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    facet_wrap(~gene, scales = 'free_y')
}

box_prot(c('lon', 'sucA'))





# volcano plots -----------------------------------------

# volcano plots

# a simple volcano plot
stats %>% 
  filter(group1 == 'clpX' & group2 == 'Control') %>% 
  mutate(logpval = -log10(p.adj)) %>% 
  ggplot(aes(x = estimate, y = logpval)) +
  geom_point() +
  labs(x = 'log2 Fold Change',
       y = '-log10(p-value)')


# modify within the plot
library(ggrepel)
log2F_thres = 0.5
logpval_thres = -log10(0.05)
stats %>% 
  filter(group1 == 'clpX' & group2 == 'Control') %>% 
  mutate(logpval = -log10(p.adj)) %>% 
  ggplot(aes(x = estimate, y = logpval)) +
  geom_point(aes(color = ifelse(abs(estimate) > log2F_thres & 
                                  logpval > logpval_thres, 'red', 'grey'))) +
  labs(x = 'log2 Fold Change',
       y = '-log10(p-value)') +
  scale_color_manual(values = c('grey', 'red'),
                     name = 'Significant \n proteins') +
  geom_text_repel(aes(label = ifelse(abs(estimate) > log2F_thres & 
                                       logpval > logpval_thres, gene, NA))) +
  labs(
    x = 'log2(FC)',
    y = '-log10(FDR)',
    title = 'clpX vs Control'
  )



# modify the data frame

log2F_thres = 0.5
logpval_thres = -log10(0.05)

stats %>% 
  filter(group1 == 'clpX' & group2 == 'Control') %>% 
  mutate(logpval = -log10(p.adj)) %>% 
  mutate(col_point = ifelse(logpval > logpval_thres & abs(estimate) > log2F_thres, 
                            'black', 'grey70')) %>%
  mutate(gene_labs = ifelse(logpval > logpval_thres & abs(estimate) > log2F_thres, 
                            gene, '')) %>% 
  ggplot(aes(x = estimate, y = logpval)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed',
             color = 'grey50', alpha = 0.8) +
  geom_point(aes(color = col_point), show.legend = F) +
  scale_color_manual(values = c('black', 'grey70'),
                     breaks = c('black', 'grey70')) + 
  geom_text_repel(aes(label = gene_labs)) +
  labs(
    x = 'log2(FC)',
    y = '-log10(FDR)',
    title = 'clpX vs Control'
  ) +
  theme(
    plot.title = element_text(hjust=0.5)
  )





# heatmap -----------------------------------------------------------------

library(ComplexHeatmap)

prot %>% separate_rows(kegg, sep = ';') %>% distinct(kegg) %>% view

tca_sum = prot %>% 
  filter(str_detect(kegg,'Citrate cycle')) %>% 
  group_by(gene, sample) %>% 
  summarise(mean_int = mean(intensity, na.rm = T))

tca_matrix = tca_sum %>% 
  pivot_wider(names_from = sample, values_from = mean_int) %>% 
  arrange(gene) %>% 
  column_to_rownames('gene') %>% as.matrix()



Heatmap(tca_matrix)


## calculate z-scores

Heatmap(t(scale(t(tca_matrix))))

tca_matrix = prot %>% 
  filter(str_detect(kegg,'Citrate cycle')) %>% 
  group_by(gene, sample) %>% 
  summarise(mean_int = mean(intensity, na.rm = T)) %>% 
  group_by(gene) %>% 
  mutate(scale_int = scale(mean_int)[,1]) %>% 
  select(-mean_int) %>% 
  pivot_wider(names_from = sample, values_from = scale_int) %>% 
  arrange(gene) %>% 
  column_to_rownames('gene') %>% as.matrix()

Heatmap(tca_matrix,
        name = 'Z-score',
        column_title = 'Samples',
        row_title = 'Proteins') 



Heatmap(tca_matrix,
        name = 'Z-score',
        column_title = 'Samples',
        row_title = 'Proteins',
        row_km = 3, 
        column_km = 2) 



## example with two categories


two_cat = prot %>% 
  filter(str_detect(kegg,'Citrate cycle|Purine metabolism')) %>% 
  separate_rows(kegg, sep = ';') %>% 
  distinct(gene, sample, replicate, protein_id, .keep_all = T) %>% 
  filter(str_detect(kegg,'Citrate cycle|Purine metabolism')) %>% 
  group_by(gene, sample, kegg) %>% 
  summarise(mean_int = mean(intensity, na.rm = T)) %>% 
  ungroup

n_tca = two_cat %>% distinct(gene, .keep_all = T) %>% count(kegg)

