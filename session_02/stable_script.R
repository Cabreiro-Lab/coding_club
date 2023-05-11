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

stats %>% 
  filter(group1 == 'clpX', group2 == 'Control')
  



