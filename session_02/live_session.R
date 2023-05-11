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


names(prot)

prot = prot %>% 
  rename(gene = Gene_names,
         sample = Sample,
         intensity = Intensity,
         protein_id = Protein_IDs,
         protein_names = Protein_names, 
         kegg = KEGG_name) %>% 
  select(-Majority_protein_IDs)
  

# how many samples are 
unique(prot$sample)

# how many genes are in the dataset
unique(prot$gene)

length(unique(prot$gene))

prot %>% distinct(gene) %>% count()


# replicates
prot %>% 
  group_by(gene, protein_id, sample) %>% 
  count() %>% 
  arrange(desc(n))
  
# create the replicates info
prot = prot %>% 
  group_by(gene, protein_id, sample) %>% 
  mutate(replicate = 1:n(),
         .before = protein_id) %>% 
  ungroup




# clean the data ----------------------------------------------------------

# 1st step of data cleaning, removing single instances of genes/samples

simple_instances = prot %>% 
  drop_na(intensity) %>% 
  group_by(gene, protein_id, sample) %>% 
  count() %>% arrange(n) %>% 
  filter(n == 1) %>% 
  select(-n)

prot_filt = prot %>% 
  # drop again the NAs so you clean the target data
  drop_na(intensity) %>% 
  # anti_join to remove the groups filtered in the previous step
  anti_join(simple_instances) 

# 2nd step: remove single groups

sample_removals = prot_filt %>% 
  distinct(gene, protein_id, sample) %>% 
  group_by(gene, protein_id) %>% 
  count() %>% 
  filter(n == 1) %>% 
  select(-n)

prot_filt = prot_filt %>% 
  anti_join(sample_removals) 





# stats -------------------------------------------------------------------

proteins = unique(prot_filt$gene)

prot_filt %>% 
  filter(gene %in% proteins[1:10]) %>% 
  group_by(gene) %>% 
  t_test(intensity ~ sample, detailed = TRUE) %>% view

# entire dataset
# THIS CAN TAKE A BIT OF TIME
stats = prot_filt %>% 
  group_by(gene) %>% 
  t_test(intensity ~ sample, detailed = TRUE)


stats = stats %>% 
  adjust_pvalue(method = 'fdr') %>% 
  add_significance('p.adj')



# check data with boxplots

prot_filt %>% 
  filter(gene == 'clpX') %>% 
  ggplot(aes(x = sample, y = intensity, fill = sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  facet_wrap(~gene)


boxprot = function(genes) {
  prot_filt %>% 
    filter(gene %in% genes) %>% 
    ggplot(aes(x = sample, y = intensity, fill = sample)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    facet_wrap(~gene, scales = 'free_y')
}



boxprot('gltA')

boxprot(c('gltA', 'clpX'))




# volcano plots


stats %>% 
  filter(group1 == 'clpX' & group2 == 'Control') %>% 
  mutate(logpval = -log10(p.adj)) %>% 
  ggplot(aes(x = estimate, y = logpval)) +
  geom_point() +
  labs(x = 'log2 Fold Change',
       y = '-log10(p-value)')












