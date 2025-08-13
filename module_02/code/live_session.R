library(tidyverse)
library(readxl)
library(here)
library(glue)
library(rstatix)
library(cowplot)
library(ComplexHeatmap)
library(ggrepel)

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



boxprot('cheA')

boxprot(c('gltA', 'clpX'))



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


logpval_thres = -log10(0.01)
log2F_thres = 1

stats %>% 
  filter(group1 == 'clpX' & group2 == 'Control') %>% 
  mutate(logpval = -log10(p.adj)) %>% 
  ggplot(aes(x = estimate, y = logpval)) +
  geom_hline(yintercept = logpval_thres, linetype = 'dashed', color = 'grey80') + 
  geom_vline(xintercept = log2F_thres, linetype = 'dashed', color = 'grey80') +
  geom_vline(xintercept = -log2F_thres, linetype = 'dashed', color = 'grey80') +
  geom_point(aes(colour = ifelse(logpval > logpval_thres &
                                   abs(estimate) > log2F_thres, 
                                 'blue', 'grey'))) +
  scale_colour_manual(name = 'Significant\nProteins',
                      values = c('blue' = 'blue',
                                 'grey' = 'grey')) +
  geom_text_repel(aes(label = ifelse(logpval > logpval_thres &
                                       abs(estimate) > log2F_thres, 
                                     gene, NA)
                      )
                  ) +
  labs(
    x = 'log2(FC)',
    y = '-log10(FDR)',
    title = 'clpX vs Control'
  )



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

# str detect is very useful, comes from stringr 
prot %>% 
  filter(str_detect(kegg, 'Citrate cycle'))  %>% 
  select(-protein_names)

# separates rows by their separator and multiply rows
prot %>% 
  separate_rows(kegg, sep = ';') %>% 
  filter(kegg == 'Citrate cycle (TCA cycle)')


tca_mat = prot %>% 
  filter(str_detect(kegg, 'Citrate cycle'))  %>% 
  select(-protein_names) %>% 
  group_by(gene, sample) %>% 
  summarise(int_mean = mean(intensity, na.rm = T)) %>% 
  ungroup %>% 
  pivot_wider(names_from = sample, values_from = int_mean) %>% 
  column_to_rownames('gene') %>% as.matrix


typeof(tca_mat)

Heatmap(tca_mat)

# wrong scaling!
Heatmap(scale(tca_mat))

scale(tca_mat)

Heatmap(t(scale(t(tca_mat))))


# tidyverse solution
tca_matrix = prot %>% 
  filter(str_detect(kegg, 'Citrate cycle'))  %>% 
  select(-protein_names) %>% 
  group_by(gene, sample) %>% 
  summarise(int_mean = mean(intensity, na.rm = T)) %>% 
  ungroup %>% 
  group_by(gene) %>% 
  mutate(z_score = scale(int_mean)[,1]) %>% 
  select(-int_mean) %>% 
  pivot_wider(names_from = sample, values_from = z_score) %>% 
  column_to_rownames('gene') %>% as.matrix

Heatmap(tca_matrix,
        name = 'Z-score',
        column_title = 'Samples',
        row_title = 'Proteins',
        row_km = 3, 
        column_km = 2)

# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html

dev.copy2pdf(device = cairo_pdf,
             file = 'exploration/heatmap_session_5.pdf',
             height = 10, width = 9, useDingbats = FALSE)





# enrichment analysis -----------------------------------------------------


stats


stats %>% 
  filter(group1 == 'Control', group2 == 'DM') %>% 
  filter(p.adj < 0.05, abs(estimate) > 1)


prot_clean %>% 
  filter(gene == 'ibpB') %>% 
  ggplot(aes(x = sample, y = intensity, fill = sample)) +
  geom_boxplot() +
  geom_point()



DM_sig_prots = stats %>% 
  filter(group1 == 'Control', group2 == 'DM') %>% 
  filter(p.adj < 0.05) %>% 
  filter(abs(estimate) > 1)



DM_sig_prots %>% 
  filter(estimate > 0) %>% 
  select(gene) %>% 
  write_delim("data/DM_UP.txt")


DM_UP = DM_sig_prots %>% 
  filter(estimate > 0) %>% 
  select(gene) %>% 
  rename(genes = gene)

DM_DOWN = DM_sig_prots %>% 
  filter(estimate < 0) %>% 
  select(gene) %>% 
  rename(genes = gene)


list_of_datasets = list(
  'DM_UP' = DM_UP,
  'DM_DOWN' = DM_DOWN
)


library(openxlsx)

write.xlsx(file = 'data/DM_enrich.xlsx', list_of_datasets)



