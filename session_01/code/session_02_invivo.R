# libraries ---------------------------------------------------------------

# how to load variables? 
library(tidyverse)
library(here)
library(cowplot)
library(viridis)
library(openxlsx)

# install.packages('viridis')


theme_set(theme_cowplot(14))


# BACTERIAL GROWTH --------------------------------------------------------

# load data ---------------------------------------------------------------

# loading the data with a specific function
auc = read_csv(here('data', 'Output_595', 'Summary.csv'))  %>% 
  mutate(Metformin_mM = factor(Metformin_mM, 
                               levels = c(0,50,100,200)))
# there are some standards for the coding, as for example not passing certain
# length per line if you can, as it can make the code very difficult to read

auc

# data cleaning -----------------------------------------------------------

auc = auc %>% 
  drop_na(Strain) %>% # this removes any empty well from our analysis
  # this selects the columns we want to keep
  select(Well, Replicate, Metformin_mM, PG, Strain, `595nm_f_AUC`) %>% 
  rename(AUC = `595nm_f_AUC`) # to rename 


auc
# A tibble: 6,888 × 6

as.data.frame(auc)


# this is a matrix
test_matrix = matrix(0, 3, 3)

col_names = c('Jen', 'Andy', 'Kristin')

colnames(test_matrix) = col_names
test_matrix


# filtering

auc %>% 
  filter(Strain == 'NT12001')

auc[1,5]

auc[auc$Strain == 'NT12001',]


# access to an entire column
auc$Strain

auc[,5]


auc[auc$Strain == 'NT12001',]


# you can do more than one action at a time
auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(Metformin_mM)

auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(desc(Metformin_mM))

auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(Metformin_mM, Replicate)

auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(Replicate, Metformin_mM)


auc %>% 
  arrange(Well, Strain, Replicate, Metformin_mM) 


# exporting tables
auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(Replicate, Metformin_mM) %>% 
  write_csv('example_table.csv')


auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(Replicate, Metformin_mM) %>% 
  write.xlsx('example_table.xlsx')


# let's filter by several things



# booleans: AND, OR, NOT
# AND
auc %>% 
  filter(Strain == 'NT12060', Metformin_mM == 200)

auc %>% 
  filter(Strain == 'NT12060' & Metformin_mM == 200)


# OR 
auc %>% 
  filter(Strain == 'NT12001' | Strain == 'NT12060')

auc %>% 
  filter((Strain == 'NT12001' | Strain == 'NT12060') & 
           Metformin_mM == 200)

auc %>% 
  filter(Strain == 'NT12001' | Strain == 'NT12060') %>% 
  filter(Metformin_mM == 200)


# NOT
auc %>% 
  filter(Strain == 'NT12001')

auc %>% 
  filter(Strain != 'NT12001')


## filtering several strains

strain = unique(auc$Strain)

strain

strain[456]

strain[1:3]

# DONT USE THIS WAY
auc %>% 
  filter(Strain == strain[1:3]) 

# THIS IS THE WAY, use %in%
auc %>% 
  filter(Strain %in% strain[1:3]) %>% 
  arrange(Metformin_mM, Replicate, Strain) 

auc %>% 
  filter(Strain %in% c('NT12001', 'NT12060'))





# summary statistics ------------------------------------------------------

auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(Metformin_mM) %>% 
  group_by(Metformin_mM)

# Groups:   Metformin_mM [4]


auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(Metformin_mM) %>% 
  group_by(Metformin_mM) %>% 
  summarise(AUC_mean = mean(AUC),
            AUC_sd = sd(AUC),
            AUC_median = median(AUC),
            AUC_se = AUC_sd / sqrt(n()),
            number_elements = n())


auc %>% 
  filter(Strain == 'NT12001') %>% 
  arrange(Metformin_mM) %>% 
  group_by(Strain, Metformin_mM, Well) %>% 
  summarise(AUC_mean = mean(AUC),
            AUC_sd = sd(AUC),
            AUC_median = median(AUC),
            AUC_se = AUC_sd / sqrt(n()),
            number_elements = n())


# auc %>% 
#   filter(Strain == 'NT12001') %>% 
#   arrange(Metformin_mM) %>% 
#   group_by(Metformin_mM) %>% 
#   mutate(max_val = max(AUC)) %>% 
#   mutate(AUC_ratio = AUC / max_val)


# group at a bigger scale
auc %>% 
  filter(Strain %in% c('NT12001', 'NT12002')) %>% 
  arrange(Metformin_mM) %>% 
  group_by(Strain, Metformin_mM) %>% 
  summarise(AUC_mean = mean(AUC),
            AUC_sd = sd(AUC),
            AUC_median = median(AUC),
            AUC_se = AUC_sd / sqrt(n()),
            number_elements = n())


auc.sum = auc %>% 
  arrange(Metformin_mM) %>% 
  group_by(Strain, Metformin_mM) %>% 
  summarise(AUC_mean = mean(AUC),
            AUC_sd = sd(AUC),
            AUC_median = median(AUC),
            AUC_se = AUC_sd / sqrt(n()),
            number_elements = n())


auc.sum %>% 
  write_csv('summary_stats.csv')



# let's plot! ================

auc.sum %>% 
  filter(Strain == 'NT12345') %>% 
  ggplot(aes(x = Metformin_mM, y = AUC_mean, color = Metformin_mM)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = AUC_mean - AUC_sd  ,
                    ymax = AUC_mean + AUC_sd),
                width = 0) +
  scale_color_viridis_d() +
  labs(
    x = 'Metformin (mM)',
    y = 'AUC mean (a.u.)',
    color = 'Metformin (mM)'
  ) 

ggsave('range_plot_NT12345.pdf', height = 5, width = 6)




# growth ------------------------------------------------------------------

time_data = read_csv(here('data', 'Output_595', 'Timeseries.csv')) %>%
  filter(Data == '595nm_f') %>%
  drop_na(Strain) %>% # 6,888 × 96
  pivot_longer(cols = matches('\\d'), names_to = 'Time_s', values_to = 'OD') %>% 
  # 599,256 × 11
  drop_na(OD) %>% 
  mutate(
    Time_s = as.numeric(Time_s),
    Time_h = Time_s/3600) %>%
  select(c(-File, -Data, -Reader, -Pattern)) %>%
  # rename(Strain = Str, Replicate = Replicate_y) %>%
  mutate_at(c('Well', 'Strain', 'Metformin_mM'), as.factor)


time.sum = time_data %>% 
  group_by(Strain, PG,  Metformin_mM, Well, Time_h) %>% 
  summarise(Mean = mean(OD),
            SD = sd(OD))


time.sum %>% 
  filter(Strain == 'NT12001') %>% 
  ggplot(aes(x = Time_h, y = Mean, color = Metformin_mM)) +
  geom_line()


# ------- FROM HERE -------

time.sum %>% 
  ungroup %>% 
  filter(Strain == 'NT12005') %>% 
  ggplot(aes(x = Time_h, y = Mean, 
             color = Metformin_mM, fill = Metformin_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), 
              color = NA , alpha = 0.4) +
  geom_line() +
  scale_color_viridis_d(name = 'Metformin\n (in mM)') +
  scale_fill_viridis_d(name = 'Metformin\n (in mM)') +
  labs(
    y = 'OD (a.u.)',
    x = 'Time (in hours)'
  )


plot_single_growth = function(strain = 'NT12001') {
  time.sum %>% 
    ungroup %>% 
    filter(Strain == strain) %>% 
    ggplot(aes(x = Time_h, y = Mean, 
               color = Metformin_mM, fill = Metformin_mM)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), 
                color = NA , alpha = 0.4) +
    geom_line() +
    scale_color_viridis_d(name = 'Metformin\n (in mM)') +
    scale_fill_viridis_d(name = 'Metformin\n (in mM)') +
    labs(
      y = 'OD (a.u.)',
      x = 'Time (in hours)',
      title = strain
    )
}

plot_single_growth(strain = 'NT12061')



# plot an entire PG

time.sum %>% 
  ungroup %>% 
  # filter(Strain == 'NT12005') %>% 
  filter(PG == 7) %>% 
  ggplot(aes(x = Time_h, y = Mean, 
             color = Metformin_mM, fill = Metformin_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), 
              color = NA , alpha = 0.4) +
  geom_line() +
  scale_color_viridis_d(name = 'Metformin\n (in mM)') +
  scale_fill_viridis_d(name = 'Metformin\n (in mM)') +
  labs(
    y = 'OD (a.u.)',
    x = 'Time (in hours)'
  ) +
  facet_wrap(vars(Well), nrow = 6, ncol = 10)




## simple stats ---------



plot_single_growth(strain = 'NT12061')

stats_test_single = auc %>% 
  filter(Strain == 'NT12061') %>% 
  filter(Metformin_mM %in% c(100, 50)) %>% 
  arrange(Metformin_mM)


stats_test_single %>% 
  ggplot(aes(x = Metformin_mM, y = AUC)) +
  geom_point() +
  ylim(0, 9)

model = lm(AUC ~ Metformin_mM, data = stats_test_single)
summary(model)


# comparison of all means

plot_single_growth(strain = 'NT12061')

stats_test = auc %>% 
  filter(Strain == 'NT12061') %>% 
  arrange(Metformin_mM)


stats_test %>% 
  ggplot(aes(x = Metformin_mM, y = AUC)) +
  geom_point() +
  ylim(0, 15)

stats_test$Metformin_mM

model = lm(AUC ~ Metformin_mM, data = stats_test)
summary(model)

# changing the reference level for comparisons
stats_test$Metformin_mM = relevel(stats_test$Metformin_mM, ref = '50')

stats_test$Metformin_mM

model = lm(AUC ~ Metformin_mM, data = stats_test)
summary(model)



# pairwise comparisons

model_anova = aov(AUC ~ Metformin_mM, data = stats_test)
summary(model_anova)

tukey_anova = TukeyHSD(model_anova)

library(broom)

tidy(tukey_anova) %>% 
  mutate(p.stars = gtools::stars.pval(adj.p.value))



# stats using tidyverse

library(rstatix)

stats_test %>% 
  group_by(Strain) %>% 
  t_test(AUC ~ Metformin_mM, detailed = TRUE) %>% 
  adjust_pvalue(method = 'fdr') %>% 
  add_significance()



strains_stats = auc %>% 
  group_by(Strain) %>% 
  t_test(AUC ~ Metformin_mM, detailed = TRUE) %>% 
  adjust_pvalue(method = 'fdr') %>% 
  add_significance()

strains_stats %>% 
  write_csv('table/strains_stats.csv')











  






