
# libraries ---------------------------------------------------------------

# how to load variables? 
library(tidyverse)
library(here)
library(cowplot)
library(viridis)


theme_set(theme_cowplot(14))



# BACTERIAL GROWTH --------------------------------------------------------

# load data ---------------------------------------------------------------

# loading the data with a specific function
auc = read_csv(here('data', 'Output_595', 'Summary.csv'))  %>% 
  mutate(Metformin_mM = factor(Metformin_mM, 
                               levels = c(0,50,100,200)))
# there are some standards for the coding, as for example not passing certain
# length per line if you can, as it can make the code very difficult to read




# data cleaning -----------------------------------------------------------

auc = auc %>% 
  drop_na(Strain) %>% # this removes any empty well from our analysis
  # this selects the columns we want to keep
  select(Well, Replicate, Metformin_mM, PG, Strain, `595nm_f_AUC`) %>% 
  rename(AUC = `595nm_f_AUC`) # to rename 


# some basic data transformation
auc %>% 
  mutate(AUC_log = log2(AUC))




# filtering data ---------------------------------------

## simple filtering ------------------
auc %>% 
  filter(Strain == 'NT12060') %>% 
  arrange(Metformin_mM, Replicate)

# Create a subset of the 'auc' data frame where Strain equals 'NT12060'
subset_auc = auc[auc$Strain == 'NT12060', ]

# Order the subset by Metformin_mM and Replicate
ordered_subset_auc = subset_auc[order(subset_auc$Metformin_mM, subset_auc$Replicate), ]


## multiple filtering ----------------
strains = unique(auc$Strain)
strains[1:3]

# Using tidyverse (dplyr)
auc %>% 
  filter(Strain %in% strains[1:3]) %>% 
  arrange(Metformin_mM, Replicate)

# Create a subset of the 'auc' data frame where Strain is in the first three 
#elements of 'strains'
subset_auc = auc[auc$Strain %in% strains[1:3], ]

# Order the subset by Metformin_mM and Replicate
ordered_subset_auc = subset_auc[order(subset_auc$Metformin_mM, 
                                      subset_auc$Replicate), ]



## filtering with operators ------------

# filtering one strain and one replicate
auc %>% 
  filter(Strain == 'NT12060' & Replicate == 1)

# same effect
auc %>% 
  filter(Strain == 'NT12060', Replicate == 1)

# multiple filtering and operators
auc %>% 
  filter(Strain == 'NT12060' & Replicate %in% c(1,2))

# filtering  by AUC values, and then also by strains
auc %>% 
  filter(AUC > 15)

auc %>% 
  filter(Strain == 'NT12060' & AUC > 10)

auc %>% 
  filter(Strain == 'NT12060' | AUC > 10)


## Special operators to filter data
auc %>% 
  filter(str_detect(Strain, "NT120"))



# PLOTS -------------------------------------------------------------------


# A very basic plot with ggplot

simple_data = auc %>% 
  filter(Strain == 'NT12060', Replicate == 1)

simple_data

simple_data %>% 
  ggplot(aes(x = Metformin_mM, y = AUC)) +
  geom_point()


auc %>% 
  filter(Strain == 'NT12060', Replicate == 1) %>% 
  ggplot(aes(x = Metformin_mM, y = AUC)) +
  geom_point(color = 'red', 
             size = 4)



# let's make an example of a boxplot by isolating one of the Strains

test_boxplot = auc %>% 
  filter(Strain == 'NT12060')

# a bit more advanced
test_boxplot %>% 
  ggplot(aes(x = Metformin_mM, y = AUC, fill = Metformin_mM)) +
  geom_boxplot(color = 'grey20') +
  geom_point(size = 2, position = position_jitterdodge()) +
  labs(
    x = 'Metformin (mM)',
    y = 'Area Under the Curve (AUC, arbitrary units)'
  ) +
  scale_fill_manual(values = c(
    '#1911F0',
             '#14A5F7',
             '#1DE09E',
             '#24F714'
  ), 
  name = 'Metformin') +
  theme_half_open(14) +
  background_grid()




## summarise ---------------------------------------------------------------

single_strain =  auc %>% 
  filter(Strain == 'NT12001')

single_strain %>% 
  group_by(Metformin_mM) %>% 
  summarise(cosa = n())

single_strain %>% 
  group_by(Metformin_mM) %>% 
  summarise(AUC_mean = mean(AUC),
            AUC_sd = sd(AUC))

mean(c(12.5, 8.95, 11.9))

### plot mean and sd

single_strain %>% 
  group_by(Metformin_mM) %>% 
  summarise(AUC_mean = mean(AUC),
            AUC_sd = sd(AUC)) %>% 
  ggplot(aes(x = Metformin_mM, y = AUC_mean, colour = Metformin_mM)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = AUC_mean - AUC_sd, 
                    ymax = AUC_mean + AUC_sd),
                width = 0) +
  scale_color_viridis_d()



# growth ------------------------------------------------------------------

time_data = read_csv(here('data', 'Output_595', 'Timeseries.csv')) %>%
  filter(Data == '595nm_f') %>%
  drop_na(Strain)  %>% 
  pivot_longer(cols = matches('\\d'), names_to = 'Time_s', values_to = 'OD') %>% 
  drop_na(OD) %>% 
  mutate(
    Time_s = as.numeric(Time_s),
    Time_h = Time_s/3600) %>%
  select(c(-File, -Data, -Reader)) %>%
  # rename(Strain = Str, Replicate = Replicate_y) %>%
  mutate_at(c('Well', 'Strain', 'Metformin_mM'), as.factor)









  



