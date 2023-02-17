
# libraries ---------------------------------------------------------------

library(tidyverse)
library(here)
library(cowplot)


# load data ---------------------------------------------------------------

auc = read_csv(here('data', 'Output_595', 'Summary.csv'))


head(auc)

view(auc)


# data cleaning -----------------------------------------------------------

auc = auc %>% 
  drop_na(Strain) %>% 
  select(Well, Replicate, Metformin_mM, PG, Strain, `595nm_f_AUC`) %>% 
  rename(AUC = `595nm_f_AUC`)


# some basic data transformation
auc %>% 
  mutate(AUC_log = log2(AUC))

auc['AUC_log'] = log2(auc$AUC)

# some basic data exploration
## counting how many strains we have
auc %>% 
  distinct(Strain) %>% 
  count()

length(unique(auc$Strain))


## counting how many plates are in the dataset
auc %>% 
  distinct(PG) %>% 
  count()

length(unique(auc$PG))

## how many replicates
auc %>% 
  distinct(Replicate)

length(unique(auc$Replicate))




# filtering data ---------------------------------------

## simple filtering
auc %>% 
  filter(Strain == 'NT12060') %>% 
  arrange(Metformin_mM, Replicate)

nt12060 = auc[auc$Strain == 'NT12060',]
nt12060[order(nt12060$Metformin_mM),]



## multiple filtering

strains = unique(auc$Strain)

strains[1:3]

auc %>% 
  filter(Strain %in% strains[1:3])

auc[auc$Strain %in% strains[1:3],]

















