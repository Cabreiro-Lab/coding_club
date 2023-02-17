
# libraries ---------------------------------------------------------------

# how to load variables? 
library(tidyverse)
library(here)
library(cowplot)



# BASICS ------------------------------------------------------------------

## VECTORS ---------------------------

# Create a vector of numbers from 1 to 5
my_vector = c(1, 2, 3, 4, 5)

# Print the vector
print(my_vector)
my_vector

# Assign a new value to the third element of the vector
my_vector[3] = 10

# Print the modified vector
print(my_vector)

# Calculate the mean of the vector
mean_value = mean(my_vector)

# Print the mean
print(mean_value)

## DATAFRAME ----------------------------

# Create a vector for the first column
col1 <- c("A", "B", "C", "D")

# Create a vector for the second column
col2 <- c(1, 2, 3, 4)

# Create a vector for the third column
col3 <- c(TRUE, FALSE, FALSE, TRUE)

# Combine the vectors into a data frame
my_df <- data.frame(col1, col2, col3)

# Print the data frame
print(my_df)

## INDEXING -----------------------------

# Create a sample data frame
my_df <- data.frame(col1 = c("A", "B", "C", "D"), 
                    col2 = c(1, 2, 3, 4), 
                    col3 = c(TRUE, FALSE, FALSE, TRUE))
my_df

# Index the first row and second column using row and column numbers
my_df[1, 2]

# Index the entire second column using the column name
my_df$col2

# Index the third row using row number
my_df[3, ]

# Index the first and second columns
my_df[, c(1,2)]

# Index the rows where col3 is TRUE
my_df[my_df$col3 == TRUE, ]

# Index the rows where col2 is greater than 2 and col3 is FALSE
my_df[my_df$col2 > 2 & my_df$col3 == FALSE, ]





# BACTERIAL GROWTH --------------------------------------------------------

# load data ---------------------------------------------------------------

# loading the data with a specific function
auc = read_csv(here('data', 'Output_595', 'Summary.csv'))  %>% 
  mutate(Metformin_mM = factor(Metformin_mM, 
                               levels = c(0,50,100,200)))
# there are some standards for the coding, as for example not passing certain
# length per line if you can, as it can make the code very difficult to read


head(auc)

view(auc)


# data cleaning -----------------------------------------------------------

auc = auc %>% 
  drop_na(Strain) %>% # this removes any empty well from our analysis
  # this selects the columns we want to keep
  select(Well, Replicate, Metformin_mM, PG, Strain, `595nm_f_AUC`) %>% 
  rename(AUC = `595nm_f_AUC`) # to rename 


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
subset_auc <- auc[auc$Strain %in% strains[1:3], ]

# Order the subset by Metformin_mM and Replicate
ordered_subset_auc <- subset_auc[order(subset_auc$Metformin_mM, 
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

# a bit more complex
# auc %>% 
#   filter(Strain == 'NT12060', Replicate == 1) %>% 
#   ggplot(aes(x = Metformin_mM, y = AUC)) +
#   geom_point(color = 'red', 
#              size = 4) +
#   geom_label(aes(label = AUC))



# let's make an example of a boxplot by isolating one of the Strains

test_boxplot = auc %>% 
  filter(Strain == 'NT12060')

# the most basic boxplot
test_boxplot %>% 
  ggplot(aes(x = Metformin_mM, y = AUC)) +
  geom_boxplot()


# a bit more colorful
test_boxplot %>% 
  ggplot(aes(x = Metformin_mM, y = AUC, color = Metformin_mM)) +
  geom_boxplot()


# let's plot the points
test_boxplot %>% 
  ggplot(aes(x = Metformin_mM, y = AUC, color = Metformin_mM)) +
  geom_point() 

test_boxplot %>% 
  ggplot(aes(x = Metformin_mM, y = AUC, color = Metformin_mM)) +
  geom_point() +
  geom_boxplot()
# where are the points? !!!

# personal choices...
test_boxplot %>% 
  ggplot(aes(x = Metformin_mM, y = AUC, fill = Metformin_mM)) +
  geom_boxplot() +
  geom_point(size = 2)
# points look a bit... meh

# change the x-position of the points
test_boxplot %>% 
  ggplot(aes(x = Metformin_mM, y = AUC, fill = Metformin_mM)) +
  geom_boxplot() +
  geom_point(size = 2, position = position_jitterdodge())




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


test_boxplot %>% 
  ggplot(aes(x = Metformin_mM, y = AUC, fill = Metformin_mM)) +
  geom_boxplot(color = 'black') +
  geom_point(size = 2, position = position_jitterdodge()) +
  labs(
    x = 'Metformin (mM)',
    y = 'Area Under the Curve (AUC, arbitrary units)'
  ) +
  viridis::scale_fill_viridis(begin = 0.3, 
                              discrete = TRUE) +
  theme_half_open(14) +
  background_grid()

#### save the plot!!! -------
ggsave('plots/boxplot_NT12060.pdf', height = 5, width = 6)


