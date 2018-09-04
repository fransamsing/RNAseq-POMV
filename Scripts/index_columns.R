# create a column for each index
install.packages("tidyverse")
library(tidyverse)
raw_data <- read_csv('index_columns.csv')

separate(raw_data, col, sep = ",", remove = FALSE)
