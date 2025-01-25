library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)
library(ggthemes)
library(openxlsx)
library(dplyr)

# import data
data <-read.csv("data/GSE21942_normalized_processed.csv")


count_data <- as.data.frame(count_data)

before_normalized_long <- as.data.frame(count_data) |> 
  rownames_to_column("Gene") |> 
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") |> 
  mutate(Status = "Before Normalization")

normalized_count_data <- log2(count_data + 1)

after_normalized_long <- as.data.frame(normalized_count_data) |> 
  rownames_to_column("Gene") |> 
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") |> 
  mutate(Status = "After Normalization")

# Combine both datasets
merged_data <- bind_rows(before_normalized_long, after_normalized_long) 

# Boxplot visualization
ggplot(merged_data, aes(x = Expression, y = Status)) +
  geom_density() 
ggsave("figure/expression_before_after_normalization.png")

  
