# Reference genome analysis

install.packages("tidyverse")
library("tidyverse")
install.packages(RColorBrewer)
library(RColorBrewer)
getwd()
setwd("~/Documents/Dissertation/Data/Unprocessed Data")

DC8data <- read.csv("normalised.fulldataset.DC8_metadata_genes.csv")

allgenomes <- read.csv("all_genomes_metadata.csv")

# First need to load in a list of the reference genomes

referencegenomes <- read.csv("Kenya_Rask2010_allgenomes.csv")

referencegenomelist <- referencegenomes %>%
  select(Genome) %>%
  group_by(Genome)
  
DC8_genes_reference <- DC8data %>%
  inner_join(referencegenomelist, by = "Genome")

# View results
head(DC8_genes_reference)

allgenomesreferencelist <- allgenomes %>%
  inner_join(referencegenomelist,by = "Genome")

# Calculate proportions
proportion_reference <- (nrow(DC8_genes_reference) / nrow(allgenomesreferencelist)) *100
proportion_DC8 <- (nrow(DC8data) / nrow(allgenomes)) *100

# Create a data frame for plotting
proportion_df <- data.frame(
  Category = c("Reference Genomes", "DC8 Data"),
  Proportion = c(proportion_reference, proportion_DC8)
)

# Plot the proportions
ggplot(proportion_df, aes(x = Category, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of Reference Genomes and DC8 Data",
       y = "Proportion",
       x = "Category") +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme(legend.position = "none")

# Summarize counts of all genomes per country
allgenomes_country_count <- allgenomes %>%
  group_by(Country) %>%
  summarise(total_genomes = n())

# Summarize counts of DC8 data per country (including Continent if needed)
DC8_country_count <- DC8data %>%
  group_by(Continent, Country) %>%
  summarise(DC8_genomes = n())

# Join the two summaries by Country and calculate the proportion
DC8proportions_country <- DC8_country_count %>%
  left_join(allgenomes_country_count, by = "Country") %>%
  mutate(Proportion = (DC8_genomes / total_genomes) *100)

# View the result
print(DC8proportions_country)


combined_df <- bind_rows(proportion_df, DC8proportions_country)

# Create country labels with sample size
country_labels <- allgenomes_country_count %>%
  mutate(Country_Label = paste0(Country, " (n=", total_genomes, ")")) %>%
  select(Country, Country_Label)

# Add "Overall" label separately
overall_label <- data.frame(
  Country = "Overall",
  Country_Label = "Overall (n=" %>% paste0(nrow(allgenomes), ")")
)

# Combine labels
country_labels <- bind_rows(country_labels, overall_label)



# Merge with combined_df to replace Country with Country_Label
# Plus add in reference genome proportion (manual)
reference_annotation_manual <- data.frame(
  Proportion = as.numeric("65"),  # Replace with actual reference genome names
  Category = "Reference Genomes (manual)"
)


combined_df <- combined_df %>%
  left_join(country_labels, by = "Country")


combined_df <- bind_rows(combined_df, reference_annotation_manual)

combined_df <- combined_df %>%
  mutate(Category = case_when(
    Continent == "Africa" ~ "African",
    Continent == "Asia" ~ "Asian",
    TRUE ~ Category  # Keep existing categories for other values
  ))

preserved_rows <- combined_df %>%
  filter(Category %in% c("Reference Genomes", "Reference Genomes (manual)", "DC8 Data"))

cleaned_combined_df <- combined_df %>%
  drop_na() %>%
  filter(total_genomes >= 20)
# Combine the cleaned data with the preserved rows
final_combined_df <- bind_rows(cleaned_combined_df, preserved_rows)



ggplot(final_combined_df, aes(x = Country_Label, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge",width = 0.8) +
  theme_minimal() +
  labs(title = "Proportion of DC8 var genes in..",
       y = "Proportion",
       x = "Country") +
  scale_fill_manual(values = c(
    "African" = "purple",       # Dark purple for Africa  
    "Asian" = "plum",           # Lilac/plum for Asia  
    "Reference Genomes" = "steelblue",
    "Reference Genomes (manual)" = "tomato",  # Blue for reference genomes  
    "DC8 Data" = "orange"     # Tomato red for total data  
  )) +
  facet_wrap(~ Continent, scales = "free_x") +  # Facet by continent
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(strip.text.x = element_text(size = 12)) +  # Adjust facet titles if needed
  scale_x_discrete(drop = TRUE)   # Drop unused levels, preventing duplicate labels

