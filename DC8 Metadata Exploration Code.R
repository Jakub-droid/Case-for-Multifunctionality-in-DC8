# ---------------------------------------------
# Load Libraries
# ---------------------------------------------
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("RColorBrewer")) install.packages("RColorBrewer")

library(tidyverse)
library(RColorBrewer)

# ---------------------------------------------
# Set Working Directory & Load Data
# ---------------------------------------------
setwd("~/Documents/Dissertation/Data/Unprocessed Data")

DC8data <- read_csv("normalised.fulldataset.DC8_metadata_genes.csv") %>%
  rename(Subtype = updated_Subtype)

allgenomes <- read_csv("all_genomes_metadata.csv")

# ---------------------------------------------
# Country Sample Size Summary
# ---------------------------------------------
country_sample_size <- allgenomes %>%
  count(Country, name = "Sample_Count")

# ---------------------------------------------
# Geographic Distribution of DC8 Samples
# ---------------------------------------------
DC8countries <- DC8data %>%
  count(Country, Continent, name = "count") %>%
  mutate(total = nrow(DC8data),
         percentage = (count / total) * 100)

# ---------------------------------------------
# Helper Function: Subtype Summary by Country
# ---------------------------------------------
summarise_subtypes_by_country <- function(data, country, total_rows) {
  data %>%
    filter(Country == country) %>%
    mutate(total_country_rows = n()) %>%
    group_by(Region, Subtype) %>%
    summarise(
      TotalN = n(),
      Percentage_Of_Overall = TotalN / total_rows * 100,
      Percentage_Of_Country = TotalN / unique(total_country_rows) * 100,
      Subtype_Count = n_distinct(Subtype),
      Subtypes = toString(unique(Subtype)),
      Genome_Count = n_distinct(Genome),
      Genomes = toString(unique(Genome)),
      .groups = "drop"
    )
}

Ghana_summary <- summarise_subtypes_by_country(DC8data, "Ghana", nrow(DC8data))
Cambodia_summary <- summarise_subtypes_by_country(DC8data, "Cambodia", nrow(DC8data))

# ---------------------------------------------
# Subtype Distribution (Overall)
# ---------------------------------------------
DC8subtypes <- DC8data %>%
  group_by(Subtype) %>%
  summarise(
    number = n(),
    TotalN = nrow(DC8data),
    percentage = (number / TotalN) * 100,
    genome_count = n_distinct(Genome),
    genomes = toString(unique(Genome)),
    .groups = "drop"
  )

ggplot(DC8subtypes, aes(x = "", y = percentage, fill = Subtype)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Proportion of DC8 Subtypes") +
  scale_fill_brewer(palette = "Set3")

# ---------------------------------------------
# Subtypes per Genome - Filtered by Country Sample Size
# ---------------------------------------------
filter_genome_subtypes <- function(data, region_filter = NULL) {
  data %>%
    filter(if (!is.null(region_filter)) Continent == region_filter else TRUE) %>%
    group_by(Country, Genome) %>%
    summarise(Subtype_Count = n_distinct(Subtype), .groups = "drop") %>%
    add_count(Country, name = "country_genome_count") %>%
    filter(country_genome_count >= 20) %>%
    count(Country, Subtype_Count, name = "genomecount") %>%
    mutate(Subtype_Count = as.factor(Subtype_Count))
}

DC8_genomes_africa <- filter_genome_subtypes(DC8data, "Africa")
DC8_genomes_asia <- filter_genome_subtypes(DC8data, "Asia")

# ---------------------------------------------
# Plot: Subtypes per Genome by Region
# ---------------------------------------------
plot_subtype_per_genome <- function(data, region_name) {
  ggplot(data, aes(x = Subtype_Count, y = genomecount)) +
    geom_col(position = "dodge", aes(fill = Country)) +
    labs(
      title = paste("Number of DC8 Subtypes per Genome in", region_name),
      x = "Number of DC8 Subtypes",
      y = "Number of Genomes"
    ) +
    facet_wrap(~Country) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_subtype_per_genome(DC8_genomes_africa, "Africa")
plot_subtype_per_genome(DC8_genomes_asia, "Asia")

# ---------------------------------------------
# Subtype by Geography
# ---------------------------------------------
DC8SubtypeGeography <- DC8data %>%
  count(Subtype, Country, name = "count") %>%
  mutate(percentage = (count / nrow(DC8data)) * 100)

# ---------------------------------------------
# Age-Based Summary
# ---------------------------------------------
DC8Agedata <- DC8data %>%
  group_by(Year, Country) %>%
  summarise(
    Subtypes = toString(unique(Subtype)),
    Count = n(),
    .groups = "drop"
  ) %>%
  group_by(Year) %>%
  summarise(
    Countries = toString(unique(Country)),
    SubtypesByCountry = paste0(Country, ": ", Subtypes, collapse = "; "),
    TotalCountries = n_distinct(Country),
    TotalSubtypes = length(unique(unlist(strsplit(Subtypes, ", ")))),
    .groups = "drop"
  )

# ---------------------------------------------
# Proportion of DC8 Genomes in Asia & Africa
# ---------------------------------------------
calculate_proportion_plot <- function(region) {
  dc8_genomes <- DC8data %>%
    select(-Sample_ID, -Subtype) %>%
    distinct() %>%
    count(Continent, Country, name = "DC8genome_count") %>%
    filter(Continent == region, Country != "")
  
  total_genomes <- allgenomes %>%
    count(Continent, Country, name = "allgenomes_count") %>%
    filter(Continent == region)
  
  proportions <- dc8_genomes %>%
    inner_join(total_genomes, by = c("Continent", "Country")) %>%
    mutate(
      proportion = DC8genome_count / allgenomes_count * 100,
      label = paste(Country, allgenomes_count, sep = "_")
    ) %>%
    filter(allgenomes_count >= 20)
  
  average <- mean(proportions$proportion)
  
  ggplot(proportions, aes(x = label, y = proportion)) +
    geom_col(aes(fill = Country)) +
    geom_hline(yintercept = average, color = "red", linetype = "dashed", size = 1) +
    scale_fill_brewer(palette = "Set3") +
    theme_classic() +
    labs(title = paste("Proportion of DC8 in", region), x = "Country", y = "Proportion of DC8") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
}

calculate_proportion_plot("Asia")
calculate_proportion_plot("Africa")

# ---------------------------------------------
# Country-Level DC8 Counts
# ---------------------------------------------
ggplot(DC8countries, aes(x = reorder(Country, -count), y = count, fill = Country)) +
  geom_col() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "DC8 Instances per Country", x = "Country", y = "DC8 Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
