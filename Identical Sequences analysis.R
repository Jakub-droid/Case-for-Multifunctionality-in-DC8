# Analysis of identical sequences and Region.xs
# Load the datasets
library(tidyverse)
# -----------------------
# IDENTICAL SEQUENCE ANALYSIS SCRIPT
# -----------------------

# Load required libraries
library(dplyr)
library(ggplot2)

# Load datasets
allgenomes <- read.csv("all_genomes_metadata.csv")
DC8data <- read.csv("normalised.fulldataset.DC8_metadata_genes.csv")
DC8identical <- read.csv("DC8_identical_sequences.csv")

# Merge metadata with identical sequence info
DC8identical <- DC8identical %>%
  rename(Sample_ID = sequence_id) %>%
  rename(identical_sequences = count)

merged_DC_identical_full <- full_join(DC8data, DC8identical, by = "Sample_ID") %>%
  select(!ends_with(".y")) %>%
  rename_with(~ gsub(".x$", "", .), ends_with(".x")) %>%
  mutate(
    identical_sequences = ifelse(is.na(identical_sequences), 0, identical_sequences),
    label = ifelse(identical_sequences > 0, "Identical", "Unique")
  ) %>%
  select(-c(Unnamed..0, mds_3d, mds_3d_x, mds_3d_y, mds_3d_z, filename, Subtype, domain_structure, domain_count))

merged_DC_identical_full <- merged_DC_identical_full %>%
  rename(Subtype = updated_Subtype) %>%
  mutate(Subtype = case_when(
    Subtype == "DBLa2-CIDRa1.8 head structure with DBLg4" ~ "TIM180var1-like (EPCR + Rosetting)",
    Subtype == "DBLa2-CIDRa1.1 head structure with DBLg6" ~ "IT4var20-like (EPCR)",
    Subtype == "DBLa1.2-CIDRa1.1 head structure with DBLg6" ~ "PFD0020c-like (EPCR + gC1qR)",
    Subtype == "DBLa1.2-CIDRa1.1 head structure with DBLg4" ~ "Novel Subtype B (EPCR)",
    Subtype == "DBLa2-CIDRa1.1 head structure with DBLg4" ~ "Novel Subtype A (EPCR)",
    Subtype == "DBLa2-CIDRa1.8 head structure with DBLg6" ~ "Novel Subtype C (EPCR + Rosetting)",
    TRUE ~ "Other"
  ))

# -----------------------
# Summary: Identical vs Unique Sequences by Dataset
# -----------------------
merged_DC_identical_count <- merged_DC_identical_full %>%
  count(label) %>%
  group_by(label) %>%
  mutate(proportion = n / sum(n))

# -----------------------
# Time Trends of Identical Sequences
# -----------------------
merged_DC_identical_count_time <- merged_DC_identical_full %>%
  filter(Year != 'n/a', Year >= 2007) %>%
  count(Year, label) %>%
  group_by(Year) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

ggplot(merged_DC_identical_count_time, aes(x = Year, y = proportion, color = label, group = label)) +
  geom_line(size = 1) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Time Trends of Identical Sequences by Dataset",
    x = "Year",
    y = "Proportion",
    color = "Sequence Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -----------------------
# Country Comparison: Ghana vs Cambodia
# -----------------------
merged_DC_identical_count_time_country <- merged_DC_identical_full %>%
  filter(Country %in% c("Cambodia", "Ghana"), Year != 'n/a', Year >= 2007) %>%
  count(Country, Year, label) %>%
  group_by(Country) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

merged_DC_identical_count_time_continent <- merged_DC_identical_full %>%
  filter(Continent %in% c("Africa", "Asia"), Year != 'n/a', Year >= 2007, Year <=2013) %>%
  count(Continent, Year, label) %>%
  group_by(Continent) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

ggplot(merged_DC_identical_count_time_continent, aes(x = Year, y = n, color = label, group = label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~ Continent) +
  labs(
    title = "Identical vs Unique Sequences Over Time by Country",
    x = "Year",
    y = "Number of Sequences",
    color = "Sequence Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Plot counts
ggplot(merged_DC_identical_count_time_country, aes(x = Year, y = n, color = label, group = label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~ Country) +
  labs(
    title = "Identical vs Unique Sequences Over Time by Country",
    x = "Year",
    y = "Number of Sequences",
    color = "Sequence Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# -----------------------
# Subtype Analysis Over Time
# -----------------------
merged_DC_identical_count_time_subtype <- merged_DC_identical_full %>%
  filter(Year >= 2007) %>%
  count(updated_Subtype, Year, Source, label) %>%
  group_by(updated_Subtype) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

ggplot(merged_DC_identical_count_time_subtype, aes(x = Year, y = n, color = updated_Subtype, group = updated_Subtype)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~label) +
  labs(
    title = "Subtype Distribution of Identical vs Unique Sequences Over Time",
    x = "Year",
    y = "Number of Sequences",
    color = "Subtype"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -----------------------
# Domain Architecture Breakdown (DC8 only)
# -----------------------
subtype_structures_time_DC8 <- merged_DC_identical_full %>%
  filter(label == 'Identical', Source == 'DC8', Year >= 2007,
         identical_sequences >= 10, updated_Subtype != 'n/a') %>%
  count(updated_arc, updated_Subtype, Year, Source, label, identical_sequences) %>%
  group_by(updated_arc) %>%
  mutate(proportion = n / sum(n), total_n = sum(n)) %>%
  filter(total_n >= 8) %>%
  ungroup()

ggplot(subtype_structures_time_DC8, aes(x = Year, y = n, color = updated_arc, group = updated_arc)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~updated_Subtype) +
  labs(
    title = "Identical DC8 PfEMP1 Structures Over Time by Subtype",
    x = "Year",
    y = "Number of Sequences",
    color = "Domain Architecture"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

library(dplyr)
library(ggplot2)

library(dplyr)
library(ggplot2)

# -----------------------------
# 1. Countries with the highest number of distinct identical sequences
# -----------------------------

identical_by_country <- merged_DC_identical_full %>%
  filter(label == "Identical", !is.na(Country)) %>%
  distinct(Sample_ID, Country) %>%
  count(Country, sort = TRUE) %>%
  rename(distinct_identical_sequences = n)


sum(identical_by_country$distinct_identical_sequences)

library(dplyr)
library(ggplot2)
library(sf)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Filter for Africa and Asia only
world_filtered <- world %>%
  filter(region_un %in% c("Africa", "Asia"))

# Prepare country-level sequence data from identical_by_country
# Rename 'Country' to 'name' to match world map naming
identical_by_country_clean <- identical_by_country %>%
  rename(name = Country)

# Join sequence data with filtered world map
world_data <- world_filtered %>%
  left_join(identical_by_country_clean, by = "name")

# Plot map
ggplot() +
  geom_sf(data = world_data, aes(fill = distinct_identical_sequences), color = "white", size = 0.3) +  # Regular borders
  geom_sf(data = world_data %>% filter(!is.na(distinct_identical_sequences)), 
          aes(fill = distinct_identical_sequences), color = "black", size = 1.45) +  # Emphasized borders
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 150), na.value = "lightgrey") + 
  labs(
    title = "Distribution of Distinct Identical Sequences",
    subtitle = "Africa and Asia",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = NULL,
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )



# -----------------------------
# 2. Identical sequences with the highest number of copies
# -----------------------------

top_identical_sequences <- merged_DC_identical_full %>%
  filter(label == "Identical") %>%
  group_by(sequence) %>%
  summarise(
    identical_sequences = first(identical_sequences),
    Sample_IDs = paste(unique(Sample_ID), collapse = ", "),
    Countries = paste(unique(Country), collapse = ", "),
    Subtypes = paste(unique(Subtype), collapse = ", ")
  ) %>%
  arrange(desc(identical_sequences)) %>%
  mutate(sequence_ID = paste0("Seq_", row_number()))




# Get the number of rows (appearances) for each sequence per country
sequence_by_country <- merged_DC_identical_full %>%
  filter(label == "Identical", !is.na(Country)) %>%
  group_by(sequence, Country, identical_sequences,Subtype) %>%
  summarise(num_occurrences = n(), .groups = "drop") %>%
  mutate(sequence_label = paste("Seq", dense_rank(sequence))) %>%
  filter(identical_sequences >= 10) 




ggplot(sequence_by_country, aes(x = reorder(sequence_label, -identical_sequences), y = num_occurrences, fill = Country)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +  # Black outlines around segments
  labs(
    title = "",
    x = "",
    y = "Number of Copies",
    fill = "Country"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",  # Hide legend — you’ll handle this manually
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
  )

