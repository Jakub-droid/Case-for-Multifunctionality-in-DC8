# Looking at Cambodia specifically for both DC8 and DC13

install.packages("tidyverse")
library("tidyverse")
install.packages(RColorBrewer)
library(RColorBrewer)
getwd()
setwd("~/Documents/Dissertation/Data/Unprocessed Data")
library(tidyverse)

allgenomes <- read.csv("all_genomes_metadata.csv")
DC8data <- read.csv("normalised.fulldataset.DC8_metadata_genes.csv")
DC8identical <- read.csv("DC8_identical_sequences.csv")

# Merge metadata with identical sequence info
DC8identical <- DC8identical %>%
  rename(Sample_ID = sequence_id) %>%
  rename(identical_sequences = count)

# Trimming datasets down
merged_DC_identical_full <- full_join(DC8data, DC8identical, by = "Sample_ID") %>%
  select(!ends_with(".y")) %>%
  rename_with(~ gsub(".x$", "", .), ends_with(".x")) %>%
  mutate(
    identical_sequences = ifelse(is.na(identical_sequences), 0, identical_sequences),
    label = ifelse(identical_sequences > 0, "Identical", "Unique")
  ) %>%
  select(-c(Unnamed..0, mds_3d, mds_3d_x, mds_3d_y, mds_3d_z, filename, Subtype, domain_structure, domain_count))

# Renaming with upodated subtypes
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


# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Step 1: Filter for Cambodia
cambodia_df <- merged_DC_identical_full %>%
  filter(Country == "Cambodia", !is.na(Region), !is.na(Year))


# Looking at the proportions of all sampled Pf3k genomes that are DC8 in each region
# Step 1: Total genomes per region (all subtypes)
regionall <- allgenomescambodia %>%
  group_by(Region) %>%
  summarise(TotalSamples = n()) %>%
  filter(TotalSamples > 2)

# Step 2: DC8 counts per region
dc8_region_summary <- cambodia_df %>%
  group_by(Region) %>%
  summarise(
    DC8count = n(),  # Total DC8 genomes
    IdenticalDC8 = sum(label == "Identical")  # Identical DC8s
  )

# Step 1: Regional proportions
region_proportions <- dc8_region_summary %>%
  inner_join(regionall, by = "Region") %>%
  mutate(
    DC8_proportion = (DC8count / TotalSamples) * 100,
    Identical_proportion_within_DC8 = (IdenticalDC8 / DC8count) * 100,
    Identical_proportion_overall = (IdenticalDC8 / TotalSamples) * 100
  )

# Step 2: Add a total row for all regions combined
total_row <- region_proportions %>%
  summarise(
    Region = "Total",
    DC8count = sum(DC8count),
    IdenticalDC8 = sum(IdenticalDC8),
    TotalSamples = sum(TotalSamples)
  ) %>%
  mutate(
    DC8_proportion = (DC8count / TotalSamples) * 100,
    Identical_proportion_within_DC8 = (IdenticalDC8 / DC8count) * 100,
    Identical_proportion_overall = (IdenticalDC8 / TotalSamples) * 100
  )

# Step 3: Combine total with regional data
region_proportions_combined <- bind_rows(region_proportions, total_row)

# Step 2: Pivot the data into long format
region_proportions_long <- region_proportions %>%
  pivot_longer(cols = c("DC8_proportion", "Identical_proportion_within_DC8", "Identical_proportion_overall"),
               names_to = "Proportion_Type",
               values_to = "Proportion_Value")

# Step 3: Create the plot with faceting by Region
ggplot(region_proportions_long, aes(x = Proportion_Type, y = Proportion_Value, fill = Proportion_Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Dodge for side-by-side bars
  scale_fill_manual(values = c(
    "DC8_proportion" = "#F1A340",  # Color for DC8 Proportion
    "Identical_proportion_within_DC8" = "#4575B4",  # Color for Identical Proportion within DC8
    "Identical_proportion_overall" = "#313695"  # Color for Identical Proportion Overall
  )) +
  labs(
    title = "Proportions of Identical Sequences vs DC8 in Cambodia by Region",
    x = "Proportion Type",
    y = "Proportion (%)",
    fill = "Proportion Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    legend.position = "none",  # Remove legend as we have labels
    panel.grid = element_blank()
  ) +
  facet_wrap(~ Region, scales = "fixed")  # Facet by Region


 # Summarising the number of each subtype that is identical or not
subtypes_summary <- cambodia_df %>%
  group_by(Region, Subtype, label) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Region, Subtype) %>%
  mutate(prop = count / sum(count))
library(scales)
ggplot(subtypes_summary, aes(x = Subtype, y = prop, fill = label)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  facet_wrap(~Region, ncol = 3) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = c("Identical" = "#440154", "Unique" = "#21908C")) +
  labs(
    title = "Subtype Composition per Region",
    x = "Subtype",
    y = "Proportion",
    fill = "Sequence Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

dc8_temporal_summary <- cambodia_df %>%
  filter(!is.na(Year), label %in% c("Identical", "Unique")) %>%
  mutate(Year = as.numeric(as.character(Year))) %>%  # Convert Year to numeric if it's a factor/character
  group_by(Year, label) %>%
  summarise(num_sequences = n(), .groups = "drop") %>%
  filter(num_sequences > 1)

comparison_data <- cambodia_df %>%
  filter(label %in% c("Identical", "Unique"), !is.na(Year)) %>%
  mutate(Year_group = ifelse(Year == 2011, "2011", "Other")) %>%
  group_by(Year_group, label) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = n, values_fill = 0)

yearly_summary <- cambodia_df %>%
  filter(label %in% c("Identical", "Unique"), !is.na(Year)) %>%
  group_by(Year, label) %>%
  summarise(num_sequences = n(), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = num_sequences, values_fill = list(num_sequences = 0))

# Fisher's exact to determine signficance 
fisher_test <- fisher.test(comparison_data[, c("Identical", "Unique")])

ggplot(dc8_temporal_summary, aes(x = Year, y = num_sequences, color = label, group = label)) +
  geom_line(size = 1.8) +
  geom_point(size = 3.5, shape = 21, stroke = 1.2, fill = "white") +
  scale_color_manual(values = c("Identical" = "#440154", "Unique" = "#21908C")) +
  labs(
    title = "",
    x = "Year",
    y = "Number of Sequences",
    color = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    axis.line = element_line(color = "grey40"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )




##### MAP PLOT
# Install and load if not already
install.packages(c("rnaturalearth", "rnaturalearthdata"))
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggplot2)
library(dplyr)

# Step 1: Get world admin-1 regions
adm1 <- ne_states(country = "Cambodia", returnclass = "sf")

# Step 2: Summarise identical sequence data
cambodia_summary <- merged_DC_identical_full %>%
  filter(Country == "Cambodia", label == "Identical", !is.na(Region)) %>%
  group_by(Region) %>%
  summarise(
    total_identical_sequences = n(),
    .groups = "drop"
  ) %>%
  filter(total_identical_sequences > 1)%>%
  mutate(Region = case_when(
    Region == "Ratanakiri" ~ "Rôtânôkiri",
    Region == "Preah Vihear" ~ "Preah Vihéar",
    Region == "Pursat" ~ "Pouthisat",
    Region == "Battambang" ~ "Batdâmbâng",
    Region == "Pailin" ~ "Krong Pailin",
    Region == "Western Region" ~ "Western Region"
  ))

# Step 3: Merge map with data
map_data <- adm1 %>%
  left_join(cambodia_summary, by = c("name" = "Region"))

# Step 4: Plot
ggplot(map_data) +
  geom_sf(aes(fill = total_identical_sequences), color = "black", size = 0.3) +
  scale_fill_viridis_c(option = "viridis", direction = -1, na.value = "lightgrey") +
  labs(
    title = "",
    fill = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )


###### Identical sequences with the most copies

top_identical_sequences <- cambodia_df %>%
  filter(label == "Identical") %>%  # Focus on Identical sequences
  group_by(Region,sequence,Subtype) %>%  # Group by sequence, region, and subtype (you can adjust based on your data)
  summarise(
    total_copies = n(), na.rm = TRUE)  # Summing the identical sequences
    .groups = "drop"
  ) %>%
  arrange(desc(total_copies))  # Ar

regions <- Cambodia_data %>%
  group_by(Region) %>%
  summarise(n())
