library(tidyverse)

setwd("~/Documents/Dissertation/Data/Unprocessed Data/Domain Identical Sequences")

# Initial Code to combine files together
# Reading files
DBLa2 <- read.csv("DBLa2_identical_sequences.csv")
DBLa1.2 <- read.csv("DBLa1.2_identical_sequences.csv") 
CIDRa1.1 <- read.csv("CIDRa1.1_identical_sequences.csv") 
CIDRa1.8 <- read.csv("CIDRa1.8_identical_sequences.csv") 
DBLb12 <- read.csv("DBLb12_identical_sequences.csv") 
DBLg4 <- read.csv("DBLg4_identical_sequences.csv") 
DBLg6 <- read.csv("DBLg6_identical_sequences.csv") 
DBLd1 <- read.csv("DBLd1_identical_sequences.csv")
CIDRb1 <- read.csv("CIDRb1_identical_sequences.csv")

# Function to label and clean datasets
clean_and_label <- function(df, domain_label) {
  df %>%
    mutate(Domain = domain_label) %>%
    select(-c(X, filename))
}

# Apply function to all datasets
DBLa2 <- clean_and_label(DBLa2, 'DBLa2') 
# DBLa2 has had 3DMDS also, so more cluttered and needs its own pruning
DBLa2 <- DBLa2 %>%
  select(-Unnamed..0.1,-Unnamed..0,-mds_3d, -mds_3d_x,-mds_3d_y,-mds_3d_z)
DBLa1.2 <- clean_and_label(DBLa1.2, 'DBLa1.2')
CIDRa1.1 <- clean_and_label(CIDRa1.1, 'CIDRa1.1')
CIDRa1.8 <- clean_and_label(CIDRa1.8, 'CIDRa1.8')
DBLb12 <- clean_and_label(DBLb12, 'DBLb12')
DBLg4 <- clean_and_label(DBLg4, 'DBLg4')
DBLg6 <- clean_and_label(DBLg6, 'DBLg6')
DBLd1 <- clean_and_label(DBLd1, 'DBLd1')
CIDRb1 <- clean_and_label(CIDRb1, 'CIDRb1')

# Combine all datasets into one
DC8_domain_combined_data <- bind_rows(DBLa2, DBLa1.2, CIDRa1.1, CIDRa1.8, DBLb12, DBLg4, DBLg6, DBLd1, CIDRb1) 


colnames(DC8_domain_combined_data)[colnames(DC8_domain_combined_data) == "count"] <- "identical_sequences"

# Calculate the number of identical sequences for each group
identical_DC_distribution <- DC8_domain_combined_data %>%
  group_by(Domain, identical_sequences) %>%
  summarise(identical_sequences_count = n(), .groups = "drop") %>%  # Calculate count of identical sequences
  mutate(Actual_Number_of_identical_sequences = identical_sequences_count / identical_sequences)  # Calculate actual sequences


# Calculate sample size for each domain
samplesize <- DC8_domain_combined_data %>%
  group_by(Domain) %>%
  summarise(sample_size = n()) %>%
  mutate(Domain_label = paste(Domain, "\nn =", sample_size))  # Create a new label with sample size

# Plot the histogram and include sample size in facet labels
ggplot(identical_DC_distribution, aes(x = identical_sequences, weight = Actual_Number_of_identical_sequences)) +
  geom_histogram(binwidth = 5, fill = "blue", alpha = 0.7, color = "black") +
  scale_y_log10() +  # Log scale for better visibility
  facet_wrap(~ Domain, labeller = labeller(Domain = setNames(samplesize$Domain_label, samplesize$Domain))) +  # Use custom labels with sample sizes
  labs(
    title = "Distribution of Identical Sequences by Domain",
    x = "Identical Sequence Copies (log)",
    y = "Number of sequences"
  ) +
  theme_minimal()  # Apply minimal theme


ggplot(identical_DC_distribution, aes(x = identical_sequences, y = Actual_Number_of_identical_sequences, fill = Domain)) +
  geom_col(alpha = 0.7, color = "black", position = "dodge") + # Dodge for better visibility
  scale_fill_brewer(palette = "Dark2") +  # Use a colorblind-friendly palette
  labs(
    title = "Distribution of Identical Sequences Across Domains",
    x = "Number of Identical Sequences",
    y = "Actual Count of Identical Sequences",
    fill = "Domain"
  ) +
  theme_minimal()

ggplot(identical_DC_distribution, aes(x = Actual_Number_of_identical_sequences,weight=identical_sequences, fill = Domain)) +
  geom_histogram(binwidth = 1, position = "dodge", color = "black", alpha = 0.7) +  
  scale_x_log10() +  # Log scale for better visualization
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Distribution of Identical Sequences Across Domains",
    x = "Number of Identical Sequences (log scale)",
    y = "Count of Identical Sequences",
    fill = "Domain"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 Want to look at subtype distribution #######
 
 library(dplyr)
 
 # trying whatever the hell this is
 # Ideally, want to see how many identical domains occur on average per subtye and per domain. e.g CIDRa1.1 having identical sequences across multiple PfEMP1 structures would mean it is significantly more conserved than others
 
 
 identical_DC_distribution_subtype <- DC8_domain_combined_data %>%
   group_by(sequence,updated_arc,updated_Subtype, Domain, identical_sequences) %>%
   summarise(identical_sequences_count = n(), .groups = "drop") %>%  # Calculate count of identical sequences
   mutate(Actual_Number_of_identical_sequences = identical_sequences_count / identical_sequences) %>%
   group_by(updated_Subtype, Domain) %>%
   summarise(Identical_Domain_Count = n())# Calculate actual sequences 
 
 identical_DC_distribution_domain <- DC8_domain_combined_data %>%
   group_by(sequence,updated_arc,updated_Subtype, Domain, identical_sequences) %>%
   summarise(identical_sequences_count = n(), .groups = "drop") %>%  # Calculate count of identical sequences
   mutate(Actual_Number_of_identical_sequences = identical_sequences_count / identical_sequences) %>%
   group_by(updated_arc,Domain,updated_Subtype) %>%
   summarise(Identical_Domain_Count = n()) %>%
   group_by(Domain,updated_Subtype)
 
 
 ggplot(identical_DC_distribution_subtype, aes(x = updated_Subtype, y = Identical_Domain_Count)) +
   geom_boxplot(alpha = 0.7, fill = "lightgray") +  # Overall boxplot in gray
   geom_jitter(aes(color = Domain), shape = 16, position = position_jitter(0.0), alpha = 0.7) +  # Jitter points colored by Domain
   theme_minimal() +
   labs(title = "Distribution of Identical Sequences Across Domains and Subtypes",
        x = "Subtype",
        y = "Number of Identical Sequences",
        color = "Domain") +  # Legend title should be "color", not "fill"
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 ggplot(identical_DC_distribution_domain, aes(x = Domain, y = Identical_Domain_Count)) +
   geom_violin(aes(fill = Domain), alpha = 0.5) +  # Violin plot with slight transparency
  # Jittered points colored by subtype
   theme_minimal() +
   labs(title = "Distribution of Identical Domains across Subtypes",
        x = "Domain",
        y = "Number of Identical Sequences",
        color = "Subtype",
        fill = "Domain") +  # Adjust legend titles
   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
 
 library(patchwork)  # Load patchwork for combining plots
 
 # First plot: Identical Sequences Across Subtypes
 p1 <- ggplot(identical_DC_distribution_subtype, aes(x = updated_Subtype, y = Identical_Domain_Count)) +
   geom_boxplot(alpha = 0.7, fill = "lightgray") +  
   geom_jitter(aes(color = Domain), shape = 16, position = position_jitter(0.1), alpha = 0.7) +  
   theme_minimal() +
   labs(title = "Identical Sequences Across Subtypes",
        x = "Subtype",
        y = "Number of Identical Sequences",
        color = "Domain") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 # Second plot: Identical Domains Across Subtypes
 p2 <- ggplot(identical_DC_distribution_domain, aes(x = Domain, y = Identical_Domain_Count)) +
   geom_boxplot(alpha = 0.7, fill = "lightgray") +  
   geom_jitter( shape = 16, position = position_jitter(0.1), alpha = 0.7) +  
   theme_minimal() +
   labs(title = "Identical Domains Across Subtypes",
        x = "Domain",
        y = "Number of Identical Sequences",
        color = "Subtype") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 # Combine Plots
 p1 + p2 + plot_layout(ncol = 2)  # Side-by-side layout
 
 write_csv(identical_DC_distribution_domain,"~/Downloads/DC identical Domain.csv", row.names = FALSE)
 
 write.csv(identical_DC_distribution_subtype,"~/Downloads/DC identical subtype.csv", row.names = FALSE)
 

   


 
# checking by each unique sequence, how many copies does it have in the dataset, and then how many unique PfEMP1 structures is this shared across???

identical_DC_distribution_structure <- DC8_domain_combined_data %>%
  group_by(sequence, Domain, updated_arc, updated_Subtype) %>%
  summarise(identical_sequences_count = n(), .groups = "drop") %>%  # Count identical sequences within each updated_arc
  group_by(sequence, Domain) %>%  # Now group by Sequence and Domain
  summarise(Number_of_Copies = sum(identical_sequences_count),  # Total occurrences of identical sequences across updated_arcs
            Distinct_Updated_Arcs = n_distinct(updated_arc),  # Number of distinct structures the identical sequence is found in
            Distinct_Updated_Subtypes = n_distinct(updated_Subtype)) %>%  # Number of distinct subtypes the identical sequence is found in
  arrange(sequence, Domain)  # Optional: Order for easier inspection

# Violin plot for Number of Copies
ggplot(identical_DC_distribution_structure, aes(x = Domain, y = Number_of_Copies)) +
  geom_violin(fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Violin Plot: Number of Copies of Identical Sequences Across Domains",
       x = "Domain", 
       y = "Number of Copies") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Summarize the number of identical sequences per domain
identical_counts <- DC8_domain_combined_data %>%
  group_by(Domain) %>%
  summarise(Number_of_Identical_Sequences = n())

# Create the violin plot with labels for identical sequences
ggplot(identical_DC_distribution_structure, aes(x = Domain, y = Distinct_Updated_Arcs)) +
  geom_violin(fill = "lightgreen", color = "black") +
  geom_text(data = identical_counts, aes(x = Domain, y = 0, label = paste("n=", Number_of_Identical_Sequences)),
            color = "black", size = 3, vjust = -0.5) +  # Adjust vjust to position text
  theme_minimal() +
  labs(title = "Distinct Updated Arcs (Unique PfEMP1s) for Each Domain",
       x = "Domain", 
       y = "Unique PfEMP1s") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
   
