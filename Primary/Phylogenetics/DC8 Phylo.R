
# Install ggtree
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtreeExtra")

# bootstrap extractor

extract_bootstrap_values <- function(tree_df) {
  # Extract internal nodes with bootstrap values
  df_tree_bootstrap <- tree_df %>%
    mutate(bootstrap = as.numeric(label)) %>%
    filter(!is.na(bootstrap))  # Keep only rows with valid bootstrap values
  
  # Return just the bootstrap values
  return(df_tree_bootstrap$bootstrap)
}

# Consistent color scheme
custom_colors <- c(
  'IT4var20-like (EPCR)' = '#1f77b4',
  'Novel Subtype A (EPCR)' = '#ff7f0e',
  'PFD0020c-like (EPCR + gC1qR)' = '#2ca02c',
  'TIM180var1-like (EPCR + Rosetting)' = '#d62728',
  'Novel Subtype B (EPCR)' = '#98df8a',
  'Novel Subtype C (EPCR + Rosetting)' = '#9467bd',
  'Other' = 'gray70'
)


######
#load packages

library(tidyverse) #select() #separate() #mutate() #%>% #gather()
library(ape) #read.tree()
library(ggtree) #ggtree()
library(ggtreeExtra)
library(treeio) #read.newick()
library(tidytree)
library(ggnewscale)

#set directory

setwd("~/Documents/Dissertation/Data/Unprocessed Data/Phylogenetics")
getwd()

# make tree for full length sequences -----------------------------------------------------------------

### load sequences

#upload tree
pfemp1_tree <-
  read.newick(
    "DC8_filtered_sequences.noATS.aligned.fasta.treefile")

#upload metadata
geo <- read.csv("DC8_df_filtered.csv") 

geo <- geo %>%
  mutate(Subtype = case_when(
    Subtype == "DBLa2-CIDRa1.8 head with DBLg4" ~ "TIM180var1-like",
    Subtype == "DBLa2-CIDRa1.1 head with DBLg6" ~ "IT4var20-like",
    Subtype == "DBLa1.2-CIDRa1.1 head with DBLg6" ~ "PFD0020c-like",
    Subtype == "DBLa1.2-CIDRa1.1 head with DBLg4" ~ "DBLa1.2-CIDRa1.1-DBLb12-DBLg4",
    Subtype == "DBLa2-CIDRa1.1 head with DBLg4" ~ "DBLa2-CIDRa1.1-DBLb12-DBLg4",
    Subtype == "DBLa2-CIDRa1.8 head with DBLg6" ~ "DBLa2-CIDRa1.8-DBLb12-DBLg6",
    TRUE ~ "Other"
  ))

#create a dataframe based on the tree (giving tip labels, branch lengths, nodes, bootstrap values)
df_tree <- tidytree::as_tibble(pfemp1_tree) 

#make a few changes: sort by label, remove bootstrap values 
df_tree_sorted <-
  df_tree %>%
  arrange(label) %>%
  filter(is.na(as.numeric(label)) & label != "") #might give a warning but that's okay

# Hello future me put this for all domains please

DC8_bootstrap_df <- data.frame(DC8_bootstrap = extract_bootstrap_values(df_tree))

# Summary stats for DC8!

DC8_bootstrap_summary <- DC8_bootstrap_df %>%
  filter(!is.na(DC8_bootstrap)) %>%   # Exclude NA values if any
  summarise(
    Mean = mean(DC8_bootstrap),         # Mean of bootstrap values
    Median = median(DC8_bootstrap),     # Median of bootstrap values
    Standard_Deviation = sd(DC8_bootstrap),  # Standard deviation of bootstrap values
    Min = min(DC8_bootstrap),           # Minimum value
    Max = max(DC8_bootstrap),           # Maximum value
    Count = n(),                    # Count of bootstrap values
    NA_Count = sum(is.na(DC8_bootstrap)) # Count of NAs (just in case)
  )

extra_labels <- setdiff(df_tree_sorted$label, geo$sequence_id)
print(extra_labels)
df_tree_sorted <- df_tree_sorted %>%
  filter(!(label %in% extra_labels))

# For whatever reason here TIM80var1 was an extra label sequence. Had to delete here unfortunately to get it to work

#pull out metadata columns of interest + add labels to metadata to compare whether metadata is the same. Here jessie is choosing country, continent, subtype etc
geo_selected <- 
  geo %>%
  select(sequence_id,Continent,Subtype,Country) %>%
  arrange(sequence_id)%>%
  mutate(label = df_tree_sorted$label) 

#test if Sample_ID and label match -> NEED TO MAKE SURE LABELS MATCH OR YOU WILL HAVE a mess
test_if_labels_match <- 
  geo_selected %>%
  select(sequence_id, label) %>%
  mutate(test = ifelse(sequence_id == label, "match", "mismatch")) %>%
  arrange(desc(test))
# Success!


#modify mismatches (either manually or in R)
geo_selected <- 
  geo_selected %>%
  mutate(Sample_ID = ifelse(Sample_ID == "11019var1_exon1_aa","11019var1_aa", Sample_ID))
#then go back and check test_if_labels_match again

# Above isn't needed! Labels match completely

df_tree_geo <-
  df_tree_sorted %>%
  merge(geo_selected, by = "label", all = TRUE)

# Now to get to the juicy tree
#######




#generate tree
phylo <- 
  ggtree(pfemp1_tree, layout = "rectangular")

phylo

#look at node numbers
phylo + geom_text(aes(label = node), size = 2) # add node numbers

#make adjustments and set midpoint - midrooting tree
# also can look at subsets using this! 
adj_phylo <- 
  viewClade(phylo, node=1080, xmax_adjust = 2) 

adj_phylo

adj_phylo + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12) +
  geom_tiplab(align = T, linesize = 0.05, size = 0.5) +
  
  ## Heatmap for Subtype
  new_scale_fill() +
  geom_fruit(
    data = df_tree_geo,
    geom = geom_tile,
    mapping = aes(
      y = label, 
      x = "",
      fill = Subtype,
      width = 0.4
    ),
    offset = 0.95,
    pwidth = 2.0,
    grid.params = list(alpha = 0),
    position = position_identityx(hexpand = 2.1)
  ) +
  scale_fill_manual(
    values = c(
      'IT4var20-like' = '#1f77b4',
      'DBLa2-CIDRa1.1-DBLb12-DBLg4' = '#ff7f0e',
      'PFD0020c-like' = '#2ca02c',
      'TIM180var1-like' = '#d62728',
      'DBLa1.2-CIDRa1.1-DBLb12-DBLg4' = '#98df8a',
      'DBLa2-CIDRa1.8-DBLb12-DBLg6' = '#9467bd',
      'Other' = 'gray70'  # default for anything else
    )
  ) +
  labs(fill = "Subtype") +
  theme(
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  ) +
  
  ## Heatmap for Continent
  new_scale_fill() +
  geom_fruit(
    data = df_tree_geo,
    geom = geom_tile,
    mapping = aes(
      y = label, 
      x = "",
      fill = Continent,
      width = 0.4
    ),
    offset = 1.0,
    pwidth = 2.0,
    grid.params = list(alpha = 0),
    position = position_identityx(hexpand = 2.1)
  ) +
  scale_fill_manual(
    breaks = c("Africa", "Asia", "S.America"),
    labels = c("Africa", "Asia", "South America"),
    values = c("Africa" = "#8679b5", "Asia" = "#c781b1", "S.America" = "#baace6")
  ) +
  labs(fill = "Continent") +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15)
  )


# Continent v Country overall
adj_phylo + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ #add scale and adjust appearance
  #geom_nodelab(size = 1.5, hjust = .05)+ #add bootstrap values
  geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ # add sample names (can remove this to reduce file size if needed)
  
  ##heatmaps
  
  new_scale_fill() + # Reset Colours
  
  geom_fruit(
    data=df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Subtype, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.95, #how far apart is the heatmap from the tree?
    pwidth = 2.0, #how narrow/wide is the heatmap?
    
    #leave this as it is (keeps the heatmap in position)
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  
  #customise legend
  labs(fill = "Subtype") + #name your legend
  theme(
    legend.title = element_text(size = 18), #change font size
    legend.text = element_text(size = 15)) + #change font size
  
  
  
  #add another layer
  new_scale_fill() +
  
  geom_fruit(
    data=df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Continent, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 1.0, #how far apart is the heatmap from the tree?
    pwidth = 2.0, #how narrow/wide is the heatmap?
    
    # hide grid lines in rectangular grids
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1))+
  
  #customise legend
  labs(fill = "Continent") + #name your legend
  theme(
    legend.title = element_text(size = 18), #change font size
    legend.text = element_text(size = 15))+ #change font size
  scale_fill_manual(
    breaks = c("Africa","Asia","S.America"), #indicate ALL groups with the exact names - for ordering
    labels = c("Africa","Asia","South America"), #rename groups (optional)
    values = c("#8679b5","#c781b1","#baace6"))

#####
# Continent v country

adj_phylo + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ #add scale and adjust appearance
  #geom_nodelab(size = 1.5, hjust = .05)+ #add bootstrap values
  geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ # add sample names (can remove this to reduce file size if needed)
  
  ##heatmaps
  
  new_scale_fill() + # Reset Colours
  
  geom_fruit(
    data=df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Country, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.95, #how far apart is the heatmap from the tree?
    pwidth = 2.0, #how narrow/wide is the heatmap?
    
    #leave this as it is (keeps the heatmap in position)
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  
  #customise legend
  labs(fill = "Country") + #name your legend
  theme(
    legend.title = element_text(size = 18), #change font size
    legend.text = element_text(size = 15)) + #change font size
  
  
  
  #add another layer
  new_scale_fill() +
  
  geom_fruit(
    data=df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Continent, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 1.0, #how far apart is the heatmap from the tree?
    pwidth = 2.0, #how narrow/wide is the heatmap?
    
    # hide grid lines in rectangular grids
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1))+
  
  #customise legend
  labs(fill = "Continent") + #name your legend
  theme(
    legend.title = element_text(size = 18), #change font size
    legend.text = element_text(size = 15))+ #change font size
  scale_fill_manual(
    breaks = c("Africa","Asia","S.America"), #indicate ALL groups with the exact names - for ordering
    labels = c("Africa","Asia","South America"), #rename groups (optional)
    values = c("#8679b5","#c781b1","#baace6"))


###### 
#Now want to specfically look at the outgrouping at the top!


adj_phylo_top <- 
  viewClade(phylo, node=1268, xmax_adjust = 2)
  
adj_phylo_top + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ #add scale and adjust appearance
  #geom_nodelab(size = 1.5, hjust = .05)+ #add bootstrap values
  geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ # add sample names (can remove this to reduce file size if needed)
  
  ##heatmaps
  
  new_scale_fill() + # Reset Colours
  
  geom_fruit(
    data=df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Country, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.95, #how far apart is the heatmap from the tree?
    pwidth = 2.0, #how narrow/wide is the heatmap?
    
    #leave this as it is (keeps the heatmap in position)
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  
  #customise legend
  labs(fill = "Subtype") + #name your legend
  theme(
    legend.title = element_text(size = 18), #change font size
    legend.text = element_text(size = 15)) + #change font size
  
  
  
  #add another layer
  new_scale_fill() +
  
  geom_fruit(
    data=df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Continent, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 1.0, #how far apart is the heatmap from the tree?
    pwidth = 2.0, #how narrow/wide is the heatmap?
    
    # hide grid lines in rectangular grids
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1))+
  
  #customise legend
  labs(fill = "Continent") + #name your legend
  theme(
    legend.title = element_text(size = 18), #change font size
    legend.text = element_text(size = 15))+ #change font size
  scale_fill_manual(
    breaks = c("Africa","Asia","S.America"), #indicate ALL groups with the exact names - for ordering
    labels = c("Africa","Asia","South America"), #rename groups (optional)
    values = c("#8679b5","#c781b1","#baace6"))


DBla2 ########

setwd("/Users/Jakub/Documents/Dissertation/Data/Unprocessed Data/Phylogenetics/Domain Sequences Fasta Files")

DBla2_tree <-
  read.newick(
    "DC8_DBla2_wholeseq.aligned.fasta.treefile")

DBla2geo <- read.csv("DC8_DBLa2_3D_metadata.csv") 

#create a dataframe based on the tree (giving tip labels, branch lengths, nodes, bootstrap values)
DBla2_df_tree <- tidytree::as_tibble(DBla2_tree) 
DBla2_df_tree_sorted <-
  DBla2_df_tree %>%
  arrange(label) %>%
  filter(is.na(as.numeric(label)) & label != "")

DBLa2geo_selected <- 
  DBla2geo %>%
  select(sequence_id, Continent, updated_Subtype, Country) %>%
  arrange(sequence_id) %>%
  mutate(label = DBla2_df_tree_sorted$label) 

test_if_labels_match <- 
  DBLa2geo_selected %>%
  select(sequence_id, label) %>%
  mutate(test = ifelse(sequence_id == label, "match", "mismatch")) %>%
  arrange(desc(test))
# Splendid

DBLa2_df_tree_geo <-
  DBla2_df_tree_sorted %>%
  merge(DBLa2geo_selected, by = "label", all = TRUE)


phylo <- 
  ggtree(DBla2_tree, layout = "rectangular")

phylo

#look at node numbers
phylo + geom_text(aes(label = node), size = 2) # add node numbers

#make adjustments and set midpoint - midrooting tree
# also can look at subsets using this! 
adj_phylo <- 
  viewClade(phylo, node=1205, xmax_adjust = 0.4) 

adj_phylo

adj_phylo + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ #add scale and adjust appearance
  #geom_nodelab(size = 1.5, hjust = .05)+ #add bootstrap values
  geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ # add sample names (can remove this to reduce file size if needed)
  
  ##heatmaps
  
  new_scale_fill() + # Reset Colours
  
  geom_fruit(
    data=DBLa2_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = updated_Subtype, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.65, #how far apart is the heatmap from the tree?
    pwidth = 1.0, #how narrow/wide is the heatmap?
    
    #leave this as it is (keeps the heatmap in position)
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  
  #customise legend
  labs(fill = "Subtype") + #name your legend
  theme(
    legend.title = element_text(size = 5), #change font size
    legend.text = element_text(size = 5)) + #change font size
  
  
  
  #add another layer
  new_scale_fill() +
  
  geom_fruit(
    data=DBLa2_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Continent, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.5, #how far apart is the heatmap from the tree?
    pwidth = 0.3, #how narrow/wide is the heatmap?
    
    # hide grid lines in rectangular grids
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1))+
  
  #customise legend
  labs(fill = "Continent") + #name your legend
  theme(
    legend.title = element_text(size = 18), #change font size
    legend.text = element_text(size = 15))+ #change font size
  scale_fill_manual(
    breaks = c("Africa","Asia","S.America"), #indicate ALL groups with the exact names - for ordering
    labels = c("Africa","Asia","South America"), #rename groups (optional)
    values = c("#8679b5","#c781b1","#baace6"))


CIDRa1.1 #####


setwd("~/Documents/Dissertation/Data/Unprocessed Data/Phylogenetics/Domain Sequences Fasta Files")

CIDRa1.1_tree <-
  read.newick(
    "DC8_CIDRa1.1_wholeseq.aligned.fasta.treefile")

CIDRa1.1geo <- read.csv("CIDRa1.1_geo.csv") %>%
  filter(filename != 	
           'normalised.fulldataset.all_DC8_sequences_domains_subclasses.fasta')

CIDRa1.1geo <- CIDRa1.1geo %>%
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

#create a dataframe based on the tree (giving tip labels, branch lengths, nodes, bootstrap values)
CIDRa1.1_df_tree <- tidytree::as_tibble(CIDRa1.1_tree) 
CIDRa1.1_df_tree_sorted <-
  CIDRa1.1_df_tree %>%
  arrange(label) %>%
  filter(is.na(as.numeric(label)) & label != "")

CIDRa1.1_bootstrap_df <- data.frame(CIDRa1.1_bootstrap = extract_bootstrap_values(CIDRa1.1_df_tree))

# Summary stats for DC8!

CIDRa1.1bootstrap_summary <- CIDRa1.1_bootstrap_df %>%
  filter(!is.na(CIDRa1.1_bootstrap)) %>%   # Exclude NA values if any
  summarise(
    Mean = mean(CIDRa1.1_bootstrap),         # Mean of bootstrap values
    Median = median(CIDRa1.1_bootstrap),     # Median of bootstrap values
    Standard_Deviation = sd(CIDRa1.1_bootstrap),  # Standard deviation of bootstrap values
    Min = min(CIDRa1.1_bootstrap),           # Minimum value
    Max = max(CIDRa1.1_bootstrap),           # Maximum value
    Count = n(),                    # Count of bootstrap values
    NA_Count = sum(is.na(CIDRa1.1_bootstrap)) # Count of NAs (just in case)
  )

CIDRa1.1geo_selected <- 
  CIDRa1.1geo %>%
  select(sequence_id, Continent,Subtype, Country) %>%
  arrange(sequence_id) %>%
  mutate(label = CIDRa1.1_df_tree_sorted$label) 

test_if_labels_match <- 
  CIDRa1.1geo_selected %>%
  select(sequence_id, label) %>%
  mutate(test = ifelse(sequence_id == label, "match", "mismatch")) %>%
  arrange(desc(test))
# Splendid

CIDRa1.1_df_tree_geo <-
  CIDRa1.1_df_tree_sorted %>%
  merge(CIDRa1.1geo_selected, by = "label", all = TRUE)


phylo <- 
  ggtree(CIDRa1.1_tree, layout = "circular")

phylo

#look at node numbers
phylo + geom_text(aes(label = node), size = 2) # add node numbers

#make adjustments and set midpoint - midrooting tree
# also can look at subsets using this! 
adj_phylo <- 
  viewClade(phylo, node=799, xmax_adjust = 0.4) 

adj_phylo

# Add labels. Very disparate suprisingly
  adj_phylo + 
    geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ 
    geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ 
    
    ## Heatmaps ##
    new_scale_fill() + 
    geom_fruit(
      data=CIDRa1.1_df_tree_geo,
      geom=geom_tile,
      mapping = aes(y=label, x="1", fill = Subtype),
      offset = 0.16, 
      pwidth = 0.2755, 
      grid.params=list(alpha = 0),
      position = position_identityx(hexpand = 2.1)) + labs(fill='Subtype') +  scale_fill_manual(values = custom_colors)
    

  
  

  
CIDRa1.8 ######


CIDRa1.8_tree <-
  read.newick(
    "DC8_CIDRa1.8_wholeseq.aligned.fasta.treefile")

CIDRa1.8geo <- read.csv("CIDRa1.8_geo.csv") %>%
  filter(filename != 	
           'normalised.fulldataset.all_DC8_sequences_domains_subclasses.fasta')
CIDRa1.8geo <- CIDRa1.8geo %>%
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
#create a dataframe based on the tree (giving tip labels, branch lengths, nodes, bootstrap values)
CIDRa1.8_df_tree <- tidytree::as_tibble(CIDRa1.8_tree) 
CIDRa1.8_df_tree_sorted <-
  CIDRa1.8_df_tree %>%
  arrange(label) %>%
  filter(is.na(as.numeric(label)) & label != "")

CIDRa1.8_bootstrap_df <- data.frame(CIDRa1.8_bootstrap = extract_bootstrap_values(CIDRa1.8_df_tree))

# Summary stats for DC8!

CIDRa1.8bootstrap_summary <- CIDRa1.8_bootstrap_df %>%
  filter(!is.na(CIDRa1.8_bootstrap)) %>%   # Exclude NA values if any
  summarise(
    Mean = mean(CIDRa1.8_bootstrap),         # Mean of bootstrap values
    Median = median(CIDRa1.8_bootstrap),     # Median of bootstrap values
    Standard_Deviation = sd(CIDRa1.8_bootstrap),  # Standard deviation of bootstrap values
    Min = min(CIDRa1.8_bootstrap),           # Minimum value
    Max = max(CIDRa1.8_bootstrap),           # Maximum value
    Count = n(),                    # Count of bootstrap values
    NA_Count = sum(is.na(CIDRa1.8_bootstrap)) # Count of NAs (just in case)
  )

CIDRa1.8geo_selected <- 
  CIDRa1.8geo %>%
  select(sequence_id, Continent, Subtype, Country) %>%
  arrange(sequence_id) %>%
  mutate(label = CIDRa1.8_df_tree_sorted$label) 

test_if_labels_match <- 
  CIDRa1.8geo_selected %>%
  select(sequence_id, label) %>%
  mutate(test = ifelse(sequence_id == label, "match", "mismatch")) %>%
  arrange(desc(test))
# Splendid

CIDRa1.8_df_tree_geo <-
  CIDRa1.8_df_tree_sorted %>%
  merge(CIDRa1.8geo_selected, by = "label", all = TRUE)

CIDRa1.8phylo <- 
  ggtree(CIDRa1.8_tree, layout = "rectangular")

CIDRa1.8phylo

CIDRa1.8phylo_circ <- ggtree(CIDRa1.8_tree, layout = 'circular')

CIDRa1.8phylo_circ
#look at node numbers
phylo + geom_text(aes(label = node), size = 2) # add node numbers

#make adjustments and set midpoint - midrooting tree
# also can look at subsets using this! 
CIDRa1.8_adj_phylo <- 
  viewClade(CIDRa1.8phylo, node = 477, xmax_adjust = 1.5)

CIDRa1.8_adj_phylo

CIDRa1.8_adj_phylo + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ #add scale and adjust appearance
  #geom_nodelab(size = 1.5, hjust = .05)+ #add bootstrap values
  geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ # add sample names (can remove this to reduce file size if needed)
  
  ##heatmaps
  
  new_scale_fill() + # Reset Colours
  
  geom_fruit(
    data=CIDRa1.8_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Subtype, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.52, #how far apart is the heatmap from the tree?
    pwidth = 0.5, #how narrow/wide is the heatmap?
    
    #leave this as it is (keeps the heatmap in position)
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  
  #customise legend
  labs(fill = "Subtype") + scale_fill_manual(values = custom_colors)  #name your legend
  #change font size






DBLb12########

DBlb12_tree <-
  read.newick(
    "DC8_DBlb12_wholeseq.aligned.fasta.treefile")

DBlb12geo <- read.csv("DBLb12_geo.csv") %>%
  filter(filename!='normalised.fulldataset.all_DC8_sequences_domains_subclasses.fasta')

#create a dataframe based on the tree (giving tip labels, branch lengths, nodes, bootstrap values)
DBlb12_df_tree <- tidytree::as_tibble(DBlb12_tree) 
DBlb12_df_tree_sorted <-
  DBlb12_df_tree %>%
  arrange(label) %>%
  filter(is.na(as.numeric(label)) & label != "")

DBLb12geo_selected <- 
  DBlb12geo %>%
  select(sequence_id, Continent, updated_Subtype, Country) %>%
  arrange(sequence_id) %>%
  mutate(label = DBlb12_df_tree_sorted$label) 

test_if_labels_match <- 
  DBLb12geo_selected %>%
  select(sequence_id, label) %>%
  mutate(test = ifelse(sequence_id == label, "match", "mismatch")) %>%
  arrange(desc(test))
# Splendid

DBLb12_df_tree_geo <-
  DBlb12_df_tree_sorted %>%
  merge(DBLb12geo_selected, by = "label", all = TRUE)


DBLb12_phylo <- 
  ggtree(DBlb12_tree, layout = "rectangular")

DBLb12_phylo

#look at node numbers
phylo + geom_text(aes(label = node), size = 2) # add node numbers

#make adjustments and set midpoint - midrooting tree
# also can look at subsets using this! 
DBLb12_adj_phylo <- 
  viewClade(DBLb12_phylo, node=1275, xmax_adjust = 1.5) 

DBLb12_adj_phylo

DBLb12_adj_phylo + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ #add scale and adjust appearance
  #geom_nodelab(size = 1.5, hjust = .05)+ #add bootstrap values
  geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ # add sample names (can remove this to reduce file size if needed)
  
  ##heatmaps
  
  new_scale_fill() + # Reset Colours
  
  geom_fruit(
    data=DBLb12_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = updated_Subtype, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.2, #how far apart is the heatmap from the tree?
    pwidth = 0.5, #how narrow/wide is the heatmap?
    
    #leave this as it is (keeps the heatmap in position)
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  
  #customise legend
  labs(fill = "Subtype") + #name your legend
  theme(
    legend.title = element_text(size = 5), #change font size
    legend.text = element_text(size = 5)) + #change font size
  
  
  
  #add another layer
  new_scale_fill() +
  
  geom_fruit(
    data=DBLb12_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Continent, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.3, #how far apart is the heatmap from the tree?
    pwidth = 0.3, #how narrow/wide is the heatmap?
    
    # hide grid lines in rectangular grids
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1))+
  
  #customise legend
  labs(fill = "Continent") + #name your legend
  theme(
    legend.title = element_text(size = 10), #change font size
    legend.text = element_text(size = 5))+ #change font size
  scale_fill_manual(
    breaks = c("Africa","Asia","S.America"), #indicate ALL groups with the exact names - for ordering
    labels = c("Africa","Asia","South America"), #rename groups (optional)
    values = c("#8679b5","#c781b1","#baace6"))

DBLg4#######


DBlg4_tree <-
  read.newick(
    "DC8_DBlg4_wholeseq.aligned.fasta.treefile")

DBlg4geo <- read.csv("DBLg4_geo.csv") %>%
  filter(filename!='normalised.fulldataset.all_DC8_sequences_domains_subclasses.fasta')

DBlg4geo <- DBlg4geo %>%
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

#create a dataframe based on the tree (giving tip labels, branch lengths, nodes, bootstrap values)
DBlg4_df_tree <- tidytree::as_tibble(DBlg4_tree) 
DBlg4_df_tree_sorted <-
  DBlg4_df_tree %>%
  arrange(label) %>%
  filter(is.na(as.numeric(label)) & label != "")

DBLg4geo_selected <- 
  DBlg4geo %>%
  select(sequence_id, Continent, Subtype, Country) %>%
  arrange(sequence_id) %>%
  mutate(label = DBlg4_df_tree_sorted$label) 

test_if_labels_match <- 
  DBLg4geo_selected %>%
  select(sequence_id, label) %>%
  mutate(test = ifelse(sequence_id == label, "match", "mismatch")) %>%
  arrange(desc(test))
# Splendid

DBLg4_df_tree_geo <-
  DBlg4_df_tree_sorted %>%
  merge(DBLg4geo_selected, by = "label", all = TRUE)


DBLg4_phylo <- 
  ggtree(DBlg4_tree, layout = "rectangular")

DBLg4_phylo

#look at node numbers
phylo + geom_text(aes(label = node), size = 2) # add node numbers

#make adjustments and set midpoint - midrooting tree
# also can look at subsets using this! 
DBLg4_adj_phylo <- 
  viewClade(DBLg4_phylo, node=822, xmax_adjust = 1) 

DBLg4_adj_phylo

DBLg4_adj_phylo + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ #add scale and adjust appearance
  #geom_nodelab(size = 1.5, hjust = .05)+ #add bootstrap values
  geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ # add sample names (can remove this to reduce file size if needed)
  
  ##heatmaps
  
  new_scale_fill() + # Reset Colours
  
  geom_fruit(
    data=DBLg4_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Subtype, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.4, #how far apart is the heatmap from the tree?
    pwidth = 0.5, #how narrow/wide is the heatmap?
    
    #leave this as it is (keeps the heatmap in position)
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  
  #customise legend
  labs(fill = "Subtype") + #name your legend
  theme(
    legend.title = element_text(size = 5), #change font size
    legend.text = element_text(size = 5))+ scale_fill_manual(values = custom_colors) #change font size
  
  
  
  #add another layer
  new_scale_fill() +
  
  geom_fruit(
    data=DBLg4_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Continent, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.5, #how far apart is the heatmap from the tree?
    pwidth = 0.3, #how narrow/wide is the heatmap?
    
    # hide grid lines in rectangular grids
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1))+
  
  #customise legend
  labs(fill = "Continent") + #name your legend
  theme(
    legend.title = element_text(size = 10), #change font size
    legend.text = element_text(size = 5))+ #change font size
  scale_fill_manual(
    breaks = c("Africa","Asia","S.America"), #indicate ALL groups with the exact names - for ordering
    labels = c("Africa","Asia","South America"), #rename groups (optional)
    values = c("#8679b5","#c781b1","#baace6"))

DBLg6#######


DBlg6_tree <-
  read.newick(
    "DC8_DBlg6_wholeseq.aligned.fasta.treefile")

DBlg6geo <- read.csv("DBLg6_geo.csv") %>%
  filter(filename!='normalised.fulldataset.all_DC8_sequences_domains_subclasses.fasta')

DBlg6geo <- DBlg6geo %>%
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

#create a dataframe based on the tree (giving tip labels, branch lengths, nodes, bootstrap values)
DBlg6_df_tree <- tidytree::as_tibble(DBlg6_tree) 
DBlg6_df_tree_sorted <-
  DBlg6_df_tree %>%
  arrange(label) %>%
  filter(is.na(as.numeric(label)) & label != "")

DBLg6geo_selected <- 
  DBlg6geo %>%
  select(sequence_id, Continent, Subtype, Country) %>%
  arrange(sequence_id) %>%
  mutate(label = DBlg6_df_tree_sorted$label) 

test_if_labels_match <- 
  DBLg6geo_selected %>%
  select(sequence_id, label) %>%
  mutate(test = ifelse(sequence_id == label, "match", "mismatch")) %>%
  arrange(desc(test))
# Splendid

DBLg6_df_tree_geo <-
  DBlg6_df_tree_sorted %>%
  merge(DBLg6geo_selected, by = "label", all = TRUE)


DBLg6_phylo <- 
  ggtree(DBlg6_tree, layout = "rectangular")

DBLg6_phylo

#look at node numbers
phylo + geom_text(aes(label = node), size = 2) # add node numbers

#make adjustments and set midpoint - midrooting tree
# also can look at subsets using this! 
DBLg6_adj_phylo <- 
  viewClade(DBLg6_phylo, node=455, xmax_adjust = 1) 

DBLg6_adj_phylo

DBLg6_adj_phylo + 
  geom_treescale(x=0.2, y=300, fontsize = 20, linesize = 2.5, offset = 12)+ #add scale and adjust appearance
  #geom_nodelab(size = 1.5, hjust = .05)+ #add bootstrap values
  geom_tiplab(align = T, linesize = 0.05, size = 0.5)+ # add sample names (can remove this to reduce file size if needed)
  
  ##heatmaps
  
  new_scale_fill() + # Reset Colours
  
  geom_fruit(
    data=DBLg6_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Subtype, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.49, #how far apart is the heatmap from the tree?
    pwidth = 0.5, #how narrow/wide is the heatmap?
    
    #leave this as it is (keeps the heatmap in position)
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  
  #customise legend
  labs(fill = "Subtype") + scale_fill_manual(values = custom_colors) #name your legend
  #change font size
  
  
  
  #add another layer
  new_scale_fill() +
  
  geom_fruit(
    data=DBLg6_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Continent, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.6 , #how far apart is the heatmap from the tree?
    pwidth = 0.3, #how narrow/wide is the heatmap?
    
    # hide grid lines in rectangular grids
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1))+
  
  
  #customise legend
  labs(fill = "Continent") + #name your legend #change font size
  scale_fill_manual(
    breaks = c("Africa","Asia","S.America"), #indicate ALL groups with the exact names - for ordering
    labels = c("Africa","Asia","South America"), #rename groups (optional)
    values = c("#8679b5","#c781b1","#baace6")) + new_scale_fill() +
  
  geom_fruit(
    data=DBLg6_df_tree_geo, #set the metadata df
    geom=geom_tile,
    mapping = aes(
      y=label, 
      x="",
      fill = Country, #set variable to use for heatmap
      width = 0.4 #how narrow/wide is the heatmap?
    ),
    offset = 0.7 , #how far apart is the heatmap from the tree?
    pwidth = 0.3, #how narrow/wide is the heatmap?
    
    # hide grid lines in rectangular grids
    grid.params=list(alpha = 0),
    position = position_identityx(hexpand = 2.1)) +
  theme(
    legend.position = "right",  
    legend.box = "vertical",  
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white", color = "black"),  # Ensure rectangular legend keys
    legend.key.height = unit(0.6, "cm"),  # Control height
    legend.key.width = unit(1.2, "cm"),  # Control width
    legend.spacing.y = unit(0.3, "cm")  # Space between legend items
  ) +
  guides(
    fill = guide_legend(ncol = 3,  # Arrange legends horizontally (3 columns here, you can adjust)
                        override.aes = list(size = 5))  # Keeps rectangular layout for items
  )




Comparing topologies #######

# Extracting node distances as approximation of topologies

node_distances_CIDRa1.1 <- dist.nodes(CIDRa1.1_tree)
node_distances_CIDRa1.8 <- dist.nodes(CIDRa1.8_tree)

node_distances_DBLg4 <- dist.nodes(DBlg4_tree)
node_distances_DBLg6 <- dist.nodes(DBlg6_tree)


# Then normalising with density rather than frequency as both CIDRa1.8 and DBLg6 have less incidence
# Calculate histogram objects without plotting, using density (freq=FALSE)
hist_DBLg4 <- hist(node_distances_DBLg4, breaks = 50, plot = FALSE, freq = FALSE)
hist_DBLg6 <- hist(node_distances_DBLg6, breaks = 50, plot = FALSE, freq = FALSE)

# Set x-axis limits (0 to 4) and determine y-axis limit from density values
xlim_range <- c(0, 4)
ylim_max <- max(c(hist_DBLg4$density, hist_DBLg6$density))

# Plot the first histogram (density plot)
hist(node_distances_DBLg4,
     breaks = 50,
     freq = FALSE,                           # Plot density instead of counts
     col = rgb(0, 0, 1, 0.5),                  # semi-transparent blue
     xlab = "Node Distance",
     main = "Normalized Node Distance Distributions",
     xlim = xlim_range,
     ylim = c(0, ylim_max)) +

# Overlay the second histogram (density plot)
hist(node_distances_DBLg6,
     breaks = 50,
     freq = FALSE,                           # Plot density
     col = rgb(1, 0, 0, 0.5),                  # semi-transparent red
     add = TRUE) +

# Add legend
legend("topright",
       legend = c("DBLg4", "DBLg6"),
       fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))





# Calculate histogram objects without plotting, using density (freq=FALSE)
hist_CIDRa1.1 <- hist(node_distances_CIDRa1.1, breaks = 50, plot = FALSE, freq = FALSE)
hist_CIDRa1.8 <- hist(node_distances_CIDRa1.8, breaks = 50, plot = FALSE, freq = FALSE)

# Set x-axis limits (0 to 4) and determine y-axis limit from density values
xlim_range <- c(0, 4)
ylim_range <- c(0,1.5)

# Plot the first histogram (density plot)
hist(node_distances_CIDRa1.1,
     breaks = 50,
     freq = FALSE,                           # Plot density instead of counts
     col = (186/255, 110/255, 64/255),                 # semi-transparent blue
     xlab = "Node Distance",
     main = "Normalized Node Distance Distributions",
     xlim = xlim_range,
     ylim =  ylim_range) +

# Overlay the second histogram (density plot)
hist(node_distances_CIDRa1.8,
     breaks = 50,
     freq = FALSE,                           # Plot density
     col = (186/255, 110/255, 64/255),               # semi-transparent red
     add = TRUE) +

# Add legend
legend("topright",
       legend = c("DBLg4", "DBLg6"),
       fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))

# Define your custom color
custom_color <- rgb(196/255, 160/255, 0/255, alpha = 0.5)  # warm taupe with transparency

# Plot the first histogram (density plot)
hist(node_distances_CIDRa1.1,
     breaks = 50,
     freq = FALSE,
     col = custom_color,
     xlab = "Node Distance",
     main = "Normalized Node Distance Distributions",
     xlim = xlim_range,
     ylim = ylim_range)

# Overlay the second histogram (density plot)
hist(node_distances_CIDRa1.8,
     breaks = 50,
     freq = FALSE,
     col = rgb(112/255, 128/255, 144/255, alpha = 0.5),  # slate gray with transparency
     add = TRUE)

# Add legend
legend("topright",
       legend = c("CIDRa1.1", "CIDRa1.8"),
       fill = c(custom_color, rgb(112/255, 128/255, 144/255, alpha = 0.5)))



# Distribution of CIDRa.18 are not normal, so simple wilcox will do:
wilcox.test(node_distances_CIDRa1.1, node_distances_CIDRa1.8, exact = FALSE, alternative = "two.sided")

wilcox.test(node_distances_DBLg4, node_distances_DBLg6, exact = FALSE, alternative = "two.sided")

var2csa#####


#upload tree
var2csa_tree <-
  read.newick(
    "VAR2CSA.aligned.fasta.treefile")

#upload metadata
geo <- read.csv("vardb_var2csa.cleaned.isolates.population.csv") 

geo <- geo %>%
  filter(!(duplicated(sequence_id) & sequence_id == "PC0083-C.g148"))


#create a dataframe based on the tree (giving tip labels, branch lengths, nodes, bootstrap values)
var2csa_df_tree <- tidytree::as_tibble(var2csa_tree) 

# Step 1: Extract internal nodes (bootstrap values)
df_var2csa_bootstrap <- var2csa_df_tree %>%
  mutate(bootstrap = as.numeric(label)) %>%              # Convert labels to numeric
  filter(!is.na(bootstrap)) %>%                          # Keep only those that converted
  arrange(desc(bootstrap))                               # Optional: sort by strength

# Step 2: Clean and sort tip labels (non-numeric = sequence IDs)
df_tree_sorted <- var2csa_df_tree %>%
  arrange(label) %>%
  filter(is.na(as.numeric(label)) & label != "")         # Keep only tip labels

# Step 3: Remove stray labels not in your metadata
extra_labels <- setdiff(df_tree_sorted$label, geo$sequence_id)
print(extra_labels)

df_tree_sorted <- df_tree_sorted %>%
  filter(!(label %in% extra_labels))                     # Drop unmatched labels

# Step 4: Pull metadata of interest and check alignment
geo_selected <- geo %>%
  select(sequence_id, Region, Country) %>%
  arrange(sequence_id) %>%
  mutate(label = df_tree_sorted$label)                   # Assign label from tree

# Step 5: Confirm that labels match
test_if_labels_match <- geo_selected %>%
  select(sequence_id, label) %>%
  mutate(test = ifelse(sequence_id == label, "match", "mismatch")) %>%
  arrange(desc(test))
print(test_if_labels_match)

# Step 6: Add metadata to tips
df_tree_tips <- df_tree_sorted %>%
  mutate(is_tip = TRUE) %>%
  merge(geo_selected, by = "label", all.x = TRUE)

# Step 7: Prepare internal nodes with bootstrap values
df_tree_bootstrap <- var2csa_df_tree %>%
  mutate(bootstrap = as.numeric(label)) %>%
  filter(!is.na(bootstrap)) %>%
  mutate(is_tip = FALSE)

# Step 8: Combine tips and internal nodes into full tree
df_tree_combined <- bind_rows(df_tree_tips, df_tree_bootstrap)


# Visualize the tree
library(ggtree)
library(ggplot2)
library(ggnewscale)

# Generate tree
phylo <- ggtree(var2csa_tree, layout = "rectangular")

# Look at node numbers
phylo + geom_text(aes(label = node), size = 2)  # Add node numbers

# Adjust root position (midrooting tree)
adj_phylo <- viewClade(phylo, node = 1366, xmax_adjust = 2)

# Visualize tree with bootstrap values and metadata (subtypes and continents)
adj_phylo + 
  geom_treescale(x = 0.2, y = 300, fontsize = 20, linesize = 2.5, offset = 12) +
  geom_tiplab(align = T, linesize = 0.05, size = 0.5) +  # Add sample names
  

  # Heatmap for Continent
  new_scale_fill() +
  geom_fruit(
    data = df_tree_combined,
    geom = geom_tile,
    mapping = aes(
      y = label, 
      x = "",
      fill = Region,
      width = 0.4
    ),
    offset = 1.0,
    pwidth = 2.0
  ) +
  scale_fill_manual(
    breaks = c("Africa", "Asia", "S.America"),
    labels = c("Africa", "Asia", "South America"),
    values = c("Africa" = "#8679b5", "Asia" = "#c781b1", "S.America" = "#baace6")
  ) +
  labs(fill = "Continent") +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15)
  )

bootstrap_summary <- df_tree_combined %>%
  filter(!is.na(bootstrap)) %>%   # Exclude NA values if any
  summarise(
    Mean = mean(bootstrap),         # Mean of bootstrap values
    Median = median(bootstrap),     # Median of bootstrap values
    Standard_Deviation = sd(bootstrap),  # Standard deviation of bootstrap values
    Min = min(bootstrap),           # Minimum value
    Max = max(bootstrap),           # Maximum value
    Count = n(),                    # Count of bootstrap values
    NA_Count = sum(is.na(bootstrap)) # Count of NAs (just in case)
  )

# Print the summary statistics
print(bootstrap_summary)



