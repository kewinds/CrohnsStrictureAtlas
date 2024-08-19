setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/newHuman"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages_updated",.libPaths()))

library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(tidyr)
library(Seurat)

#### Load controlled dataset ####
# Define the colors (assuming you have the correct color codes)
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", 
                     "#cfe2f3", "#6fa8dc", "#0b5394", "#073763")
                     
# Load the dataset
sub <- readRDS("fibrocontrolled_dataset.rds")

# Create new labels based on disease state
sub$newlabels <- "NA"
sub$newlabels[sub$diseasestate == "Normal MAT"] <- "Normal MAT"
sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
sub$newlabels[sub$diseasestate == "Involved"] <- "CF"
sub$newlabels[sub$diseasestate == "Normal Bowel"] <- "Normal Bowel"
sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"

#### Disease state ####

# Create a table and calculate proportions
tab <- table(sub$patientID, sub$newlabels, sub$fibro.annotation)
ptab <- prop.table(tab, 1) * 100

# Convert table to a data frame
ptab_df <- as.data.frame(ptab)
colnames(ptab_df) <- c("PatientID", "DiseaseState", "Fibroblast", "Frequency")

# Ensure the order of DiseaseState
desired_order <- c("Normal Bowel", "Uninvolved Bowel", "Stricture", 
                   "Normal MAT", "Uninvolved MAT", "CF")
ptab_df$DiseaseState <- factor(ptab_df$DiseaseState, levels = desired_order)

# Filter out rows where all frequencies for a PatientID within a DiseaseState are zero
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  filter(sum(Frequency) > 0) %>%
  ungroup()

# Recalculate frequencies so that they sum to 100% within each PatientID and DiseaseState
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  mutate(Frequency = Frequency / sum(Frequency) * 100) %>%
  ungroup()

# Calculate mean, standard deviation, and confidence intervals
ptab_stats <- ptab_df %>%
  group_by(DiseaseState, Fibroblast) %>%
  summarise(
    Mean = mean(Frequency),
    SD = sd(Frequency),
    n = n(),
    SE = SD / sqrt(n),
    CI_Lower = pmax(0, Mean - qt(1 - (0.05 / 2), df = n - 1) * SE),  # Ensure lower bound is at least 0
    CI_Upper = pmin(100, Mean + qt(1 - (0.05 / 2), df = n - 1) * SE)  # Cap upper bound at 100
  )
# Plot
plot <- ggplot(ptab_df, aes(x = Fibroblast, y = Frequency, color = Fibroblast)) + 
  geom_beeswarm(cex = 0.8, size = 1.5) +
  geom_errorbar(data = ptab_stats, aes(y = Mean, ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
  geom_point(data = ptab_stats, aes(y = Mean), size = 3, color = "black", shape = 18) +
  theme_classic() +
  labs(title = "Fibroblast Disease State Enrichment", 
       x = "Fibroblast Cluster", y = "Percent (%)", color = "Cluster") + 
  scale_color_manual(values = colors) + 
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.placement = "outside", 
        strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title.x = element_blank(),
        title = element_blank()) + 
  facet_grid(~DiseaseState, scales = "free_x", space = "free_x", switch = "x")
# Perform paired t-tests
t_test_results <- ptab_df %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"),
         DiseaseState %in% c("Stricture", "Uninvolved Bowel", "Normal Bowel", "CF", "Uninvolved MAT", "Normal MAT")) %>%
  group_by(Fibroblast) %>%
  summarise(
    p_value_CF_NormalMAT = t.test(Frequency[DiseaseState == "CF"], Frequency[DiseaseState == "Normal MAT"])$p.value,
    p_value_CF_UninvolvedMAT = t.test(Frequency[DiseaseState == "CF"], Frequency[DiseaseState == "Uninvolved MAT"])$p.value,
    p_value_Str_NormalBowel = t.test(Frequency[DiseaseState == "Stricture"], Frequency[DiseaseState == "Normal Bowel"])$p.value,
    p_value_Str_UninvolvedBowel = t.test(Frequency[DiseaseState == "Stricture"], Frequency[DiseaseState == "Uninvolved Bowel"])$p.value,
    p_value_UninvolvedBowel_NormalBowel = t.test(Frequency[DiseaseState == "Uninvolved Bowel"], Frequency[DiseaseState == "Normal Bowel"])$p.value,
    p_value_UninvolvedMAT_NormalMAT = t.test(Frequency[DiseaseState == "Uninvolved MAT"], Frequency[DiseaseState == "Normal MAT"])$p.value
  ) %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"))

annotations <- t_test_results %>%
  select(Fibroblast, p_value_CF_NormalMAT, p_value_CF_UninvolvedMAT, p_value_Str_NormalBowel, p_value_Str_UninvolvedBowel, p_value_UninvolvedBowel_NormalBowel,  p_value_UninvolvedMAT_NormalMAT) %>%
  pivot_longer(cols = starts_with("p_value_"), names_to = "Comparison", values_to = "p_value") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
annotations #view significance

#### Therapy ####
#Generate proptable and barplots of fibroblasts in drug treatments 
treat <- sub
Idents(treat) <- "diseasestate"
treat <- subset(treat, idents = c("Normal Bowel", "Stricture", "Normal MAT", "Involved"))
# Change labels of normal to "None"
treat$treated[treat$diseasestate == "Normal Bowel"] <- "Non-CD"
treat$treated[treat$diseasestate == "Normal MAT"] <- "Non-CD"
treat <- subset(treat, treated == "NA", invert = T) #Remove cells where treatment unknown
treat$treated <- factor(treat$treated, levels = c("Non-CD", "IM", "Biologic", "IM + Biologic"))
tab = table(treat$patientID, treat$treated, treat$fibro.annotation)
ptab = prop.table(tab, 1)*100

# Convert table to a data frame
ptab_df <- as.data.frame(ptab)
head(ptab_df)
colnames(ptab_df) <- c("PatientID", "Treatment", "Fibroblast", "Frequency")
head(ptab_df)
# Ensure the order of therapies
desired_order <- c("Non-CD", "IM", "Biologic", 
                   "IM + Biologic")
ptab_df$Treatment <- factor(ptab_df$Treatment, levels = desired_order)

# Filter out rows where all frequencies for a PatientID within a Treatment are zero
ptab_df <- ptab_df %>%
  group_by(PatientID, Treatment) %>%
  filter(sum(Frequency) > 0) %>%
  ungroup()

# Recalculate frequencies so that they sum to 100% within each PatientID and DiseaseState
ptab_df <- ptab_df %>%
  group_by(PatientID, Treatment) %>%
  mutate(Frequency = Frequency / sum(Frequency) * 100) %>%
  ungroup()

# Calculate mean, standard deviation, and confidence intervals
ptab_stats <- ptab_df %>%
  group_by(Treatment, Fibroblast) %>%
  summarise(
    Mean = mean(Frequency),
    SD = sd(Frequency),
    n = n(),
    SE = SD / sqrt(n),
    CI_Lower = pmax(0, Mean - qt(1 - (0.05 / 2), df = n - 1) * SE),  # Ensure lower bound is at least 0
    CI_Upper = pmin(100, Mean + qt(1 - (0.05 / 2), df = n - 1) * SE)  # Cap upper bound at 100
  )
# Perform t-tests and prepare p-values only for specific fibroblast populations
t_test_results <- ptab_df %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"),
         Treatment %in% c("Non-CD", "IM", "Biologic", "IM + Biologic")) %>%
  group_by(Fibroblast) %>%
  summarise(
    p_value_Bio_IMBio = t.test(Frequency[Treatment == "Biologic"], Frequency[Treatment == "IM + Biologic"])$p.value,
    p_value_Bio_NonCD = t.test(Frequency[Treatment == "Biologic"], Frequency[Treatment == "Non-CD"])$p.value,
    p_value_IMBio_NonCD = t.test(Frequency[Treatment == "IM + Biologic"], Frequency[Treatment == "Non-CD"])$p.value,
   ) %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"))

annotations <- t_test_results %>%
  select(Fibroblast, p_value_Bio_IMBio, p_value_Bio_NonCD, p_value_IMBio_NonCD) %>%
  pivot_longer(cols = starts_with("p_value_"), names_to = "Comparison", values_to = "p_value") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
annotations #view significance
# Plot
ggplot(ptab_df, aes(x = Fibroblast, y = Frequency, color = Fibroblast)) + 
  geom_beeswarm(cex = 0.8, size = 1.5) +
  geom_errorbar(data = ptab_stats, aes(y = Mean, ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
  geom_point(data = ptab_stats, aes(y = Mean), size = 3, color = "black", shape = 18) +
  theme_classic() +
  labs(title = "Fibroblast Disease State Enrichment", 
       x = "Fibroblast Cluster", y = "Percent (%)", color = "Cluster") + 
  scale_color_manual(values = colors) + 
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.placement = "outside", 
        strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title.x = element_blank(),
        title = element_blank()) + 
  facet_grid(~Treatment, scales = "free_x", space = "free_x", switch = "x")
ggsave("SuppFigures/beeswarmtreatment.svg")


#### Age (Adult) ####
adult <- sub
Idents(adult) <- "age"
adult <- subset(adult, idents = c("Adult"))
tab = table(adult$patientID, adult$treated, adult$diseasestate, adult$fibro.annotation)
ptab = prop.table(tab, 1)*100

# Create new labels based on disease state
adult$newlabels <- "NA"
adult$newlabels[adult$diseasestate == "Normal MAT"] <- "Normal MAT"
adult$newlabels[adult$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
adult$newlabels[adult$diseasestate == "Involved"] <- "CF"
adult$newlabels[adult$diseasestate == "Normal Bowel"] <- "Normal Bowel"
adult$newlabels[adult$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
adult$newlabels[adult$diseasestate == "Stricture"] <- "Stricture"

# Create a table and calculate proportions
tab <- table(adult$patientID, adult$newlabels, adult$fibro.annotation)
ptab <- prop.table(tab, 1) * 100

# Convert table to a data frame
ptab_df <- as.data.frame(ptab)
colnames(ptab_df) <- c("PatientID", "DiseaseState", "Fibroblast", "Frequency")
head(ptab_df)
# Ensure the order of DiseaseState
desired_order <- c("Normal Bowel", "Uninvolved Bowel", "Stricture", 
                   "Normal MAT", "Uninvolved MAT", "CF")
ptab_df$DiseaseState <- factor(ptab_df$DiseaseState, levels = desired_order)

# Filter out rows where all frequencies for a PatientID within a DiseaseState are zero
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  filter(sum(Frequency) > 0) %>%
  ungroup()

# Recalculate frequencies so that they sum to 100% within each PatientID and DiseaseState
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  mutate(Frequency = Frequency / sum(Frequency) * 100) %>%
  ungroup()

# Calculate mean, standard deviation, and confidence intervals
ptab_stats <- ptab_df %>%
  group_by(DiseaseState, Fibroblast) %>%
  summarise(
    Mean = mean(Frequency),
    SD = sd(Frequency),
    n = n(),
    SE = SD / sqrt(n),
    CI_Lower = pmax(0, Mean - qt(1 - (0.05 / 2), df = n - 1) * SE),  # Ensure lower bound is at least 0
    CI_Upper = pmin(100, Mean + qt(1 - (0.05 / 2), df = n - 1) * SE)  # Cap upper bound at 100
  )
# Perform t-tests and prepare p-values only for specific fibroblast populations
t_test_results <- ptab_df %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"),
         DiseaseState %in% c("Stricture", "Uninvolved Bowel", "Normal Bowel", "CF", "Uninvolved MAT", "Normal MAT")) %>%
  group_by(Fibroblast) %>%
  summarise(
    p_value_CF_NormalMAT = t.test(Frequency[DiseaseState == "CF"], Frequency[DiseaseState == "Normal MAT"])$p.value,
    p_value_CF_UninvolvedMAT = t.test(Frequency[DiseaseState == "CF"], Frequency[DiseaseState == "Uninvolved MAT"])$p.value,
    p_value_Str_NormalBowel = t.test(Frequency[DiseaseState == "Stricture"], Frequency[DiseaseState == "Normal Bowel"])$p.value,
    p_value_Str_UninvolvedBowel = t.test(Frequency[DiseaseState == "Stricture"], Frequency[DiseaseState == "Uninvolved Bowel"])$p.value,
    p_value_UninvolvedBowel_NormalBowel = t.test(Frequency[DiseaseState == "Uninvolved Bowel"], Frequency[DiseaseState == "Normal Bowel"])$p.value,
    p_value_UninvolvedMAT_NormalMAT = t.test(Frequency[DiseaseState == "Uninvolved MAT"], Frequency[DiseaseState == "Normal MAT"])$p.value
  ) %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"))

annotations <- t_test_results %>%
  select(Fibroblast, p_value_CF_NormalMAT, p_value_CF_UninvolvedMAT, p_value_Str_NormalBowel, p_value_Str_UninvolvedBowel, p_value_UninvolvedBowel_NormalBowel,  p_value_UninvolvedMAT_NormalMAT) %>%
  pivot_longer(cols = starts_with("p_value_"), names_to = "Comparison", values_to = "p_value") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
annotations #view significance
# Plot
ggplot(ptab_df, aes(x = Fibroblast, y = Frequency, color = Fibroblast)) + 
  geom_beeswarm(cex = 0.8, size = 1.5) +
  geom_errorbar(data = ptab_stats, aes(y = Mean, ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
  geom_point(data = ptab_stats, aes(y = Mean), size = 3, color = "black", shape = 18) +
  theme_classic() +
  labs(title = "Fibroblast Disease State Enrichment", 
       x = "Fibroblast Cluster", y = "Percent (%)", color = "Cluster") + 
  scale_color_manual(values = colors) + 
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.placement = "outside", 
        strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title.x = element_blank(),
        title = element_blank()) + 
  facet_grid(~DiseaseState, scales = "free_x", space = "free_x", switch = "x")
ggsave("SuppFigures/beeswarmadult.svg")
#### Age (Pediatric) ####
adult <- sub
Idents(adult) <- "age"
adult <- subset(adult, idents = c("Pediatric"))
tab = table(adult$patientID, adult$treated, adult$diseasestate, adult$fibro.annotation)
ptab = prop.table(tab, 1)*100

# Create new labels based on disease state
adult$newlabels <- "NA"
adult$newlabels[adult$diseasestate == "Normal MAT"] <- "Normal MAT"
adult$newlabels[adult$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
adult$newlabels[adult$diseasestate == "Involved"] <- "CF"
adult$newlabels[adult$diseasestate == "Normal Bowel"] <- "Normal Bowel"
adult$newlabels[adult$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
adult$newlabels[adult$diseasestate == "Stricture"] <- "Stricture"

# Create a table and calculate proportions
tab <- table(adult$patientID, adult$newlabels, adult$fibro.annotation)
ptab <- prop.table(tab, 1) * 100

# Convert table to a data frame
ptab_df <- as.data.frame(ptab)
colnames(ptab_df) <- c("PatientID", "DiseaseState", "Fibroblast", "Frequency")
head(ptab_df)
# Ensure the order of DiseaseState
desired_order <- c("Normal Bowel", "Uninvolved Bowel", "Stricture", 
                   "Normal MAT", "Uninvolved MAT", "CF")
ptab_df$DiseaseState <- factor(ptab_df$DiseaseState, levels = desired_order)

# Filter out rows where all frequencies for a PatientID within a DiseaseState are zero
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  filter(sum(Frequency) > 0) %>%
  ungroup()

# Recalculate frequencies so that they sum to 100% within each PatientID and DiseaseState
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  mutate(Frequency = Frequency / sum(Frequency) * 100) %>%
  ungroup()

# Calculate mean, standard deviation, and confidence intervals
ptab_stats <- ptab_df %>%
  group_by(DiseaseState, Fibroblast) %>%
  summarise(
    Mean = mean(Frequency),
    SD = sd(Frequency),
    n = n(),
    SE = SD / sqrt(n),
    CI_Lower = pmax(0, Mean - qt(1 - (0.05 / 2), df = n - 1) * SE),  # Ensure lower bound is at least 0
    CI_Upper = pmin(100, Mean + qt(1 - (0.05 / 2), df = n - 1) * SE)  # Cap upper bound at 100
  )
# Perform t-tests and prepare p-values only for specific fibroblast populations
t_test_results <- ptab_df %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"),
         DiseaseState %in% c("Stricture", "Uninvolved Bowel", "CF", "Uninvolved MAT")) %>%
  group_by(Fibroblast) %>%
  summarise(
    p_value_CF_UninvolvedMAT = t.test(Frequency[DiseaseState == "CF"], Frequency[DiseaseState == "Uninvolved MAT"])$p.value,
    p_value_Str_UninvolvedBowel = t.test(Frequency[DiseaseState == "Stricture"], Frequency[DiseaseState == "Uninvolved Bowel"])$p.value,
  ) %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"))

annotations <- t_test_results %>%
  select(Fibroblast, p_value_CF_UninvolvedMAT, p_value_Str_UninvolvedBowel) %>%
  pivot_longer(cols = starts_with("p_value_"), names_to = "Comparison", values_to = "p_value") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
annotations #view significance

# Plot
ggplot(ptab_df, aes(x = Fibroblast, y = Frequency, color = Fibroblast)) + 
  geom_beeswarm(cex = 0.8, size = 1.5) +
  geom_errorbar(data = ptab_stats, aes(y = Mean, ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
  geom_point(data = ptab_stats, aes(y = Mean), size = 3, color = "black", shape = 18) +
  theme_classic() +
  labs(title = "Fibroblast Disease State Enrichment", 
       x = "Fibroblast Cluster", y = "Percent (%)", color = "Cluster") + 
  scale_color_manual(values = colors) + 
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.placement = "outside", 
        strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title.x = element_blank(),
        title = element_blank()) + 
  facet_grid(~DiseaseState, scales = "free_x", space = "free_x", switch = "x")
ggsave("SuppFigures/beeswarmpeds.svg")


#### Sex (Male) ####
adult <- sub
Idents(adult) <- "sex"
adult <- subset(adult, idents = c("Male"))
tab = table(adult$patientID, adult$treated, adult$diseasestate, adult$fibro.annotation)
ptab = prop.table(tab, 1)*100

# Create new labels based on disease state
adult$newlabels <- "NA"
adult$newlabels[adult$diseasestate == "Normal MAT"] <- "Normal MAT"
adult$newlabels[adult$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
adult$newlabels[adult$diseasestate == "Involved"] <- "CF"
adult$newlabels[adult$diseasestate == "Normal Bowel"] <- "Normal Bowel"
adult$newlabels[adult$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
adult$newlabels[adult$diseasestate == "Stricture"] <- "Stricture"

# Create a table and calculate proportions
tab <- table(adult$patientID, adult$newlabels, adult$fibro.annotation)
ptab <- prop.table(tab, 1) * 100

# Convert table to a data frame
ptab_df <- as.data.frame(ptab)
colnames(ptab_df) <- c("PatientID", "DiseaseState", "Fibroblast", "Frequency")
head(ptab_df)
# Ensure the order of DiseaseState
desired_order <- c("Normal Bowel", "Uninvolved Bowel", "Stricture", 
                   "Normal MAT", "Uninvolved MAT", "CF")
ptab_df$DiseaseState <- factor(ptab_df$DiseaseState, levels = desired_order)

# Filter out rows where all frequencies for a PatientID within a DiseaseState are zero
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  filter(sum(Frequency) > 0) %>%
  ungroup()

# Recalculate frequencies so that they sum to 100% within each PatientID and DiseaseState
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  mutate(Frequency = Frequency / sum(Frequency) * 100) %>%
  ungroup()

# Calculate mean, standard deviation, and confidence intervals
ptab_stats <- ptab_df %>%
  group_by(DiseaseState, Fibroblast) %>%
  summarise(
    Mean = mean(Frequency),
    SD = sd(Frequency),
    n = n(),
    SE = SD / sqrt(n),
    CI_Lower = pmax(0, Mean - qt(1 - (0.05 / 2), df = n - 1) * SE),  # Ensure lower bound is at least 0
    CI_Upper = pmin(100, Mean + qt(1 - (0.05 / 2), df = n - 1) * SE)  # Cap upper bound at 100
  )
# Perform t-tests and prepare p-values only for specific fibroblast populations
t_test_results <- ptab_df %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"),
         DiseaseState %in% c("Stricture", "Uninvolved Bowel", "Normal Bowel", "CF", "Uninvolved MAT")) %>%
  group_by(Fibroblast) %>%
  summarise(
    p_value_CF_UninvolvedMAT = t.test(Frequency[DiseaseState == "CF"], Frequency[DiseaseState == "Uninvolved MAT"])$p.value,
    p_value_Str_NormalBowel = t.test(Frequency[DiseaseState == "Stricture"], Frequency[DiseaseState == "Normal Bowel"])$p.value,
    p_value_Str_UninvolvedBowel = t.test(Frequency[DiseaseState == "Stricture"], Frequency[DiseaseState == "Uninvolved Bowel"])$p.value,
    p_value_UninvolvedBowel_NormalBowel = t.test(Frequency[DiseaseState == "Uninvolved Bowel"], Frequency[DiseaseState == "Normal Bowel"])$p.value,
  ) %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"))

annotations <- t_test_results %>%
  select(Fibroblast, p_value_CF_UninvolvedMAT, p_value_Str_NormalBowel, p_value_Str_UninvolvedBowel, p_value_UninvolvedBowel_NormalBowel) %>%
  pivot_longer(cols = starts_with("p_value_"), names_to = "Comparison", values_to = "p_value") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
annotations #view significance

# Plot
ggplot(ptab_df, aes(x = Fibroblast, y = Frequency, color = Fibroblast)) + 
  geom_beeswarm(cex = 0.8, size = 1.5) +
  geom_errorbar(data = ptab_stats, aes(y = Mean, ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
  geom_point(data = ptab_stats, aes(y = Mean), size = 3, color = "black", shape = 18) +
  theme_classic() +
  labs(title = "Fibroblast Disease State Enrichment", 
       x = "Fibroblast Cluster", y = "Percent (%)", color = "Cluster") + 
  scale_color_manual(values = colors) + 
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.placement = "outside", 
        strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title.x = element_blank(),
        title = element_blank()) + 
  facet_grid(~DiseaseState, scales = "free_x", space = "free_x", switch = "x")

ggsave("SuppFigures/beeswarmmale.svg")

#### Sex (Female) ####
adult <- sub
Idents(adult) <- "sex"
adult <- subset(adult, idents = c("Female"))
tab = table(adult$patientID, adult$treated, adult$diseasestate, adult$fibro.annotation)
ptab = prop.table(tab, 1)*100

# Create new labels based on disease state
adult$newlabels <- "NA"
adult$newlabels[adult$diseasestate == "Normal MAT"] <- "Normal MAT"
adult$newlabels[adult$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
adult$newlabels[adult$diseasestate == "Involved"] <- "CF"
adult$newlabels[adult$diseasestate == "Normal Bowel"] <- "Normal Bowel"
adult$newlabels[adult$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
adult$newlabels[adult$diseasestate == "Stricture"] <- "Stricture"

# Create a table and calculate proportions
tab <- table(adult$patientID, adult$newlabels, adult$fibro.annotation)
ptab <- prop.table(tab, 1) * 100

# Convert table to a data frame
ptab_df <- as.data.frame(ptab)
colnames(ptab_df) <- c("PatientID", "DiseaseState", "Fibroblast", "Frequency")
head(ptab_df)
# Ensure the order of DiseaseState
desired_order <- c("Normal Bowel", "Uninvolved Bowel", "Stricture", 
                   "Normal MAT", "Uninvolved MAT", "CF")
ptab_df$DiseaseState <- factor(ptab_df$DiseaseState, levels = desired_order)

# Filter out rows where all frequencies for a PatientID within a DiseaseState are zero
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  filter(sum(Frequency) > 0) %>%
  ungroup()

# Recalculate frequencies so that they sum to 100% within each PatientID and DiseaseState
ptab_df <- ptab_df %>%
  group_by(PatientID, DiseaseState) %>%
  mutate(Frequency = Frequency / sum(Frequency) * 100) %>%
  ungroup()

# Calculate mean, standard deviation, and confidence intervals
ptab_stats <- ptab_df %>%
  group_by(DiseaseState, Fibroblast) %>%
  summarise(
    Mean = mean(Frequency),
    SD = sd(Frequency),
    n = n(),
    SE = SD / sqrt(n),
    CI_Lower = pmax(0, Mean - qt(1 - (0.05 / 2), df = n - 1) * SE),  # Ensure lower bound is at least 0
    CI_Upper = pmin(100, Mean + qt(1 - (0.05 / 2), df = n - 1) * SE)  # Cap upper bound at 100
  )
# Perform t-tests and prepare p-values only for specific fibroblast populations
t_test_results <- ptab_df %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"),
         DiseaseState %in% c("Stricture", "Uninvolved Bowel", "Normal Bowel", "CF", "Uninvolved MAT", "Normal MAT")) %>%
  group_by(Fibroblast) %>%
  summarise(
    p_value_CF_NormalMAT = t.test(Frequency[DiseaseState == "CF"], Frequency[DiseaseState == "Normal MAT"])$p.value,
    p_value_CF_UninvolvedMAT = t.test(Frequency[DiseaseState == "CF"], Frequency[DiseaseState == "Uninvolved MAT"])$p.value,
    p_value_Str_UninvolvedBowel = t.test(Frequency[DiseaseState == "Stricture"], Frequency[DiseaseState == "Uninvolved Bowel"])$p.value,
    p_value_UninvolvedMAT_NormalMAT = t.test(Frequency[DiseaseState == "Uninvolved MAT"], Frequency[DiseaseState == "Normal MAT"])$p.value
  ) %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"))

annotations <- t_test_results %>%
  select(Fibroblast, p_value_CF_NormalMAT, p_value_CF_UninvolvedMAT, p_value_Str_UninvolvedBowel,  p_value_UninvolvedMAT_NormalMAT) %>%
  pivot_longer(cols = starts_with("p_value_"), names_to = "Comparison", values_to = "p_value") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
annotations #view significance


# Plot
ggplot(ptab_df, aes(x = Fibroblast, y = Frequency, color = Fibroblast)) + 
  geom_beeswarm(cex = 0.8, size = 1.5) +
  geom_errorbar(data = ptab_stats, aes(y = Mean, ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
  geom_point(data = ptab_stats, aes(y = Mean), size = 3, color = "black", shape = 18) +
  theme_classic() +
  labs(title = "Fibroblast Disease State Enrichment", 
       x = "Fibroblast Cluster", y = "Percent (%)", color = "Cluster") + 
  scale_color_manual(values = colors) + 
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.placement = "outside", 
        strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title.x = element_blank(),
        title = element_blank()) + 
  facet_grid(~DiseaseState, scales = "free_x", space = "free_x", switch = "x")

ggsave("SuppFigures/beeswarmfemale.svg")
#### Therapy ####
#Generate proptable and barplots of fibroblasts in drug treatments 
treat <- sub
Idents(treat) <- "diseasestate"
treat <- subset(treat, idents = c("Normal Bowel", "Stricture", "Normal MAT", "Involved"))
# Change labels of normal to "None"
treat$treated[treat$diseasestate == "Normal Bowel"] <- "Non-CD"
treat$treated[treat$diseasestate == "Normal MAT"] <- "Non-CD"
treat <- subset(treat, treated == "NA", invert = T) #Remove cells where treatment unknown
treat$treated <- factor(treat$treated, levels = c("Non-CD", "IM", "Biologic", "IM + Biologic"))
tab = table(treat$patientID, treat$treated, treat$fibro.annotation)
ptab = prop.table(tab, 1)*100

# Convert table to a data frame
ptab_df <- as.data.frame(ptab)
head(ptab_df)
colnames(ptab_df) <- c("PatientID", "Treatment", "Fibroblast", "Frequency")
head(ptab_df)
# Ensure the order of therapies
desired_order <- c("Non-CD", "IM", "Biologic", 
                   "IM + Biologic")
ptab_df$Treatment <- factor(ptab_df$Treatment, levels = desired_order)

# Filter out rows where all frequencies for a PatientID within a Treatment are zero
ptab_df <- ptab_df %>%
  group_by(PatientID, Treatment) %>%
  filter(sum(Frequency) > 0) %>%
  ungroup()

# Recalculate frequencies so that they sum to 100% within each PatientID and DiseaseState
ptab_df <- ptab_df %>%
  group_by(PatientID, Treatment) %>%
  mutate(Frequency = Frequency / sum(Frequency) * 100) %>%
  ungroup()

# Calculate mean, standard deviation, and confidence intervals
ptab_stats <- ptab_df %>%
  group_by(Treatment, Fibroblast) %>%
  summarise(
    Mean = mean(Frequency),
    SD = sd(Frequency),
    n = n(),
    SE = SD / sqrt(n),
    CI_Lower = pmax(0, Mean - qt(1 - (0.05 / 2), df = n - 1) * SE),  # Ensure lower bound is at least 0
    CI_Upper = pmin(100, Mean + qt(1 - (0.05 / 2), df = n - 1) * SE)  # Cap upper bound at 100
  )
# Perform t-tests and prepare p-values only for specific fibroblast populations
t_test_results <- ptab_df %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"),
         Treatment %in% c("Non-CD", "IM", "Biologic", "IM + Biologic")) %>%
  group_by(Fibroblast) %>%
  summarise(
    p_value_Bio_IMBio = t.test(Frequency[Treatment == "Biologic"], Frequency[Treatment == "IM + Biologic"])$p.value,
    p_value_Bio_NonCD = t.test(Frequency[Treatment == "Biologic"], Frequency[Treatment == "Non-CD"])$p.value,
    p_value_IMBio_NonCD = t.test(Frequency[Treatment == "IM + Biologic"], Frequency[Treatment == "Non-CD"])$p.value,
  ) %>%
  filter(Fibroblast %in% c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+"))

annotations <- t_test_results %>%
  select(Fibroblast, p_value_Bio_IMBio, p_value_Bio_NonCD, p_value_IMBio_NonCD) %>%
  pivot_longer(cols = starts_with("p_value_"), names_to = "Comparison", values_to = "p_value") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
annotations #view significance
# Plot
ggplot(ptab_df, aes(x = Fibroblast, y = Frequency, color = Fibroblast)) + 
  geom_beeswarm(cex = 0.8, size = 1.5) +
  geom_errorbar(data = ptab_stats, aes(y = Mean, ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
  geom_point(data = ptab_stats, aes(y = Mean), size = 3, color = "black", shape = 18) +
  theme_classic() +
  labs(title = "Fibroblast Disease State Enrichment", 
       x = "Fibroblast Cluster", y = "Percent (%)", color = "Cluster") + 
  scale_color_manual(values = colors) + 
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.placement = "outside", 
        strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title.x = element_blank(),
        title = element_blank()) + 
  facet_grid(~Treatment, scales = "free_x", space = "free_x", switch = "x")




