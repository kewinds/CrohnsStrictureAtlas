library(Seurat)
library(dplyr)
library(ggplot2)

# Assuming 'fibro' is your Seurat object

# Step 1: Create TRUE unique patient IDs by combining patientID with sampleID (or another unique identifier)
# (Assuming you have a column like 'sampleID' to distinguish multiple samples from the same patient)
fibro@meta.data$unique_patientID <- paste(fibro@meta.data$patientID, fibro@meta.data$orig.ident, sep = "_")

# Step 2: Properly exclude NA values in 'treated' column
fibro <- subset(fibro, subset == 'NA')

fibro$general <- NA
fibro$general[fibro$diseasestate == "Stricture"] <- "Stricture"
fibro$general[fibro$diseasestate == "Stricture Uninflamed"] <- "Uninvolved"
fibro$general[fibro$diseasestate == "Uninvolved"] <- "Uninvolved"
fibro$general[fibro$diseasestate == "Involved"] <- "Stricture"


# Step 3: Calculate proportions (corrected counting method)
data_summary <- fibro@meta.data %>%
  group_by(unique_patientID, treated, general, tissue) %>%
  summarize(
    CTHRC1_mFib_count = sum(fibro_annotation == 'CTHRC1+ mFib'),  # Count CTHRC1+ cells
    total_fib_count = n()  # Count total fibroblasts
  ) %>%
  mutate(proportion = CTHRC1_mFib_count / total_fib_count) %>%
  ungroup()

# Step 4: Proper ANCOVA model specification
# Convert to factors (ensure no NA levels remain)
data_summary <- data_summary %>%
  mutate(across(c(treated, general, tissue), as.factor))

# Run ANCOVA with interaction term and covariates
ancova_results <- aov(proportion ~ treated + general, data = data_summary)

# Reorder 'general' factor levels
data_summary$general <- factor(data_summary$general, 
                               levels = c("Uninvolved", "Stricture"))

# Calculate mean and standard deviation
summary_stats <- data_summary %>%
  group_by(general) %>%
  summarise(mean_prop = mean(proportion),
            sd_prop = sd(proportion))

# Fit ANCOVA model with treatment as covariate
ancova <- lm(proportion ~ treated + general, data = data_summary)

# Get p-value for 'treated' from ancova (with proper formatting)
treated_p_value <- format.pval(summary(ancova)$coefficients["treated", "Pr(>|t|)"], digits = 3)
cat("P-value for treatment effect:", treated_p_value, "\n")

# Get adjusted means (LS means) for disease state
library(emmeans)
adjusted_means <- emmeans(ancova, ~ general) %>% 
  as.data.frame() %>% 
  mutate(model = factor("Adjusted", 
                        levels = c("Unadjusted", "Adjusted")))

# Get unadjusted means
unadjusted_means <- data_summary %>% 
  group_by(general) %>% 
  summarise(emmean = mean(proportion),
            SE = sd(proportion)/sqrt(n())) %>% 
  mutate(model = factor("Unadjusted", 
                        levels = c("Unadjusted", "Adjusted")))

# Combine results
combined_means <- bind_rows(unadjusted_means, adjusted_means)

# Create comparison visualization
library(ggplot2)
ggplot(combined_means, aes(x = general, y = emmean, color = model)) +
  geom_pointrange(aes(ymin = emmean - SE, ymax = emmean + SE),
                  position = position_dodge(width = 0.6),
                  size = 1,
                  fatten = 4) +  # Increased from 1.5 to 4
  labs(
       y = "CTHRC1+ mFib Enrichment",
       color = "Analysis Type") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 8),
        legend.title = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"), 
        aspect.ratio = 1) +
  scale_color_manual(values = c("#1f77b4", "#8b0000")) +
  coord_cartesian(ylim = c(0, max(combined_means$emmean) * 1.2)) 
  
ggsave("Figures/ancovaresults.pdf")