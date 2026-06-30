
install.packages("gridExtra")
install.packages("vegan")
install.packages("rstatix")

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(permute)
library(vegan)
library(rstatix)
library(car)

# Read data
df <- read.csv("C:/Users/Herman/Desktop/Data_Burkmeer_Macroinvertebrates.csv", stringsAsFactors=TRUE)
df_ffg <- read.csv("C:/Users/Herman/Desktop/Data_FFG.csv", stringsAsFactors=TRUE)

{
# Summarise by Location, Method, and Taxon
df_sum <- df %>%
  group_by(Location, Method, Taxon) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(Abundance > 0)

# Unique locations and methods (Net first, Sampler second)
locations <- unique(df_sum$Location)
methods   <- c("Net", "Sampler")

# Pie chart function for a given location + method
make_pie <- function(data, loc, meth) {
  plot_data   <- data %>% filter(Location == loc, Method == meth)
  total_abund <- sum(plot_data$Abundance)
  
  ggplot(
    plot_data,
    aes(x = "", y = Abundance, fill = Taxon)
  ) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    labs(title = paste0(loc, " - ", meth, "\n(Total: ", total_abund, ")"),
         fill = "Taxon") +
    theme_void() +
    theme(
      plot.title       = element_text(size = 10, hjust = 0.5),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm"),
      plot.margin      = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

# Iterate methods first, then locations → top row = Net, bottom row = Sampler
plots <- lapply(methods, function(meth) {
  lapply(locations, function(loc) make_pie(df_sum, loc, meth))
}) |> unlist(recursive = FALSE)

# Export with fixed size
png("pie_charts.png", width = 12, height = 8, units = "in", res = 300)
grid.arrange(
  grobs   = plots,
  ncol    = 3,
  widths  = c(1, 1, 1),
  heights = c(1, 1)
)
dev.off()

}


### Total Abundance ANOVA + tukey post-hoc 
# Check assumptions of normality and equal variances
{
#Summarise total abundance per sample group
  df_abundance <- df %>%
    group_by(Method, Location, Sample) %>%
    summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
#Normality per group (Shapiro-Wilk)
  df_abundance %>%
    group_by(Method, Location) %>%
    summarise(
      p_normality = tryCatch(
        shapiro.test(Total_Abundance)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    )
  
#Homogeneity of variances (Levene's test)

  leveneTest(Total_Abundance ~ interaction(Method, Location), data = df_abundance)
  
# Visual check: QQ-plots per group
  ggplot(df_abundance, aes(sample = Total_Abundance)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~ Method + Location) +
    theme_classic() +
    labs(title = "QQ-plots per Method × Location")
}

{
# Summarise total abundance per sample
df_abund <- df %>%
  group_by(Method, Location, Sample) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

# Function: ANOVA + Tukey + compact letter display per method
y_max <- 22000  # shared y-axis limit

get_letters <- function(method_data, method_name) {
  cat("\n========================================\n")
  cat("Method:", method_name, "\n")
  cat("========================================\n")
  
  aov_result <- aov(Total_Abundance ~ Location, data = method_data)
  
  cat("\n--- ANOVA Summary ---\n")
  print(summary(aov_result))
  
  tukey <- TukeyHSD(aov_result)
  cat("\n--- Tukey HSD Post-hoc ---\n")
  print(tukey)
  
  tukey_loc <- tukey$Location
  locs      <- unique(method_data$Location)
  
  # Build a p-value matrix
  p_matrix <- matrix(1, nrow = length(locs), ncol = length(locs),
                     dimnames = list(locs, locs))
  
  for (i in seq_len(nrow(tukey_loc))) {
    pair  <- strsplit(rownames(tukey_loc)[i], "-")[[1]]
    p_val <- tukey_loc[i, "p adj"]
    p_matrix[pair[1], pair[2]] <- p_val
    p_matrix[pair[2], pair[1]] <- p_val
  }
  
  # Assign letters based on significance (p < 0.05)
  sig_diff       <- p_matrix < 0.05
  letters_vec    <- setNames(rep(NA, length(locs)), locs)
  current_letter <- 1
  
  for (loc in locs) {
    if (is.na(letters_vec[loc])) {
      letters_vec[loc] <- letters[current_letter]
      current_letter   <- current_letter + 1
    }
    for (other in locs) {
      if (loc != other && !sig_diff[loc, other] && is.na(letters_vec[other])) {
        letters_vec[other] <- letters_vec[loc]
      }
    }
  }
  
  # y position just above the max, but capped below the plot limit
  label_df <- method_data %>%
    group_by(Location) %>%
    summarise(y_pos = min(max(Total_Abundance) * 1.4, y_max * 0.92), .groups = "drop") %>%
    mutate(label = toupper(letters_vec[Location]))
  
  return(label_df)
}
# Apply per method, printing ANOVA output

letter_df <- df_abund %>%
  group_by(Method) %>%
  group_modify(~ get_letters(.x, .y$Method)) %>%
  ungroup() %>%
  mutate(Method = case_when(
    Method == "Net"     ~ "A) Net",
    Method == "Sampler" ~ "B) Sampler"
  ))

# Recode in df_abund for facet labels
df_abund <- df_abund %>%
  mutate(Method_label = dplyr::recode(Method, "Net" = "A) Net", "Sampler" = "B) Sampler"))

# Plot
ggplot(df_abund, aes(x = Location, y = Total_Abundance, fill = Method)) +
  geom_boxplot(
    outlier.shape  = 21,
    outlier.size   = 1.5,
    outlier.fill   = "white",
    width          = 0.6,
    color          = "black",
    linewidth      = 0.4
  ) +
  geom_text(
    data        = letter_df,
    aes(x = Location, y = y_pos, label = label),
    inherit.aes = FALSE,
    size        = 4,
    position    = position_nudge(x = 0.3)
  ) +
  facet_wrap(~ Method_label) +
  scale_y_log10(
    limits = c(1, y_max),
    breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 22000),
    labels = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 22000)
  ) +
  scale_fill_manual(values = c("Net" = "#F4A5A5", "Sampler" = "#85C7E0")) +
  labs(
    x    = "Location",
    y    = "Log(Total Abundance)",
    fill = "Method"
  ) +
  theme_classic() +
  theme(
    strip.text        = element_text(size = 11),
    strip.background  = element_blank(),
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9),
    axis.line         = element_line(color = "black", linewidth = 0.4),
    panel.border      = element_blank(),
    legend.position   = "none"
  )
}

# Calculate means per Location/Method
{df %>%
  group_by(LocMet, Sample) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE),
            .groups = "drop"
            ) %>%
  group_by(LocMet) %>%
  summarise(Mean_Abundance = round(mean(Total_Abundance),4))

# Perform pairwise T-tests
pairwise_results <- pairwise.t.test(df$Abundance, df$LocMet, p.adjust.method = "bonferroni")

# Print results
print(pairwise_results)
}

### Relative Abundance
{
  # Summarise and calculate relative abundance
  
  df_rel <- df %>%
    group_by(Method, Location, Taxon) %>%
    summarise(Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
              .groups = "drop")
  
  df_totals <- df %>%
    group_by(Method, Location) %>%
    summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  ggplot(df_rel, aes(x = Location, y = Relative_Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    facet_wrap(
      ~ Method,
      labeller = labeller(Method = c("Net" = "A) Net",
                                     "Sampler" = "B) Sampler"))
    ) +
    geom_text(data = df_totals,
              aes(x = Location,
              y = 100,
              label = paste0("n = ",
                             Total_Abundance
                             ),
              fill = NULL
              ),
              size = 3,
              hjust = 0.5,
              vjust = -0.5
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Location",
      y = "Relative abundance (%)",
      fill = "Taxon"
    ) +
    theme_classic() +
    theme(
      strip.text       = element_text(size = 15),
      strip.background = element_blank(),
      axis.title       = element_text(size = 15),
      axis.text        = element_text(size = 15),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      legend.text      = element_text(size = 15),
      legend.title     = element_text(size = 15),
      legend.key.size  = unit(0.4, "cm")
    )
}

### Functional Feeding Group distribution
{
  # Summarise and calculate relative abundance epr FFG
  df_rel <- df %>%
    group_by(Method, Location, FFG) %>%
    summarise(Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
              .groups = "drop")
  

  df_totals <- df %>%
    group_by(Method, Location) %>%
    summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  ggplot(df_rel,
         aes(x = Location,
             y = Relative_Abundance,
             fill = FFG
         )
  ) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    facet_wrap(
      ~ Method,
      labeller = labeller(Method = c("Net" = "A) Net",
                                     "Sampler" = "B) Sampler"
      )
      )
    ) +
    geom_text(data = df_totals,
              aes(x = Location,
                  y = 100,
                  label = paste0("n = ",
                                 Total_Abundance),
                  fill = NULL),
              size = 3,
              hjust = 0.5,
              vjust = -0.5
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = "Location",
      y = "Relative abundance (%)",
      fill = "FFG"
    ) +
    theme_classic() +
    theme(
      strip.text       = element_text(size = 15),
      strip.background = element_blank(),
      axis.title       = element_text(size = 15),
      axis.text        = element_text(size = 15),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      legend.text      = element_text(size = 15),
      legend.title     = element_text(size = 15),
      legend.key.size  = unit(0.4, "cm")
    )
}


### Alpha diversity / Shannon index 

{
  # Shannon index per Sample
  alpha_df <- df %>%
    group_by(Method, Location, Sample, Taxon) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Method, Location, Sample) %>%
    summarise(Shannon = vegan::diversity(Abundance, index = "shannon"),
              .groups = "drop")
  
  # ANOVA + Tukey + compact letter display
  get_letters_shannon <- function(method_data, method_name) {
    cat("\n========================================\n")
    cat("Method:", method_name, "\n")
    cat("========================================\n")
    
    aov_result <- aov(Shannon ~ Location, data = method_data)
    
    cat("\n--- ANOVA Summary ---\n")
    print(summary(aov_result))
    
    tukey <- TukeyHSD(aov_result)
    cat("\n--- Tukey HSD Post-hoc ---\n")
    print(tukey)
    
    tukey_loc <- tukey$Location
    locs      <- unique(method_data$Location)
    
    # Build p-value matrix
    p_matrix <- matrix(1, nrow = length(locs), ncol = length(locs),
                       dimnames = list(locs, locs))
    
    for (i in seq_len(nrow(tukey_loc))) {
      pair  <- strsplit(rownames(tukey_loc)[i], "-")[[1]]
      p_val <- tukey_loc[i, "p adj"]
      p_matrix[pair[1], pair[2]] <- p_val
      p_matrix[pair[2], pair[1]] <- p_val
    }
    
    # Assign compact letters
    sig_diff       <- p_matrix < 0.05
    letters_vec    <- setNames(rep(NA, length(locs)), locs)
    current_letter <- 1
    
    for (loc in locs) {
      if (is.na(letters_vec[loc])) {
        letters_vec[loc] <- letters[current_letter]
        current_letter   <- current_letter + 1
      }
      for (other in locs) {
        if (loc != other && !sig_diff[loc, other] && is.na(letters_vec[other])) {
          letters_vec[other] <- letters_vec[loc]
        }
      }
    }
    
    label_df <- method_data %>%
      group_by(Location) %>%
      summarise(y_pos = max(Shannon) * 1.1, .groups = "drop") %>%
      mutate(label = toupper(letters_vec[Location]))
    
    return(label_df)
  }
  
  letter_df <- alpha_df %>%
    group_by(Method) %>%
    group_modify(~ get_letters_shannon(.x, .y$Method)) %>%
    ungroup()
  
  letter_df$Method_label <- c(rep("A) Net",3),
                              rep("B) Sampler",3)
                              )

  alpha_df$Method_label <- c(rep("A) Net",12),
                             rep("B) Sampler",12)
                             )
    
  #Plot
  ggplot(alpha_df, aes(x = Location, y = Shannon, fill = Method)) +
    geom_boxplot(
      width         = 0.6,
      color         = "black",
      linewidth     = 0.4,
      outlier.shape = 21,
      outlier.size  = 1.5,
      outlier.fill  = "white"
    ) +
    geom_text(
      data        = letter_df,
      aes(x = Location, y = y_pos, label = label),
      inherit.aes = FALSE,
      size        = 4,
      position    = position_nudge(x = 0.3)
    ) +
    facet_wrap(~ Method_label) +
    scale_fill_manual(values = c("Net" = "#F4A5A5", "Sampler" = "#85C7E0")) +
    labs(
      x    = "Location",
      y    = "Shannon diversity (H')",
      fill = "Method"
    ) +
    theme_classic() +
    theme(
      strip.text       = element_text(size = 15),
      strip.background = element_blank(),
      axis.title       = element_text(size = 15),
      axis.text        = element_text(size = 15),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      legend.position   = "none"
    )


}

{
  
  # Calculate shannon index per sample, add to copy of df named Shannon
#   Moet dit per sample of per location/method? 
  
#   Per Sample
  Shannon <- df |>
    dplyr::mutate(shannon_index_Sample = vegan::diversity(Abundance), .by = ID)
  
#   Per Location
  Shannon <- Shannon |>
    dplyr::mutate(shannon_index_Location = vegan::diversity(Abundance), .by = Location)
  
#   Per location + Method
  Shannon <- Shannon |>
    dplyr::mutate(shannon_index_LocmMet = vegan::diversity(Abundance), .by = LocMet)
  
#   Store Shannon index per Location + Method
  Location <- c("Drained", "Paludiculture", "Restored", "Drained", "Paludiculture", "Restored")
  
  Method <- c(rep("Net",3), rep("Sampler",3))
  
  Shannon_index <- unique(Shannon$shannon_index_LocmMet)

  Total_Abundance_df <- df %>%
  group_by(LocMet) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE))
  Total_Abundance <- Total_Abundance_df$Total_Abundance
  
  Unique_Taxa_df <- df %>%
  group_by(LocMet) %>%
  summarise(Unique_Taxa = n_distinct(Taxon))
  Unique_Taxa <- Unique_Taxa_df$Unique_Taxa
  
  Table_Shannon <- data.frame(Location, Method,Total_Abundance,Unique_Taxa,Shannon_index)
  
#   Check for normality in each LocMet using shapiro wilk test
  Shannon %>%
    group_by(LocMet) %>%
    summarise(p_value = shapiro.test(shannon_index_LocmMet)$p.value, .groups = "drop")
#   Returns error - no normality because of groups to small
  

#Kruskal-Willis non-parametric test Between locations, per Method

    # Net only - differences between locations
  kruskal.test(shannon_index_LocmMet ~ Location,
               data = Shannon %>% 
                 filter(Method == "Net")
               )
#   Calculate the number of rows (N-value)
  Shannon %>% filter(Method == "Net") %>% nrow()
  
  dunn_test(Shannon %>% 
            filter(Method == "Net"),
            shannon_index_LocmMet ~ Location,
            p.adjust.method = "bonferroni"
            )
  
  # Sampler only - differences between locations
  kruskal.test(shannon_index_LocmMet ~ Location,
               data = Shannon %>% 
                 filter(Method == "Sampler")
               )
               
#   Calculate the number of rows (N-value)
 Shannon %>% filter(Method == "Sampler") %>% nrow()
               
  dunn_test(Shannon %>% filter
            (Method == "Sampler"),
            shannon_index_LocmMet ~ Location,
            p.adjust.method = "bonferroni"
            )
  
}


# Shannon per functional feeding group

{
  # Shannon index per sample by FFG
  alpha_df_ffg <- df_ffg %>%
    group_by(Method, Location, Sample, FFG) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Method, Location, Sample) %>%
    summarise(
      Shannon = vegan::diversity(Abundance, index = "shannon"),
      .groups = "drop"
    )
  
  get_letters_shannon <- function(method_data, method_name) {
    cat("\n========================================\n")
    cat("Method:", method_name, "\n")
    cat("========================================\n")
    
    aov_result <- aov(Shannon ~ Location, data = method_data)
    
    cat("\n--- ANOVA Summary ---\n")
    print(summary(aov_result))
    
    tukey <- TukeyHSD(aov_result)
    cat("\n--- Tukey HSD Post-hoc ---\n")
    print(tukey)
    
    tukey_loc <- tukey$Location
    locs      <- unique(method_data$Location)
    
    # Build p-value matrix
    p_matrix <- matrix(1, nrow = length(locs), ncol = length(locs),
                       dimnames = list(locs, locs))
    
    for (i in seq_len(nrow(tukey_loc))) {
      pair  <- strsplit(rownames(tukey_loc)[i], "-")[[1]]
      p_val <- tukey_loc[i, "p adj"]
      p_matrix[pair[1], pair[2]] <- p_val
      p_matrix[pair[2], pair[1]] <- p_val
    }
    
    # Assign compact letters
    sig_diff       <- p_matrix < 0.05
    letters_vec    <- setNames(rep(NA, length(locs)), locs)
    current_letter <- 1
    
    for (loc in locs) {
      if (is.na(letters_vec[loc])) {
        letters_vec[loc] <- letters[current_letter]
        current_letter   <- current_letter + 1
      }
      for (other in locs) {
        if (loc != other && !sig_diff[loc, other] && is.na(letters_vec[other])) {
          letters_vec[other] <- letters_vec[loc]
        }
      }
    }
    
    label_df <- method_data %>%
      group_by(Location) %>%
      summarise(y_pos = max(Shannon) * 1.1, .groups = "drop") %>%
      mutate(label = toupper(letters_vec[Location]))
    
    return(label_df)
  }
  
  # Letters for FFG
  letter_df_ffg <- alpha_df_ffg %>%
    group_by(Method) %>%
    group_modify(~ get_letters_shannon(.x, .y$Method)) %>%
    ungroup()
  
  # Method labels
  letter_df_ffg$Method_label <- c(rep("A) Net", 3),
                                  rep("B) Sampler", 3))
  
  alpha_df_ffg$Method_label <- c(rep("A) Net", 12),
                                 rep("B) Sampler", 12))
  
  # Plot Shannon diversity per FFG
  ggplot(alpha_df_ffg, aes(x = Location, y = Shannon, fill = Method)) +
    geom_boxplot(
      width         = 0.6,
      color         = "black",
      linewidth     = 0.4,
      outlier.shape = 21,
      outlier.size  = 1.5,
      outlier.fill  = "white"
    ) +
    geom_text(
      data        = letter_df_ffg,
      aes(x = Location, y = y_pos, label = label),
      inherit.aes = FALSE,
      size        = 4,
      position    = position_nudge(x = 0.3)
    ) +
    facet_wrap(~ Method_label) +
    scale_fill_manual(values = c("Net" = "#F4A5A5", "Sampler" = "#85C7E0")) +
    labs(
      x    = "Location",
      y    = "Shannon diversity (H')",
      fill = "Method"
    ) +
    theme_classic() +
    theme(
      strip.text       = element_text(size = 15),
      strip.background = element_blank(),
      axis.title       = element_text(size = 15),
      axis.text        = element_text(size = 15),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      legend.position  = "none"
    )
}

{
  # Per location + Method (FFG-based)
  Shannon_FFG <- df_ffg |>
    dplyr::mutate(shannon_index_LocmMet_FFG = vegan::diversity(Abundance), .by = LocMet)
  
  # Store Shannon index per Location + Method
  Location_FFG <- c("Drained", "Paludiculture", "Restored", "Drained", "Paludiculture", "Restored")
  
  Method_FFG <- c(rep("Net", 3), rep("Sampler", 3))
  
  Shannon_index_FFG <- unique(Shannon_FFG$shannon_index_LocmMet_FFG)
  
  Total_Abundance_FFG_df <- df_ffg %>%
    group_by(LocMet) %>%
    summarise(Total_Abundance = sum(Abundance, na.rm = TRUE))
  Total_Abundance_FFG <- Total_Abundance_FFG_df$Total_Abundance
  
  Unique_FFG_df <- df_ffg %>%
    group_by(LocMet) %>%
    summarise(Unique_FFG = n_distinct(FFG))
  Unique_FFG <- Unique_FFG_df$Unique_FFG
  
  Table_Shannon_FFG <- data.frame(Location_FFG, Method_FFG, Total_Abundance_FFG, Unique_FFG, Shannon_index_FFG)
  
  # Check for normality in each LocMet using Shapiro-Wilk test
  Shannon_FFG %>%
    group_by(LocMet) %>%
    summarise(p_value = shapiro.test(shannon_index_LocmMet_FFG)$p.value, .groups = "drop")
  # Returns error - no normality because of groups too small
  
  # Kruskal-Wallis non-parametric test between locations, per Method
  
  # Net only - differences between locations
  kruskal.test(shannon_index_LocmMet_FFG ~ Location,
               data = Shannon_FFG %>%
                 filter(Method == "Net")
  )
  # Calculate the number of rows (N-value)
  Shannon_FFG %>% filter(Method == "Net") %>% nrow()
  
  dunn_test(Shannon_FFG %>%
              filter(Method == "Net"),
            shannon_index_LocmMet_FFG ~ Location,
            p.adjust.method = "bonferroni"
  )
  
  # Sampler only - differences between locations
  kruskal.test(shannon_index_LocmMet_FFG ~ Location,
               data = Shannon_FFG %>%
                 filter(Method == "Sampler")
  )
  # Calculate the number of rows (N-value)
  Shannon_FFG %>% filter(Method == "Sampler") %>% nrow()
  
  dunn_test(Shannon_FFG %>%
              filter(Method == "Sampler"),
            shannon_index_LocmMet_FFG ~ Location,
            p.adjust.method = "bonferroni"
  )
}












