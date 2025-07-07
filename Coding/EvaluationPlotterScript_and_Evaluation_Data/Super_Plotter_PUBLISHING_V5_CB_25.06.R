
library(dplyr)
library(ggplot2)
library(viridis)
install.packages("ggforce")
library(ggforce)
library(grid)
library(gridExtra)
install.packages("glue")
library(glue)

# READ IN 
Master_df_GY <- read.csv("Definitive_GY_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
Master_df_TGW <- read.csv("Definitive_TGW_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
Master_df_GPC <- read.csv("Definitive_GPC_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")
Master_df_HET <- read.csv("Definitive_HET_36_GBM.fit_MASTER_METRICS_DF_27.03.2025.csv")

####  MASTER DATAFRAME FOR PLOTTING ####


Master_df <-  rbind(Master_df_GY,Master_df_GPC, Master_df_TGW, Master_df_HET)


Master_df$Phenotype[Master_df$Phenotype == "GY_1"] <- "Grain Yield"
Master_df$Phenotype[Master_df$Phenotype == "HET_1"] <- "Height to Eartip"
Master_df$Phenotype[Master_df$Phenotype == "GPC_1"] <- "Grain Protein Content"
Master_df$Phenotype[Master_df$Phenotype == "TGW_1"] <- "Thousand Grain Weight"


## 1. Palette ----------------------------------------------------------
lr_levels <- c("0.0025", "0.042", "0.0815", "0.121", "0.1605", "0.2")

# Okabe–Ito palette (modified for distinctness)
lr_cols_cb <- c(
  "#0072B2",  # Blue
  "#E69F00",  # Orange
  "#009E73",  # Bluish Green
  "#CC79A7",  # Reddish Purple
  "#D55E00",  # Vermillion
  "#56B4E9"   # Sky Blue
)

#lr_shapes <- c(16, 17, 15, 23, 4, 8)  # Circle, Triangle, Square, Diamond, X, Asterisk
lr_shapes <- c(21, 24, 22, 23, 4, 8)
pal_col_cb   <- scale_colour_manual(values = setNames(lr_cols_cb, lr_levels),
                                    breaks = lr_levels, name = "Learning rate")
pal_fill_cb  <- scale_fill_manual( values = setNames(lr_cols_cb, lr_levels), guide = "none")
pal_shape_cb <- scale_shape_manual(values = setNames(lr_shapes, lr_levels),
                                   breaks = lr_levels, name = "Learning rate")


# 2. Define an “export high-res” helper
export_hires_4_traits <- function(plot, basename,
                                  width_in  = 8,    # inches
                                  height_in = 5,    # inches
                                  dpi_png   = 600,  # dots per inch for PNG
                                  dpi_tiff  = 600)  # for TIFF if you need
{
  # —— Vector PDF (scales without loss) ——
  ggsave(
    filename = glue("{basename}.pdf"),
    plot     = plot,
    device   = cairo_pdf,     # uses Cairo for embedded fonts, anti-aliasing
    width    = width_in,
    height   = height_in,
    units    = "in",
    bg       = "white"
  )
  
  # —— High-res PNG (for slide decks that need raster) ——
  ggsave(
    filename = glue("{basename}.png"),
    plot     = plot,
    device   = png,           # base R PNG; for even crisper try ragg::agg_png
    width    = width_in,
    height   = height_in,
    units    = "in",
    dpi      = dpi_png,
    bg       = "white"
  )
  
  # —— Optional TIFF (print journals often require TIFF) ——
  ggsave(
    filename = glue("{basename}.tiff"),
    plot     = plot,
    device   = tiff,          # or ragg::agg_tiff
    width    = width_in,
    height   = height_in,
    units    = "in",
    dpi      = dpi_tiff,
    compression = "lzw",
    bg       = "white"
  )
}


#### ICC absolute  ####
# Filter the data for ICC metric
icc_data <- Master_df %>% 
  filter(Metric == "ICC_Prediction_Stability")


# Summarize the data by computing the mean, min, and max for each combination of BoostRounds, LearningRate, and Phenotype
# Summarize the data: calculating mean, standard error, and min/max
summary_df <- icc_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_ICC = mean(Value), # Mean value of ICC - 5 for each combination of BoostRounds, LearningRate
    sd_ICC = sd(Value),
    se_ICC = sd(Value) / sqrt(n()),
    min_ICC = min(Value),
    max_ICC = max(Value),
    .groups = "drop"
  )


# Use the same factor levels you defined
summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
summary_df$Phenotype    <- factor(summary_df$Phenotype, levels = c(
  "Grain Yield",
  "Thousand Grain Weight",
  "Grain Protein Content",
  "Height to Eartip"
))
icc_data$LearningRate   <- factor(as.character(icc_data$LearningRate), levels = lr_levels)
icc_data$Phenotype      <- factor(icc_data$Phenotype, levels = levels(summary_df$Phenotype))

# 1. Get original levels
phen_levels <- levels(summary_df$Phenotype)

# 2. Build labels "(a) Grain Yield", "(b) Thousand Grain Weight", …
phen_labels <- setNames(
  paste0("(", letters[seq_along(phen_levels)], ") ", phen_levels),
  phen_levels
)

icc_4_traits_absolute_fixed_y <- ggplot() +
  # summary lines & errorbars
  geom_line(data = summary_df,
            aes(x      = BoostRounds,
                y      = mean_ICC,
                colour = LearningRate),
            linewidth = .65,
            linetype  = "solid",
            alpha     = 0.7) +
  geom_errorbar(data = summary_df,
                aes(x      = BoostRounds,
                    ymin   = mean_ICC - se_ICC,
                    ymax   = mean_ICC + se_ICC,
                    colour = LearningRate),
                width = 100) +

  geom_point(data = icc_data,
             aes(x      = BoostRounds,
                 y      = Value,
                 colour = LearningRate,
                 shape  = LearningRate),
             alpha    = 0.6,
             size     = 0.8) + 
  facet_wrap(~ Phenotype,
             labeller = labeller(Phenotype = as_labeller(phen_labels))
  ) +
  pal_col_cb +
  pal_shape_cb +
  pal_fill_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "ICC",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# Export the high-res plot
export_hires_4_traits(icc_4_traits_absolute_fixed_y,
                      "ICC_vs_BoostRounds_cb_4_traits_fixed_y",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


install.packages("egg")
library(egg)
# Then tag facets:
p_tagged <- tag_facet(
  icc_4_traits_absolute_fixed_y,
  open     = "(",       # opening bracket
  close    = ")",       # closing bracket
  tag_pool = letters,   # use letters a, b, c, ...
  x        = -Inf,      # left edge
  y        =  Inf,      # top edge
  hjust    = -0.1,      
  vjust    =  1.2,
  fontface =  "bold"
)

# Display:
p_tagged

# Summarize the data by computing the mean ICC for each combination of BoostRounds, LearningRate, and Phenotype




#### ICC normalised #### -------------------------------------------------------

# 2. Normalise 0–1 within each Phenotype ---------------------------
icc_data_norm <- icc_data %>%               # start from raw replicates
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup()

summary_df_norm <- icc_data_norm %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarise(
    mean_ICC = mean(Value_norm),
    sd_ICC   = sd(Value_norm),
    se_ICC   = sd_ICC / sqrt(n()),
    min_ICC  = min(Value_norm),
    max_ICC  = max(Value_norm),
    .groups  = "drop"
  )

# (keep factor levels consistent)
summary_df_norm$LearningRate <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)
icc_data_norm$LearningRate   <- factor(as.character(icc_data_norm$LearningRate), levels = lr_levels)

# 3. Plot normalised values (0–1 axis shared by definition) --------
icc_4_traits_norm <- ggplot() +
  geom_line(data = summary_df_norm,
            aes(x = BoostRounds, y = mean_ICC, colour = LearningRate),
            linewidth = .65, linetype = "solid", alpha = 0.7) +
  geom_errorbar(data = summary_df_norm,
                aes(x = BoostRounds,
                    ymin = mean_ICC - se_ICC,
                    ymax = mean_ICC + se_ICC,
                    colour = LearningRate),
                width = 100) +
  geom_point(data = icc_data_norm,
             aes(x = BoostRounds, y = Value_norm,
                 colour = LearningRate, shape = LearningRate),
             alpha = 0.6, size = 0.8) +
  facet_wrap(~ Phenotype,
             labeller = labeller(Phenotype = as_labeller(phen_labels))
  ) +             # fixed-y (0–1) across all facets
  pal_col_cb + pal_shape_cb + pal_fill_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "BoostRounds", y = "ICC (normalized 0–1)",
       colour = "Learning Rate", shape = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(icc_4_traits_norm,
                      "ICC_vs_BoostRounds_cb_4_traits_norm",
                      width_in = 12, height_in = 7.5, dpi_png = 600)










#### FLEISS KAAPPA ####

# Filter and summarize for Fleiss Kappa metric
FLEISS_data <- Master_df %>% 
  filter(Metric == "Fleiss_Kappa_Prediction_Stability")

# FLEISS_data$Phenotype[FLEISS_data$Phenotype == "GY_1"] <- "Grain Yield"
# FLEISS_data$Phenotype[FLEISS_data$Phenotype == "HET_1"] <- "Height to Eartip"
# FLEISS_data$Phenotype[FLEISS_data$Phenotype == "GPC_1"] <- "Grain Protein Content"
# FLEISS_data$Phenotype[FLEISS_data$Phenotype == "TGW_1"] <- "Thousand Grain Weight"

summary_df <- FLEISS_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_FLEISS = mean(Value),
    sd_FLEISS = sd(Value),
    se_FLEISS = sd(Value) / sqrt(n()),
    min_FLEISS = min(Value),
    max_FLEISS = max(Value),
    .groups = "drop"
  )

# Adjusted Plot 1: Fleiss Kappa vs BoostRounds by LearningRate
ggplot(summary_df, aes(x = BoostRounds, y = mean_FLEISS, color = as.factor(LearningRate))) +
  geom_line(linewidth = 0.05, alpha = 0.5, linetype = "dashed") +  
  geom_errorbar(aes(
    ymin = mean_FLEISS - se_FLEISS,
    ymax = mean_FLEISS + se_FLEISS
  ), width = 100) +
  geom_smooth(data = FLEISS_data, 
              aes(x = as.numeric(BoostRounds), y = Value, color = as.factor(LearningRate), 
                  fill = as.factor(LearningRate)),
              method = "loess", level = 0.95, span = 1.0,
              linetype = "dotted", alpha = 0.1, linewidth = 0.5, show.legend = FALSE) +
  geom_point(data = FLEISS_data, 
             aes(x = BoostRounds, y = Value, color = as.factor(LearningRate)),
             position = position_jitter(width = 0.1), alpha = 0.8, size = 0.75) +
  facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("Grain Yield", "Thousand Grain Weight", "Grain Protein Content", "Height to Eartip"))) + # 
  labs(title = "FLEISS vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "FLEISS",
       color = "Learning Rate") +
  theme_minimal()



# Use the same factor levels you defined
summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
summary_df$Phenotype    <- factor(summary_df$Phenotype, levels = c(
  "Grain Yield",
  "Thousand Grain Weight",
  "Grain Protein Content",
  "Height to Eartip"
))
FLEISS_data$LearningRate <- factor(as.character(FLEISS_data$LearningRate), levels = lr_levels)
FLEISS_data$Phenotype    <- factor(FLEISS_data$Phenotype, levels = levels(summary_df$Phenotype))

# 1. Get original levels
phen_levels <- levels(summary_df$Phenotype)

# 2. Build labels "(a) Grain Yield", "(b) Thousand Grain Weight", …
phen_labels <- setNames(
  paste0("(", letters[seq_along(phen_levels)], ") ", phen_levels),
  phen_levels
)

fleiss_4_traits <- ggplot() +
  # summary lines & errorbars
  geom_line(data = summary_df,
            aes(x      = BoostRounds,
                y      = mean_FLEISS,
                colour = LearningRate),
            #position = position_dodge(width = 500),
            linewidth = .65,
            linetype  = "solid",
            alpha     = 0.7) +
  geom_errorbar(data = summary_df,
                aes(x      = BoostRounds,
                    ymin   = mean_FLEISS - se_FLEISS,
                    ymax   = mean_FLEISS + se_FLEISS,
                    colour = LearningRate),
                width    = 100) + #,position = position_dodge(width = 500)) +
                
  geom_point(data = summary_df,
             aes(x      = BoostRounds,
                 y      = mean_FLEISS,
                 colour = LearningRate,
                 shape  = LearningRate),
             size     = 1.5,
             stroke   = 0.5) + # ,   position = position_dodge(width = 500)
  # raw replicate points and LOESS smooth
  geom_point(data = FLEISS_data,
             aes(x      = BoostRounds,
                 y      = Value,
                 colour = LearningRate),
            # position = position_jitter(width = 100),
             alpha    = 0.6,
             size     = 0.8,
             show.legend = FALSE) +
  # geom_smooth(data = FLEISS_data,
  #             aes(x      = BoostRounds,
  #                 y      = Value,
  #                 colour = LearningRate,
  #                 fill   = LearningRate),
  #             method    = "loess",
  #             se        = TRUE,
  #             span      = 1.0,
  #             linetype  = "dotted",
  #             alpha     = 0.1,
  #             linewidth = 0.5,
  #             show.legend = FALSE) +
  facet_wrap(~ Phenotype,
             labeller = labeller(Phenotype = as_labeller(phen_labels))
  ) + #, scales = "free_y"
  pal_col_cb +
  pal_shape_cb +
  pal_fill_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "Fleiss Kappa vs BoostRounds by LearningRate",
       x      = "BoostRounds",
       y      = "Fleiss Kappa",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# Export the high-res plot
export_hires_4_traits(fleiss_4_traits,
                      "Fleiss_vs_BoostRounds_cb",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


#### FLEISS normalised #### -------------------------------------------------------
# 1. Normalise 0–1 within each Phenotype
FLEISS_data_norm <- FLEISS_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup()

# 2. Summarise normalised data: mean, SE, etc.
summary_df_norm <- FLEISS_data_norm %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarise(
    mean_FLEISS = mean(Value_norm),
    sd_FLEISS   = sd(Value_norm),
    se_FLEISS   = sd_FLEISS / sqrt(n()),
    min_FLEISS  = min(Value_norm),
    max_FLEISS  = max(Value_norm),
    .groups      = "drop"
  )

# 3. Keep factor-level consistency
summary_df_norm$LearningRate <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)
summary_df_norm$Phenotype    <- factor(summary_df_norm$Phenotype, levels = c(
  "Grain Yield",
  "Thousand Grain Weight",
  "Grain Protein Content",
  "Height to Eartip"
))
FLEISS_data_norm$LearningRate <- factor(as.character(FLEISS_data_norm$LearningRate), levels = lr_levels)
FLEISS_data_norm$Phenotype    <- factor(FLEISS_data_norm$Phenotype, levels = levels(summary_df_norm$Phenotype))


# 1. Get original levels
phen_levels <- levels(summary_df_norm$Phenotype)

# 2. Build labels "(a) Grain Yield", "(b) Thousand Grain Weight", …
phen_labels <- setNames(
  paste0("(", letters[seq_along(phen_levels)], ") ", phen_levels),
  phen_levels
)

# 4. Plot normalised Fleiss Kappa (0–1) 
fleiss_4_traits_norm <- ggplot() +
  geom_line(data = summary_df_norm,
            aes(x      = BoostRounds,
                y      = mean_FLEISS,
                colour = LearningRate),
            linewidth = .65,
            linetype  = "solid",
            alpha     = 0.7) +
  geom_errorbar(data = summary_df_norm,
                aes(x      = BoostRounds,
                    ymin   = mean_FLEISS - se_FLEISS,
                    ymax   = mean_FLEISS + se_FLEISS,
                    colour = LearningRate),
                width    = 100) +
  geom_point(data = FLEISS_data_norm,
             aes(x      = BoostRounds,
                 y      = Value_norm,
                 colour = LearningRate,
                 shape  = LearningRate),
             alpha    = 0.6,
             size     = 0.8) +
  facet_wrap(~ Phenotype,
             labeller = labeller(Phenotype = as_labeller(phen_labels))
  ) +
  pal_col_cb + pal_shape_cb + pal_fill_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "Fleiss Kappa (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 5. Export the high-res normalized plot
export_hires_4_traits(fleiss_4_traits_norm,
                      "Fleiss_vs_BoostRounds_cb_4_traits_norm",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


# # Filter the data for FLEISS metric
# FLEISS_data <- Master_df %>% 
#   filter(Metric == "Fleiss_Kappa_Prediction_Stability")
# 
# # Summarize the data by computing the mean FLEISS for each combination of BoostRounds, LearningRate, and Phenotype
# summary_df <- FLEISS_data %>%
#   group_by(BoostRounds, LearningRate, Phenotype) %>%
#   summarize(
#     mean_FLEISS = mean(Value),
#     .groups = "drop"
#   )
# 
# # Create a heatmap: x-axis is BoostRounds, y-axis is LearningRate, fill is mean FLEISS.
# ggplot(summary_df, aes(x = BoostRounds, y = as.factor(LearningRate), fill = mean_FLEISS)) +
#   geom_tile() +
#   facet_wrap(~Phenotype) +
#   scale_fill_viridis_c(option = "plasma", name = "Mean FLEISS") +
#   labs(title = "Heatmap: FLEISS vs BoostRounds by LearningRate",
#        x = "Boost Rounds",
#        y = "Learning Rate") +
#   theme_minimal()




#### Pearson ####


# Filter the data for ICC metric
pearson_data <- Master_df %>% 
  filter(Metric == "Pearson")

# pearson_data$Phenotype[pearson_data$Phenotype == "GY_1"] <- "Grain Yield"
# pearson_data$Phenotype[pearson_data$Phenotype == "HET_1"] <- "Height to Eartip"
# pearson_data$Phenotype[pearson_data$Phenotype == "GPC_1"] <- "Grain Protein Content"
# pearson_data$Phenotype[pearson_data$Phenotype == "TGW_1"] <- "Thousand Grain Weight"

# Summarize the data by computing the mean, min, and max for each combination of BoostRounds, LearningRate, and Phenotype
# Summarize the data: calculating mean, standard error, and min/max
summary_df <- pearson_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_pearson = mean(Value),
    sd_pearson = sd(Value),
    se_pearson = sd(Value) / sqrt(n()),
    min_pearson = min(Value),
    max_pearson = max(Value),
    count =  n(),
    .groups = "drop"
  ) 
summary(summary_df$count)
summary_df_subset <- pearson_data %>%
  group_by(BoostRounds, LearningRate, Phenotype, Subset) %>%
  summarize(
    mean_pearson_subset = mean(Value),
    sd_pearson_subset = sd(Value),
    se_pearson_subset = sd(Value) / sqrt(n()),
    min_pearson_subset = min(Value),
    max_pearson_subset = max(Value),
    .groups = "drop"
  )


# summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("Grain Yield", "Thousand Grain Weight", "Grain Protein Content", "Height to Eartip"))
# 
# #### PRESENT ME
# ggplot(summary_df, aes(x = BoostRounds, y = mean_pearson, color = as.factor(LearningRate))) +
#   geom_line(position = position_dodge(width = 50)) +  # Line connecting the mean values
#   geom_errorbar(aes(ymin = mean_pearson - se_pearson, ymax = mean_pearson + se_pearson), width = 500, position = position_dodge(width = 50)) +  # Error bars for the range of values
#   geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line (line of best fit)
#     facet_wrap(~Phenotype, scale = "free_y") + 
#   labs(title = "Pearson Correlation vs BoostRounds by LearningRate",
#        x = "Boost Rounds",
#        y = "Pearson Correlation",
#        color = "Learning Rate",
#        subtitle = "Mean Pearson Correlation per Subset (5)") +
#   theme_minimal()




# Use the same factor levels you defined
summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("Grain Yield", "Thousand Grain Weight", "Grain Protein Content", "Height to Eartip"))

# 1. Get original levels
phen_levels <- levels(summary_df$Phenotype)

# 2. Build labels "(a) Grain Yield", "(b) Thousand Grain Weight", …
phen_labels <- setNames(
  paste0("(", letters[seq_along(phen_levels)], ") ", phen_levels),
  phen_levels
)

pearson_4_traits <-  ggplot(summary_df,
       aes(x = BoostRounds,
           y = mean_pearson,
           colour = LearningRate,      # colour by LR
           shape  = LearningRate)) +   # shape by LR
  geom_line(position = position_dodge(width = 500) ,  linewidth = .65,
            aes(linetype = NULL)) +  # you could also map linetype if needed
  geom_errorbar(aes(ymin = mean_pearson - se_pearson,
                    ymax = mean_pearson + se_pearson),
                width = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5, position = position_dodge(width = 500)) +
  # geom_smooth(method = "loess",
  #             se = FALSE,
  #             linetype = "dotted",
  #             alpha = 0.1,
  #             linewidth = 0.5, position = position_dodge(width = 200)) +
  #facet_wrap(~ Phenotype, scales = "free_y") +
  
  facet_wrap(~ Phenotype,   scales = "free_y",
             labeller = labeller(Phenotype = as_labeller(phen_labels))
  ) +
  # your custom CB‐friendly scales:
  pal_col_cb +
  pal_shape_cb +
  # unify colour + shape into one legend:
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title    = "Pearson Correlation vs BoostRounds by LearningRate",
       x        = "BoostRounds",
       y        = "Pearson Correlation",
       colour   = "Learning Rate",
       shape    = "Learning Rate") +
       # subtitle = "Mean Pearson Correlation per Subset (5)") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )




#  Export  plot
export_hires_4_traits(pearson_4_traits, "Pearson_vs_BoostRounds_cb_free_y",
             width_in  = 12,     # tweak to your journal or slide ratio
             height_in = 7.5,
             dpi_png   = 600)


#### Pearson absolute fixed y #### -------------------------------------------------

# 1. Compute common y-limits from your summary (min/max across all traits)
# 2. Plot with fixed y-axis (same scale across facets)
pearson_4_traits_absolute_fixed_y <- ggplot(
  summary_df,
  aes(x      = BoostRounds,
      y      = mean_pearson,
      colour = LearningRate,
      shape  = LearningRate)
) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_pearson - se_pearson,
                    ymax = mean_pearson + se_pearson),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  #facet_wrap(~ Phenotype) +                   # fixed y by default
  facet_wrap(~ Phenotype,
             labeller = labeller(Phenotype = as_labeller(phen_labels))
  ) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "Pearson Correlation",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 3. Export
export_hires_4_traits(pearson_4_traits_absolute_fixed_y,
                      "Pearson_vs_BoostRounds_cb_4_traits_fixed_y",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


#### Pearson normalised #### -------------------------------------------------------

# 1. Normalise 0–1 within each Phenotype
pearson_data_norm <- pearson_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup()

# 2. Summarise normalised data
summary_df_norm <- pearson_data_norm %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarise(
    mean_pearson = mean(Value_norm),
    sd_pearson   = sd(Value_norm),
    se_pearson   = sd_pearson / sqrt(n()),
    min_pearson  = min(Value_norm),
    max_pearson  = max(Value_norm),
    .groups       = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate   <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)
summary_df_norm$Phenotype      <- factor(as.character(summary_df_norm$Phenotype),
                                         levels = c("Grain Yield",
                                                    "Thousand Grain Weight",
                                                    "Grain Protein Content",
                                                    "Height to Eartip"))
pearson_data_norm$LearningRate <- factor(as.character(pearson_data_norm$LearningRate), levels = lr_levels)
pearson_data_norm$Phenotype    <- factor(as.character(pearson_data_norm$Phenotype),
                                         levels(summary_df_norm$Phenotype))

# 1. Get original levels
phen_levels <- levels(summary_df_norm$Phenotype)

# 2. Build labels "(a) Grain Yield", "(b) Thousand Grain Weight", …
phen_labels <- setNames(
  paste0("(", letters[seq_along(phen_levels)], ") ", phen_levels),
  phen_levels
)

# 4. Plot normalised Pearson
pearson_4_traits_norm <- ggplot(
  summary_df_norm,
  aes(x      = BoostRounds,
      y      = mean_pearson,
      colour = LearningRate,
      shape  = LearningRate)
) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_pearson - se_pearson,
                    ymax = mean_pearson + se_pearson),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
    facet_wrap(~ Phenotype,              labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +                   # 0–1 axis shared by default
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "Pearson Correlation (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 5. Export
export_hires_4_traits(pearson_4_traits_norm,
                      "Pearson_vs_BoostRounds_cb_4_traits_norm",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


#### R Squared ####

# Filter the data for R_squared metric
r2_data <- Master_df %>% 
  filter(Metric == "R_squared") 

r2_data$Phenotype[r2_data$Phenotype == "GY_1"] <- "Grain Yield"
r2_data$Phenotype[r2_data$Phenotype == "HET_1"] <- "Height to Eartip"
r2_data$Phenotype[r2_data$Phenotype == "GPC_1"] <- "Grain Protein Content"
r2_data$Phenotype[r2_data$Phenotype == "TGW_1"] <- "Thousand Grain Weight"

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- r2_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_r2 = mean(Value),
    sd_r2 = sd(Value),
    se_r2 = sd(Value) / sqrt(n()),
    min_r2 = min(Value),
    max_r2 = max(Value),
    .groups = "drop"
  )

# summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("Grain Yield", "Thousand Grain Weight", "Grain Protein Content", "Height to Eartip"))
# 
# # Plot the data with error bars and a trend line
# ggplot(summary_df, aes(x = BoostRounds, y = mean_r2, color = as.factor(LearningRate))) +
#   geom_line() +  # Line connecting the mean R_squared values
#   geom_errorbar(aes(ymin = mean_r2 - se_r2, ymax = mean_r2 + se_r2), width = 100) +  # Error bars for the range of values (optional)
#   geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line (line of best fit)
#   geom_point(data = summary_df, aes(x = BoostRounds, y = mean_r2, color = as.factor(LearningRate))) + 
#   facet_wrap(~Phenotype, scale = "free_y") +
#   labs(title = "R-squared vs BoostRounds by LearningRate",
#        x = "Boost Rounds",
#        y = "R-squared",
#        color = "Learning Rate") +
#   theme_minimal()



# Use the same factor levels you defined
summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
summary_df$Phenotype    <- factor(summary_df$Phenotype, levels = c(
  "Grain Yield",
  "Thousand Grain Weight",
  "Grain Protein Content",
  "Height to Eartip"
))

r2_4_traits <- ggplot(summary_df,
                      aes(x      = BoostRounds,
                          y      = mean_r2,
                          colour = LearningRate,
                          shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_r2 - se_r2,
                    ymax = mean_r2 + se_r2),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  # (optional) trend line, if you want it uncommented:
  # geom_smooth(method     = "loess",
  #             se         = FALSE,
  #             linetype   = "dotted",
  #             alpha      = 0.1,
  #             linewidth  = 0.5,
  #             position   = position_dodge(width = 200)) +
  facet_wrap(~ Phenotype,  scales = "free_y",             labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "R-squared vs BoostRounds by LearningRate",
       x      = "BoostRounds",
       y      = expression(R^2),
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# Export the high-res plot
export_hires_4_traits(r2_4_traits,
                      "R2_vs_BoostRounds_cb_free_y",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)



#### R Squared absolute fixed y #### -------------------------------------------------

# (Assumes you’ve already created `r2_data` and `summary_df` and set factor levels as above.)

r2_4_traits_absolute_fixed_y <- ggplot(
  summary_df,
  aes(x      = BoostRounds,
      y      = mean_r2,
      colour = LearningRate,
      shape  = LearningRate)
) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_r2 - se_r2,
                    ymax = mean_r2 + se_r2),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
    facet_wrap(~ Phenotype,              labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +    # fixed y across facets
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = expression(R^2),
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# Export
export_hires_4_traits(r2_4_traits_absolute_fixed_y,
                      "R2_vs_BoostRounds_cb_4_traits_fixed_y",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


#### R Squared normalised #### -------------------------------------------------------

# 1. Normalise 0–1 within each Phenotype
r2_data_norm <- r2_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup()

# 2. Summarise normalised data
summary_df_norm <- r2_data_norm %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarise(
    mean_r2 = mean(Value_norm),
    sd_r2   = sd(Value_norm),
    se_r2   = sd_r2 / sqrt(n()),
    min_r2  = min(Value_norm),
    max_r2  = max(Value_norm),
    .groups = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate   <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)
summary_df_norm$Phenotype      <- factor(as.character(summary_df_norm$Phenotype),
                                         levels = c("Grain Yield",
                                                    "Thousand Grain Weight",
                                                    "Grain Protein Content",
                                                    "Height to Eartip"))
r2_data_norm$LearningRate      <- factor(as.character(r2_data_norm$LearningRate), levels = lr_levels)
r2_data_norm$Phenotype         <- factor(as.character(r2_data_norm$Phenotype),
                                         levels(summary_df_norm$Phenotype))

# 4. Plot normalised R²
r2_4_traits_norm <- ggplot(
  summary_df_norm,
  aes(x      = BoostRounds,
      y      = mean_r2,
      colour = LearningRate,
      shape  = LearningRate)
) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_r2 - se_r2,
                    ymax = mean_r2 + se_r2),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
    facet_wrap(~ Phenotype,              labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +    # 0–1 axis shared by default
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "R² (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# Export
export_hires_4_traits(r2_4_traits_norm,
                      "R2_vs_BoostRounds_cb_4_traits_norm",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)



#### AUC ####

# Filter the data for AUC metric
auc_data <- Master_df %>% 
  filter(Metric == "AUC")



# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- auc_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_auc = mean(Value),
    sd_auc = sd(Value),
    se_auc = sd(Value) / sqrt(n()),
    min_auc = min(Value),
    max_auc = max(Value),
    .groups = "drop"
  )

summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("Grain Yield", "Thousand Grain Weight", "Grain Protein Content", "Height to Eartip"))


# Plot the data with error bars and a trend line (Plot 1)
ggplot(summary_df, aes(x = BoostRounds, y = mean_auc, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean AUC values
  geom_errorbar(aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc), width = 100) +  # Error bars for the range of values (optional)
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_auc, color = as.factor(LearningRate))) +
  facet_wrap(~Phenotype, scale = "free_y" ) + 
  labs(title = "AUC vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "AUC",
       color = "Learning Rate") +
  theme_minimal()

# Use the same factor levels you defined
summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
summary_df$Phenotype    <- factor(summary_df$Phenotype, levels = c(
  "Grain Yield",
  "Thousand Grain Weight",
  "Grain Protein Content",
  "Height to Eartip"
))

auc_4_traits <- ggplot(summary_df,
                       aes(x      = BoostRounds,
                           y      = mean_auc,
                           colour = LearningRate,
                           shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_auc - se_auc,
                    ymax = mean_auc + se_auc),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  # (optional) add a LOESS trend line by uncommenting:
  # geom_smooth(method     = "loess",
  #             se         = FALSE,
  #             linetype   = "dotted",
  #             alpha      = 0.1,
  #             linewidth  = 0.5,
  #             position   = position_dodge(width = 200)) +
  facet_wrap(~ Phenotype,  scales = "free_y",             labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "AUC vs BoostRounds by LearningRate",
       x      = "BoostRounds",
       y      = "AUC",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# Export the high-res plot
export_hires_4_traits(auc_4_traits,
                      "AUC_vs_BoostRounds_cb_free_y",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


#### AUC absolute fixed y #### -------------------------------------------------

# (Assumes you’ve already prepared `auc_data` and `summary_df` and set factor levels as above.)

auc_4_traits_absolute_fixed_y <- ggplot(
  summary_df,
  aes(x      = BoostRounds,
      y      = mean_auc,
      colour = LearningRate,
      shape  = LearningRate)
) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_auc - se_auc,
                    ymax = mean_auc + se_auc),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
    facet_wrap(~ Phenotype,              labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +    # fixed y across facets
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "AUC",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(auc_4_traits_absolute_fixed_y,
                      "AUC_vs_BoostRounds_cb_4_traits_fixed_y",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


#### AUC normalised #### -------------------------------------------------------

# 1. Normalise 0–1 within each Phenotype
auc_data_norm <- auc_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup()

# 2. Summarise normalised data
summary_df_norm <- auc_data_norm %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarise(
    mean_auc = mean(Value_norm),
    sd_auc   = sd(Value_norm),
    se_auc   = sd_auc / sqrt(n()),
    min_auc  = min(Value_norm),
    max_auc  = max(Value_norm),
    .groups   = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate   <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)
summary_df_norm$Phenotype      <- factor(as.character(summary_df_norm$Phenotype),
                                         levels = c("Grain Yield",
                                                    "Thousand Grain Weight",
                                                    "Grain Protein Content",
                                                    "Height to Eartip"))
auc_data_norm$LearningRate     <- factor(as.character(auc_data_norm$LearningRate), levels = lr_levels)
auc_data_norm$Phenotype        <- factor(as.character(auc_data_norm$Phenotype),
                                         levels(summary_df_norm$Phenotype))

# 4. Plot normalised AUC (0–1 axis shared by default)
auc_4_traits_norm <- ggplot(
  summary_df_norm,
  aes(x      = BoostRounds,
      y      = mean_auc,
      colour = LearningRate,
      shape  = LearningRate)
) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_auc - se_auc,
                    ymax = mean_auc + se_auc),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
    facet_wrap(~ Phenotype,              labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "AUC (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(auc_4_traits_norm,
                      "AUC_vs_BoostRounds_cb_4_traits_norm",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)




#### F1_Score ####

# Filter the data for F1_Score metric
F1_Score_data <- Master_df %>% 
  filter(Metric == "F1_Score")  


# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- F1_Score_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_F1_Score = mean(Value),
    sd_F1_Score = sd(Value),
    se_F1_Score = sd(Value) / sqrt(n()),
    min_F1_Score = min(Value),
    max_F1_Score = max(Value),
    .groups = "drop"
  )

# Plot the data with error bars and a trend line (Plot 1)
ggplot(summary_df, aes(x = BoostRounds, y = mean_F1_Score, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean F1_Score values
  # geom_errorbar(aes(ymin = min_F1_Score, ymax = max_F1_Score), width = 100) +  # Error bars for the range of values (optional)
  #geom_smooth(method = "loess", se = TRUE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_F1_Score, color = as.factor(LearningRate))) +
  facet_wrap(~Phenotype) +
  labs(title = "F1 Score vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "F1_Score",
       color = "Learning Rate") +
  theme_minimal()




#### NDCG 20% ####


# Filter the data for NDCG_at_20_percent metric
ndcg_data <- Master_df %>% 
  filter(Metric == "NDCG_at_20_percent")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- ndcg_data %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarize(
    mean_ndcg = mean(Value),
    sd_ndcg = sd(Value),
    se_ndcg = sd(Value) / sqrt(n()),
    min_ndcg = min(Value),
    max_ndcg = max(Value),
    .groups = "drop"
  )

summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("Grain Yield", "Thousand Grain Weight", "Grain Protein Content", "Height to Eartip"))

# Plot the data with error bars and a trend line (Plot 1)
ggplot(summary_df, aes(x = BoostRounds, y = mean_ndcg, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean NDCG values
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg, ymax = mean_ndcg + se_ndcg ), width = 100) +  # Error bars for the range (optional)
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_ndcg, color = as.factor(LearningRate))) +
  facet_wrap(~Phenotype, scale = "free_y") + 
  labs(title = "NDCG@20% vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "NDCG@20%",
       color = "Learning Rate") +
  theme_minimal()


# Use the same factor levels you defined
summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
summary_df$Phenotype    <- factor(summary_df$Phenotype, levels = c(
  "Grain Yield",
  "Thousand Grain Weight",
  "Grain Protein Content",
  "Height to Eartip"
))

ndcg_4_traits <- ggplot(summary_df,
                        aes(x      = BoostRounds,
                            y      = mean_ndcg,
                            colour = LearningRate,
                            shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg,
                    ymax = mean_ndcg + se_ndcg),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  # (optional) add a LOESS trend line by uncommenting:
  # geom_smooth(method     = "loess",
  #             se         = FALSE,
  #             linetype   = "dotted",
  #             alpha      = 0.1,
  #             linewidth  = 0.5,
  #             position   = position_dodge(width = 200)) +
  facet_wrap(~ Phenotype,  scales = "free_y",             labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "NDCG@20% vs BoostRounds by LearningRate",
       x      = "BoostRounds",
       y      = "NDCG@20%",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# Export the high-res plot
export_hires_4_traits(ndcg_4_traits,
                      "NDCG20_vs_BoostRounds_cb_free_y",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)

 
#### NDCG@20% absolute fixed y #### -------------------------------------------------

# (Assumes you’ve already prepared `ndcg_data` and `summary_df` and set factor levels as above.)

ndcg_4_traits_absolute_fixed_y <- ggplot(
  summary_df,
  aes(x      = BoostRounds,
      y      = mean_ndcg,
      colour = LearningRate,
      shape  = LearningRate)
) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg,
                    ymax = mean_ndcg + se_ndcg),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
    facet_wrap(~ Phenotype,              labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +    # fixed y across facets
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "NDCG@20%",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(ndcg_4_traits_absolute_fixed_y,
                      "NDCG20_vs_BoostRounds_cb_4_traits_fixed_y",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)


#### NDCG@20% normalised #### -------------------------------------------------------

# 1. Normalise 0–1 within each Phenotype
ndcg_data_norm <- ndcg_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))) %>%
  ungroup()

# 2. Summarise normalised data
summary_df_norm <- ndcg_data_norm %>%
  group_by(BoostRounds, LearningRate, Phenotype) %>%
  summarise(
    mean_ndcg = mean(Value_norm),
    sd_ndcg   = sd(Value_norm),
    se_ndcg   = sd_ndcg / sqrt(n()),
    min_ndcg  = min(Value_norm),
    max_ndcg  = max(Value_norm),
    .groups    = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate   <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)
summary_df_norm$Phenotype      <- factor(as.character(summary_df_norm$Phenotype),
                                         levels = c("Grain Yield",
                                                    "Thousand Grain Weight",
                                                    "Grain Protein Content",
                                                    "Height to Eartip"))
ndcg_data_norm$LearningRate    <- factor(as.character(ndcg_data_norm$LearningRate), levels = lr_levels)
ndcg_data_norm$Phenotype       <- factor(as.character(ndcg_data_norm$Phenotype),
                                         levels(summary_df_norm$Phenotype))

# 4. Plot normalised NDCG@20%
ndcg_4_traits_norm <- ggplot(
  summary_df_norm,
  aes(x      = BoostRounds,
      y      = mean_ndcg,
      colour = LearningRate,
      shape  = LearningRate)
) +
  geom_line(position = position_dodge(width = 500),
            linewidth = .65) +
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg,
                    ymax = mean_ndcg + se_ndcg),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.5,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
    facet_wrap(~ Phenotype,              labeller = labeller(Phenotype = as_labeller(phen_labels))   ) +    # 0–1 axis shared by default
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x      = "BoostRounds",
       y      = "NDCG@20% (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(ndcg_4_traits_norm,
                      "NDCG20_vs_BoostRounds_cb_4_traits_norm",
                      width_in  = 12,
                      height_in = 7.5,
                      dpi_png   = 600)



#### ALL TRAITS ####
# ICC 
# Filter the data for ICC metric
icc_data <- Master_df %>% 
  filter(Metric == "ICC_Prediction_Stability") 



# Summarize the data by computing the mean, min, and max for each combination of BoostRounds, LearningRate, and Phenotype
# Summarize the data: calculating mean, standard error, and min/max
summary_df <- icc_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_ICC = mean(Value), # Mean value of ICC - 5 for each combination of BoostRounds, LearningRate
    sd_ICC = sd(Value),
    se_ICC = sd(Value) / sqrt(n()),
    min_ICC = min(Value),
    max_ICC = max(Value),
    count  = n(),
    .groups = "drop"
  )

summary(summary_df$count) # 20 each as expected (5 icc per trait) 


# Plot the data with error bars and a trend line
pALL_ICC <- ggplot(summary_df, aes(x = BoostRounds, y = mean_ICC, color = as.factor(LearningRate))) +
  geom_line(linewidth = 0.05, alpha= 0.5, linetype = "dashed") +  # Line connecting the mean values
  
  geom_errorbar(aes(
    ymin = mean_ICC -  se_ICC,
    ymax = mean_ICC +  se_ICC
  ), width = 100) +  
  geom_smooth(data = icc_data, aes(x = as.numeric(BoostRounds), y = Value, color = as.factor(LearningRate), fill =as.factor(LearningRate)) ,
              method = "loess", level = 0.95, span = 1.0,    linetype = "dotted", alpha = 0.1, linewidth = 0.5, show.legend = FALSE ) +
  
  geom_point(data = icc_data, aes(x = BoostRounds, y = Value, color = as.factor(LearningRate)), 
             position = position_jitter(width = 0.1), alpha = 0.8, size =0.75) + # adjust geom_point size using 
  # facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  
  labs(title = "ICC vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "ICC",
       color = "Learning Rate") +
  theme_minimal()



summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
p_all_icc <- ggplot(summary_df,
                    aes(x      = BoostRounds,
                        y      = mean_ICC,
                        colour = LearningRate,
                        shape  = LearningRate)) +
  geom_line(linewidth = 0.65,
            aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_ICC - se_ICC,
                    ymax = mean_ICC + se_ICC),
                width    = 100) +
  geom_point(size     = 1.25,
             stroke   = 0.5) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "ICC vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "ICC",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(p_all_icc,    "ICC_ALL_vs_BoostRounds_cb",    dpi_png = 600)

#### ALL TRAITS normalized ICC #### -------------------------------------------------

# 1. Normalise 0–1 across all ICC values
icc_data_norm <- icc_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE)))

# 2. Summarise normalised data by BoostRounds & LearningRate
summary_df_norm <- icc_data_norm %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_ICC = mean(Value_norm),
    sd_ICC   = sd(Value_norm),
    se_ICC   = sd_ICC / sqrt(n()),
    min_ICC  = min(Value_norm),
    max_ICC  = max(Value_norm),
    count    = n(),
    .groups  = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)

# 4. Plot normalized ICC (0–1) across all traits
p_all_icc_norm <- ggplot(summary_df_norm,
                         aes(x      = BoostRounds,
                             y      = mean_ICC,
                             colour = LearningRate,
                             shape  = LearningRate)) +
  geom_line(linewidth = 0.65, aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_ICC - se_ICC,
                    ymax = mean_ICC + se_ICC),
                width = 100) +
  geom_point(size   = 1.25,
             stroke = 0.5) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "ICC (normalized 0–1) vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "ICC (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 5. Export
export_hires_4_traits(p_all_icc_norm,
                      "ICC_ALL_vs_BoostRounds_cb_norm",
                      dpi_png = 600)



# FLEISS

# Filter and summarize for Fleiss Kappa metric
FLEISS_data <- Master_df %>% 
  filter(Metric == "Fleiss_Kappa_Prediction_Stability")

summary_df <- FLEISS_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_FLEISS = mean(Value),
    sd_FLEISS = sd(Value),
    se_FLEISS = sd(Value) / sqrt(n()),
    min_FLEISS = min(Value),
    max_FLEISS = max(Value),
    .groups = "drop"
  )

# Adjusted Plot 1: Fleiss Kappa vs BoostRounds by LearningRate
pALL_FLEISS <-ggplot(summary_df, aes(x = BoostRounds, y = mean_FLEISS, color = as.factor(LearningRate))) +
  geom_line(linewidth = 0.05, alpha = 0.5, linetype = "dashed") +  
  geom_errorbar(aes(
    ymin = mean_FLEISS - se_FLEISS,
    ymax = mean_FLEISS + se_FLEISS
  ), width = 100) +
  geom_smooth(data = FLEISS_data, 
              aes(x = as.numeric(BoostRounds), y = Value, color = as.factor(LearningRate), 
                  fill = as.factor(LearningRate)),
              method = "loess", level = 0.95, span = 1.0,
              linetype = "dotted", alpha = 0.1, linewidth = 0.5, show.legend = FALSE) +
  geom_point(data = FLEISS_data, 
             aes(x = BoostRounds, y = Value, color = as.factor(LearningRate)),
             position = position_jitter(width = 0.1), alpha = 0.8, size = 0.75) +
  # facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  labs(title = "FLEISS vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "FLEISS",
       color = "Learning Rate") +
  theme_minimal()

summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
p_all_fleiss <- ggplot(summary_df,
                       aes(x      = BoostRounds,
                           y      = mean_FLEISS,
                           colour = LearningRate,
                           shape  = LearningRate)) +
  geom_line(#position = position_dodge(width = 500),
            linewidth = 0.65,
            aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_FLEISS - se_FLEISS,
                    ymax = mean_FLEISS + se_FLEISS),
                width    = 100,
                ) + #position = position_dodge(width = 500)
  geom_point(size     = 1.25,
             stroke   = 0.5 ) +
             #position = position_dodge(width = 500)) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "Fleiss Kappa vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "Fleiss Kappa",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(p_all_fleiss, "FLEISS_ALL_vs_BoostRounds_cb", dpi_png = 600)

#### FLEISS (ALL TRAITS) normalized #### -------------------------------------------------

# 1. Normalise 0–1 across all Fleiss Kappa values
FLEISS_data_norm <- FLEISS_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE)))

# 2. Summarise normalised data by BoostRounds & LearningRate
summary_df_norm <- FLEISS_data_norm %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_FLEISS = mean(Value_norm),
    sd_FLEISS   = sd(Value_norm),
    se_FLEISS   = sd_FLEISS / sqrt(n()),
    min_FLEISS  = min(Value_norm),
    max_FLEISS  = max(Value_norm),
    .groups     = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)

# 4. Plot normalized Fleiss Kappa (0–1)
p_all_fleiss_norm <- ggplot(summary_df_norm,
                            aes(x      = BoostRounds,
                                y      = mean_FLEISS,
                                colour = LearningRate,
                                shape  = LearningRate)) +
  geom_line(linewidth = 0.65, aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_FLEISS - se_FLEISS,
                    ymax = mean_FLEISS + se_FLEISS),
                width = 100) +
  geom_point(size   = 1.25,
             stroke = 0.5) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "Fleiss Kappa (normalized 0–1) vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "Fleiss Kappa (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 5. Export the high-res normalized plot
export_hires_4_traits(p_all_fleiss_norm,
                      "Fleiss_ALL_vs_BoostRounds_cb_norm",
                      dpi_png = 600)


# Pearson

# Filter the data for ICC metric
pearson_data <- Master_df %>% 
  filter(Metric == "Pearson")

# Summarize the data by computing the mean, min, and max for each combination of BoostRounds, LearningRate, and Phenotype
# Summarize the data: calculating mean, standard error, and min/max
summary_df <- pearson_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_pearson = mean(Value),
    sd_pearson = sd(Value),
    se_pearson = sd(Value) / sqrt(n()),
    min_pearson = min(Value),
    max_pearson = max(Value),
    .groups = "drop"
  ) 

summary_df_subset <- pearson_data %>%
  group_by(BoostRounds, LearningRate, Phenotype, Subset) %>%
  summarize(
    mean_pearson_subset = mean(Value),
    sd_pearson_subset = sd(Value),
    se_pearson_subset = sd(Value) / sqrt(n()),
    min_pearson_subset = min(Value),
    max_pearson_subset = max(Value),
    .groups = "drop"
  )




# Plot the data with error bars and a trend line
#### PRESENT ME
pALL_PEARSON <-ggplot(summary_df, aes(x = BoostRounds, y = mean_pearson, color = as.factor(LearningRate))) +
  geom_line(position = position_dodge(width = 50)) +  # Line connecting the mean values
  geom_errorbar(aes(ymin = mean_pearson - se_pearson, ymax = mean_pearson + se_pearson), width = 500, position = position_dodge(width = 50)) +  # Error bars for the range of values
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line (line of best fit)
  
  labs(title = "Pearson Correlation vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "Pearson Correlation",
       color = "Learning Rate") + #,
  theme_minimal()


summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
p_all_pear <- ggplot(summary_df,
            aes(x      = BoostRounds,
                y      = mean_pearson,
                colour = LearningRate,      # colour by LR
                shape  = LearningRate)) +   # shape by LR
  geom_line(position = position_dodge(width = 500),
            linewidth = 0.65,
            aes(linetype = NULL)) +      # you could also map linetype
  geom_errorbar(aes(ymin = mean_pearson - se_pearson,
                    ymax = mean_pearson + se_pearson),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.25,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  # facet by phenotype with free y-scales
 # facet_wrap(~ Phenotype, scales = "free_y") +
  # your custom CB‐friendly scales
  
  pal_col_cb +
  pal_shape_cb +
  # unify colour + shape into one legend
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "Pearson Correlation vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "Pearson Correlation",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom")+ 
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )


# 4. Export at high res for 4-panel layout
export_hires_4_traits(p_all_pear, "Pearson_ALL_vs_BoostRounds_cb",
                      dpi_png   = 600)


#### Pearson (ALL TRAITS) normalized #### -------------------------------------------------

# 1. Normalise 0–1 across all Pearson values
pearson_data_norm <- pearson_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE)))

# 2. Summarise normalised data by BoostRounds & LearningRate
summary_df_norm <- pearson_data_norm %>%
  group_by(BoostRounds, LearningRate) %>%
  summarise(
    mean_pearson = mean(Value_norm),
    sd_pearson   = sd(Value_norm),
    se_pearson   = sd_pearson / sqrt(n()),
    min_pearson  = min(Value_norm),
    max_pearson  = max(Value_norm),
    .groups       = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate   <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)

# 4. Plot normalized Pearson Correlation (0–1)
p_all_pear_norm <- ggplot(summary_df_norm,
                          aes(x      = BoostRounds,
                              y      = mean_pearson,
                              colour = LearningRate,
                              shape  = LearningRate)) +
  geom_line(linewidth = 0.65, aes(linetype = NULL), position = position_dodge(width = 500)) +
  geom_errorbar(aes(ymin = mean_pearson - se_pearson,
                    ymax = mean_pearson + se_pearson),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.25,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "Pearson Correlation (normalized 0–1) vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "Pearson Correlation (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 5. Export the high-res normalized plot
export_hires_4_traits(p_all_pear_norm,
                      "Pearson_ALL_vs_BoostRounds_cb_norm",
                      dpi_png = 600)



# R SQUARED 

# Filter the data for R_squared metric
r2_data <- Master_df %>% 
  filter(Metric == "R_squared")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- r2_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_r2 = mean(Value),
    sd_r2 = sd(Value),
    se_r2 = sd(Value) / sqrt(n()),
    min_r2 = min(Value),
    max_r2 = max(Value),
    .groups = "drop"
  )

# Plot the data with error bars and a trend line
pALL_Rsquared <-ggplot(summary_df, aes(x = BoostRounds, y = mean_r2, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean R_squared values
  geom_errorbar(aes(ymin = mean_r2 - se_r2, ymax = mean_r2 + se_r2), width = 100) +  # Error bars for the range of values (optional)
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line (line of best fit)
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_r2, color = as.factor(LearningRate))) + 
  # facet_wrap(~Phenotype) +
  labs(title = "R-squared vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "R-squared",
       color = "Learning Rate") +
  theme_minimal()

# 1. Re‐level your LearningRate factor once
summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)

# 2. R-squared plot
p_all_r2 <- ggplot(summary_df,
                   aes(x      = BoostRounds,
                       y      = mean_r2,
                       colour = LearningRate,
                       shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500), # keep! 
            linewidth = 0.65,
            aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_r2 - se_r2,
                    ymax = mean_r2 + se_r2),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.25,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "R-squared vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "R-squared",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(p_all_r2,   "R2_ALL_vs_BoostRounds_cb",   dpi_png = 600)

#### R Squared (ALL TRAITS) normalized #### -------------------------------------------------

# 1. Normalise 0–1 across all R² values
r2_data_norm <- r2_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE)))

# 2. Summarise normalised data by BoostRounds & LearningRate
summary_df_norm <- r2_data_norm %>%
  group_by(BoostRounds, LearningRate) %>%
  summarise(
    mean_r2 = mean(Value_norm),
    sd_r2   = sd(Value_norm),
    se_r2   = sd_r2 / sqrt(n()),
    min_r2  = min(Value_norm),
    max_r2  = max(Value_norm),
    .groups = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)

# 4. Plot normalized R² (0–1)
p_all_r2_norm <- ggplot(summary_df_norm,
                        aes(x      = BoostRounds,
                            y      = mean_r2,
                            colour = LearningRate,
                            shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500),
            linewidth = 0.65,
            aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_r2 - se_r2,
                    ymax = mean_r2 + se_r2),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.25,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "R² (normalized 0–1) vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "R² (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 5. Export the high-res normalized plot
export_hires_4_traits(p_all_r2_norm,
                      "R2_ALL_vs_BoostRounds_cb_norm",
                      dpi_png = 600)


# AUC
# Filter the data for AUC metric
auc_data <- Master_df %>% 
  filter(Metric == "AUC")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- auc_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_auc = mean(Value),
    sd_auc = sd(Value),
    se_auc = sd(Value) / sqrt(n()),
    min_auc = min(Value),
    max_auc = max(Value),
    .groups = "drop"
  )

# Plot the data with error bars and a trend line (Plot 1)
pALL_AUC <-ggplot(summary_df, aes(x = BoostRounds, y = mean_auc, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean AUC values
  geom_errorbar(aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc), width = 100) +  # Error bars for the range of values (optional)
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_auc, color = as.factor(LearningRate))) +
  # facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  labs(title = "AUC vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "AUC",
       color = "Learning Rate") +
  theme_minimal()

summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
p_all_auc <- ggplot(summary_df,
                    aes(x      = BoostRounds,
                        y      = mean_auc,
                        colour = LearningRate,
                        shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500),
            linewidth = 0.65,
            aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_auc - se_auc,
                    ymax = mean_auc + se_auc),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.25,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "AUC vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "AUC",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(p_all_auc,  "AUC_ALL_vs_BoostRounds_cb",  dpi_png = 600)

#### AUC (ALL TRAITS) normalized #### -------------------------------------------------

# 1. Normalise 0–1 across all AUC values
auc_data_norm <- auc_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE)))

# 2. Summarise normalised data by BoostRounds & LearningRate
summary_df_norm <- auc_data_norm %>%
  group_by(BoostRounds, LearningRate) %>%
  summarise(
    mean_auc = mean(Value_norm),
    sd_auc   = sd(Value_norm),
    se_auc   = sd_auc / sqrt(n()),
    min_auc  = min(Value_norm),
    max_auc  = max(Value_norm),
    .groups  = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)

# 4. Plot normalized AUC (0–1)
p_all_auc_norm <- ggplot(summary_df_norm,
                         aes(x      = BoostRounds,
                             y      = mean_auc,
                             colour = LearningRate,
                             shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500),
            linewidth = 0.65,
            aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_auc - se_auc,
                    ymax = mean_auc + se_auc),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.25,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "AUC (normalized 0–1) vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "AUC (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 5. Export the high-res normalized plot
export_hires_4_traits(p_all_auc_norm,
                      "AUC_ALL_vs_BoostRounds_cb_norm",
                      dpi_png = 600)


# NDCG at 20 % 

# Filter the data for NDCG_at_20_percent metric
ndcg_data <- Master_df %>% 
  filter(Metric == "NDCG_at_20_percent")

# Summarize the data by computing the mean, standard error, and min/max for each combination of BoostRounds, LearningRate, and Phenotype
summary_df <- ndcg_data %>%
  group_by(BoostRounds, LearningRate) %>%
  summarize(
    mean_ndcg = mean(Value),
    sd_ndcg = sd(Value),
    se_ndcg = sd(Value) / sqrt(n()),
    min_ndcg = min(Value),
    max_ndcg = max(Value),
    .groups = "drop"
  )

# Plot the data with error bars and a trend line (Plot 1)
pALL_NDCG20 <-ggplot(summary_df, aes(x = BoostRounds, y = mean_ndcg, color = as.factor(LearningRate))) +
  geom_line() +  # Line connecting the mean NDCG values
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg, ymax = mean_ndcg + se_ndcg ), width = 100) +  # Error bars for the range (optional)
  geom_smooth(method = "loess", se = FALSE, linetype = "dotted", alpha = 0.1, linewidth = 0.5) +  # Overall trend line
  geom_point(data = summary_df, aes(x = BoostRounds, y = mean_ndcg, color = as.factor(LearningRate))) +
  #facet_wrap(summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("GY_1", "TGW_1", "GPC_1"))) + # "HET_1"
  labs(title = "NDCG@20% vs BoostRounds by LearningRate",
       x = "Boost Rounds",
       y = "NDCG@20%",
       color = "Learning Rate") +
  theme_minimal()

summary_df$LearningRate <- factor(as.character(summary_df$LearningRate), levels = lr_levels)
p_all_ndcg <- ggplot(summary_df,
                     aes(x      = BoostRounds,
                         y      = mean_ndcg,
                         colour = LearningRate,
                         shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500),
            linewidth = 0.65,
            aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg,
                    ymax = mean_ndcg + se_ndcg),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.25,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "NDCG@20% vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "NDCG@20%",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

export_hires_4_traits(p_all_ndcg, "NDCG_ALL_vs_BoostRounds_cb", dpi_png = 600)

#### NDCG@20% (ALL TRAITS) normalized #### -------------------------------------------------

# 1. Normalise 0–1 across all NDCG@20% values
ndcg_data_norm <- ndcg_data %>%
  group_by(Phenotype) %>%
  mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE)))

# 2. Summarise normalised data by BoostRounds & LearningRate
summary_df_norm <- ndcg_data_norm %>%
  group_by(BoostRounds, LearningRate) %>%
  summarise(
    mean_ndcg = mean(Value_norm),
    sd_ndcg   = sd(Value_norm),
    se_ndcg   = sd_ndcg / sqrt(n()),
    min_ndcg  = min(Value_norm),
    max_ndcg  = max(Value_norm),
    .groups   = "drop"
  )

# 3. Keep factor levels consistent
summary_df_norm$LearningRate <- factor(as.character(summary_df_norm$LearningRate), levels = lr_levels)

# 4. Plot normalized NDCG@20% (0–1)
p_all_ndcg_norm <- ggplot(summary_df_norm,
                          aes(x      = BoostRounds,
                              y      = mean_ndcg,
                              colour = LearningRate,
                              shape  = LearningRate)) +
  geom_line(position = position_dodge(width = 500),
            linewidth = 0.65,
            aes(linetype = NULL)) +
  geom_errorbar(aes(ymin = mean_ndcg - se_ndcg,
                    ymax = mean_ndcg + se_ndcg),
                width    = 500,
                position = position_dodge(width = 500)) +
  geom_point(size     = 1.25,
             stroke   = 0.5,
             position = position_dodge(width = 500)) +
  pal_col_cb +
  pal_shape_cb +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE),
         shape  = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(#title  = "NDCG@20% (normalized 0–1) vs BoostRounds by LearningRate",
       x      = "Boost Rounds",
       y      = "NDCG@20% (normalized 0–1)",
       colour = "Learning Rate",
       shape  = "Learning Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = c(100, 1080, 2060, 3040, 4020, 5000),
    labels = c("100", "1080", "2060", "3040", "4020", "5000")
  )

# 5. Export the high-res normalized plot
export_hires_4_traits(p_all_ndcg_norm,
                      "NDCG20_ALL_vs_BoostRounds_cb_norm",
                      dpi_png = 600)


# Plot all together 
# pALL_ICC + pALL_FLEISS + pALL_PEARSON + pALL_Rsquared + pALL_AUC + pALL_NDCG20
library(gridExtra) 
grid.arrange(pALL_ICC, pALL_FLEISS, pALL_PEARSON, pALL_Rsquared, pALL_AUC, pALL_NDCG20, ncol = 2)
# Save the plots to a PDF file
ggsave("TWR_Definitive_GY_Super_Dataset_Metrics_Plots.pdf", 
       arrangeGrob(pALL_ICC, pALL_FLEISS, pALL_PEARSON, pALL_Rsquared, pALL_AUC, pALL_NDCG20, ncol = 2), 
       width = 12, height = 10)




# Updated: all together
library(gridExtra) 
grid.arrange(p_all_icc, p_all_fleiss, p_all_pear, p_all_r2, p_all_auc, p_all_ndcg, ncol = 2)
# Save the plots to a PDF file
ggsave("TWR_Definitive_GY_Super_Dataset_Metrics_Plots.pdf", 
       arrangeGrob(p_all_icc, p_all_fleiss, p_all_pear, p_all_r2, p_all_auc, p_all_ndcg, ncol = 2), 
       width = 12, height = 10)

#  Export the high-res normalized plot
export_hires_4_traits(grid.arrange(p_all_icc, p_all_fleiss, p_all_pear, p_all_r2, p_all_auc, p_all_ndcg, ncol = 2),
                      "ALL_traits_agg", 
                      width_in  = 12,
                      height_in = 15,
                      dpi_png   = 600)


# output normalized plots together: 

# Updated: all together
library(gridExtra) 
grid.arrange(p_all_icc_norm, p_all_fleiss_norm, p_all_pear_norm, p_all_r2_norm, p_all_auc_norm, p_all_ndcg_norm, ncol = 2)
# Save the plots to a PDF file
ggsave("TWR_Definitive_GY_Super_Dataset_Metrics_Plots_NORMALIZED.pdf", 
       arrangeGrob(p_all_icc_norm, p_all_fleiss_norm, p_all_pear_norm, p_all_r2_norm, p_all_auc_norm, p_all_ndcg_norm, ncol = 2), 
       width = 12, height = 10)

  
 
  # 5. Export the high-res normalized plot
  export_hires_4_traits(grid.arrange(p_all_icc_norm, p_all_fleiss_norm, p_all_pear_norm, p_all_r2_norm, p_all_auc_norm, p_all_ndcg_norm, ncol = 2),
                        "ALL_traits_agg_normalized", 
                        width_in  = 12,
                        height_in = 15,
                        dpi_png   = 600)
                  
  
  
  
  
  ### labelled agg plot master: 
  
  
  # Load necessary libraries first
  library(gridExtra)
  library(ggplot2)
  library(grid) # For grid.draw()
  
  # 1. Manually add a centered title to each plot object, creating new variables
  p_all_icc_labelled_a    <- p_all_icc    + labs(title = "(a)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_fleiss_labelled_b <- p_all_fleiss + labs(title = "(b)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_pear_labelled_c   <- p_all_pear   + labs(title = "(c)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_r2_labelled_d     <- p_all_r2     + labs(title = "(d)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_auc_labelled_e    <- p_all_auc    + labs(title = "(e)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_ndcg_labelled_f   <- p_all_ndcg   + labs(title = "(f)") + theme(plot.title = element_text(hjust = 0.5))
  
 
  
  export_hires_4_traits(grid.arrange( p_all_icc_labelled_a, p_all_fleiss_labelled_b, p_all_pear_labelled_c, 
                                      p_all_r2_labelled_d, p_all_auc_labelled_e, p_all_ndcg_labelled_f, 
                                      ncol = 2),
                        "ALL_traits_agg_Labelled", 
                        width_in  = 12,
                        height_in = 15,
                        dpi_png   = 600)
  
  
  # 1. Manually add a centered title to each NORMALIZED plot object, creating new variables
  p_all_icc_norm_labelled_a    <- p_all_icc_norm    + labs(title = "(a)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_fleiss_norm_labelled_b <- p_all_fleiss_norm + labs(title = "(b)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_pear_norm_labelled_c   <- p_all_pear_norm   + labs(title = "(c)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_r2_norm_labelled_d     <- p_all_r2_norm     + labs(title = "(d)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_auc_norm_labelled_e    <- p_all_auc_norm    + labs(title = "(e)") + theme(plot.title = element_text(hjust = 0.5))
  p_all_ndcg_norm_labelled_f   <- p_all_ndcg_norm   + labs(title = "(f)") + theme(plot.title = element_text(hjust = 0.5))
  
  # 2. Arrange the newly labeled plots and export them using your custom function
  export_hires_4_traits(
    grid.arrange(
      p_all_icc_norm_labelled_a, p_all_fleiss_norm_labelled_b, p_all_pear_norm_labelled_c, 
      p_all_r2_norm_labelled_d, p_all_auc_norm_labelled_e, p_all_ndcg_norm_labelled_f, 
      ncol = 2
    ),
    basename = "ALL_traits_agg_normalized_Labelled", 
    width_in  = 12,
    height_in = 15,
    dpi_png   = 600
  )
  