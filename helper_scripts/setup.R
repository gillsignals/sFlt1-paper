#### Load packages
library(tidyverse)
library(patchwork)
library(R.matlab)
library(viridis)
library(ggrepel)
library(GGally)
library(latex2exp)
library(testthat)
library(ComplexHeatmap)
library(circlize)
library(scales)

#### Set options

# specify decimal resolution
options(digits = 4)

#### Define helper functions

# function to load a colorblind friendly palette
load_cb_pal <- function(){
  
  # define colorblind-friendly palette
  # black, red, sky blue, orange, green, yellow, dark blue, pink
  cb_pal <- c("#000000", "#D55E00",  "#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2",  "#CC79A7")
  
  # return color palette
  cb_pal
  
}

# define purple -> green color palette
# endpoints taken from RColorBrewer "PRGn"
pg_color_func <- colorRamp2(c(-1, 0, 1), c("#762A83", "white", "#1B7837"))

# define purple -> green color palette with more shades
# endpoints taken from RColorBrewer "PRGn"
pg_color_func_broad <- colorRamp2(c(-2.5, 0, 2.5), c("#762A83", "white", "#1B7837"))

# function to set basic plot theme
set_theme1 <- function(){
  
  # define theme basics
  theme1 <- theme_light() +
    theme(
      text = element_text(size = 16, family = "serif"),
      plot.title = element_text(hjust = 0.5),
      strip.background =element_rect(fill="gray50"),
      strip.text = element_text(color = "white", face = "bold")
    )
  
  # return theme
  theme1
  
}

#### Apply helper functions

# set theme for all plots
theme1 <- set_theme1()

# define colorblind-friendly palette
# black, red, sky blue, orange, green, yellow, dark blue, pink
cb_pal <- load_cb_pal()

# Add "big" theme
theme_big <- theme1 +
  theme(legend.position = "bottom",
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10))

# Add "smallfont" theme
theme_smallfont <- theme1 +
  theme(legend.position = "bottom",
        axis.title.y = element_text(angle = 90, vjust = 0.5, size = 16),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10))

# Add "poster" theme
theme_poster <- theme1 +
  theme(text = element_text(size = 24, family = "serif"),
        axis.text = element_text(size = 18),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 24),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom")

# Add "paper" theme
theme_paper <- theme_poster +
  theme(axis.title.y = element_text(angle=90, margin = margin(r = 10, l = 10)),
        plot.margin = margin(r = 20, t = 10))

# function to create a subdirectory in the current output directory: outfile_prefix/subdir_name/
create_subdir <- function(subdir_name){
  
  # if subdirectory doesn't exist, make it
  if(!(file.exists(paste0(outfile_prefix, subdir_name)))){
    setwd(outfile_prefix)
    dir.create(subdir_name, showWarnings = FALSE)
  }
  
  # return path to desired subdirectory
  paste0(outfile_prefix, subdir_name, "/")
}

# convert string filename to "outfile_prefix/filename.png"
out_png <- function(filename){
  paste0(outfile_prefix, filename, ".png")
}

#### Map parameter names to colors and define resulting color palette (alpha:tau)

# map parameter name to index of desired color
param <- list(alpha = 3, beta = 2, gamma = 4, delta = 5, tau = 8)

# store parameter names
param_names <- names(param)

# create palette for parameters
param_pal <- sapply(1:length(param), function(x){
  cb_pal[unlist(param[param_names[x]])]
})

# create parameter palette (list to map parameter names to values)
param_pal_list <- list()
for (i in 1:length(param_names)) {
  param_pal_list[param_names[i]] <- cb_pal[unlist(param[param_names[i]])]
}

# set labels for Greek letters in plots
color_breaks <- names(param_pal_list)
color_labels <- parse(text = color_breaks)

#### Specify geometry and conversion factor

num_cells <- 5000   # number of cells
Vx <- 3e-4          # extracellular volume (L/well)

# Specify conversion factor: (#/cell -> pmol/well)
# (#/cell * cells/well * pmol/# = pmol/well)
Nav_pmol <- 6.022e11
conv_factor_pmol <- num_cells / Nav_pmol
    
# specify conversion factor: (pmol/well to ng/mL)
# pmol/well * mol/pmol * L/mL * g/mol * ng/g = ng/well
# ng/mL = pmol/well * mol/pmol * L/mL * g/mol * ng/g * well/L
#pmol_2_ngml = 1e-12 * 1e-3 * 1.1e5 * 1e9 / Vx;      % literature MW sFlt1: 110 kDa
pmol_2_ngml <- 1e-12 * 1e-3 * 9e4 * 1e9 / Vx     # KK manuscript MW sFlt1: 90 kDa
    
# specify conversion factor: (#/cell -> ng/mL)
conv_factor_numcell_to_ngml <- conv_factor_pmol * pmol_2_ngml

#### Load experimental data

# load processed experimental data from Jung (pulse-chase, relative units)
load("../../saved-data/jung_1d.rda")

jung <- jung %>%
  filter(!(nominal_time == 10 & target == "media")) %>%
  select(target, nominal_time, sFlt) %>%
  pivot_wider(names_from = "target", values_from = sFlt)

names(jung) <- c("time", "lysate_fc", "media_frac_max")

# load control data from Kinghorn (constitutive, relative units)
kinghorn <- readMat("../../saved-data/kinghorn_ctrl_summary.mat")

kinghorn <- as.data.frame(kinghorn$kinghorn[,,1])

names(kinghorn) <- str_replace_all(names(kinghorn), "\\.", "_")

# load experimental data from Hornig (constitutive, absolute units)
load("../../saved-data/hornig.rda")

names(hornig) <- str_replace(names(hornig), "Sx", "X")

hornig <- hornig %>%
  filter(condition == "basal", nominal_time != 3)

hornig_X_ngml_24 <- hornig %>%
  filter(nominal_time == 24) %>%
  pull(X_ng_ml)

hornig_X_numcell_24 <- hornig_X_ngml_24 / conv_factor_numcell_to_ngml

hornig <- hornig %>%
  mutate(X_norm_24 = X_ng_ml / hornig_X_ngml_24,  # normalize to 24h amount
         X_numcell = X_ng_ml / conv_factor_numcell_to_ngml) # convert units from ng/ml to #/cell)






