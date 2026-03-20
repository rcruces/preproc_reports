# ----------------------------------------------------------------------
# MICs mk6240 subcortical database
# January, 2026
# 
# ----------------------------------------------------------------------

# Directory
dir_path <- "/Users/rcruces/Desktop/"

# Load the database
subcortical <- read.csv("/Users/rcruces/git_here/MICs_datasets/functions/data/subcortical_trc-mk6240.csv")
mk6240 <- read.csv("/Users/rcruces/Desktop/mk6240.csv")
mk.df <- read.csv("/Users/rcruces/Desktop/OSF_data_revision-2025/18F-MK6240_database.csv")
mk6240.2 <- read.csv("/Users/rcruces/git_here/MICs_datasets/functions/data/mk6240.csv")

# Create id column
subcortical$id <- gsub("ses-", "", gsub("sub-","",subcortical$subject_id))

# Change column type
mk6240.2$group <- factor(mk6240.2$group, levels = c("Patient", "Healthy"))
mk6240$group <- factor(mk6240$group, levels = c("Patient", "Healthy"))
mk.df$group <- factor(mk.df$group, levels = c("Patient", "Healthy"))

# Merge databases (include lateralization, matched.25)
subcortical <- merge(subcortical, mk6240[,c("id","sex","age.mk6240","group","lateralization","matched.25")], by="id")

# Slice to keep only matched.25 = 1
subcortical <- subcortical[subcortical$matched.25 == 1,]
mk6240 <- mk6240[mk6240$matched.25 == 1,]
mk6240.2 <- mk6240.2[mk6240.2$matched.25 == 1,]


# ----------------------------------------------------------------------
# Data Processing
# Boolean with the healthy controls
hc <- which(subcortical$group == "Healthy")

#### Z-scores two vectors based on the first  ####
Zscore <- function(vec1,vec2) {
  vec1 <- as.vector(vec1)
  vec2 <- as.vector(vec2)
  mu <- mean(vec1, na.rm=TRUE)
  s <- sd(vec1, na.rm=TRUE)
  Zvec <- (vec2 - mu)/s
  return(Zvec)
}

# Get the SUVR with cerebellar reference
pet.structures <- c("pet.L.Thalamus","pet.L.Caudate","pet.L.Putamen","pet.L.Pallidus","pet.L.Hippocampus","pet.L.Amygdala","pet.L.Accumbens",
                    "pet.R.Thalamus","pet.R.Caudate","pet.R.Putamen","pet.R.Pallidus","pet.R.Hippocampus","pet.R.Amygdala","pet.R.Accumbens")
subcortical[,pet.structures] <- subcortical[,pet.structures]/subcortical$cerebellarGM

# ----------------------------------------------------------------------
# Test similarity of values (all LEFT and HC must be the same)
mk6240.new <- subcortical[which(subcortical$lateralization == "L"), "pet.L.Hippocampus"]
mk6240.old <- mk6240[which(mk6240$lateralization == "L"), "suvr.ipsi.hippocampus"]
all.equal(mk6240.old, mk6240.new)

# ----------------------------------------------------------------------
# Make a copy of values
subcortical.flipped <- subcortical

# ----------------------------------------------------------------------
# FIRST Zscore the values
subcortical[,pet.structures] <- apply(subcortical[,pet.structures], 2, function(x) Zscore(x[hc], x))

# ----------------------------------------------------------------------
# Flip Right to Left function
flip_RtoL <- function(subcortical, pet.structures) {
  
  # Flip Right to left
  latR <- !is.na(subcortical$lateralization) & subcortical$lateralization == "R"
  
  # Extract unique region names (remove "pet.L." / "pet.R.")
  regions <- unique(substr(pet.structures, 7, 999))
  
  # Loop to create ipsi / contralateral columns
  for (reg in regions) {
    
    Lcol <- paste0("pet.L.", reg)
    Rcol <- paste0("pet.R.", reg)
    
    # ipsilateral SUVR
    subcortical[[paste0("suvr.ipsi.", tolower(reg))]] <-
      ifelse(latR, subcortical[[Rcol]], subcortical[[Lcol]])
    
    # contralateral SUVR
    subcortical[[paste0("suvr.cntr.", tolower(reg))]] <-
      ifelse(latR, subcortical[[Lcol]], subcortical[[Rcol]])
  }
  
  return(subcortical)
}

# zscored then flipped
subcortical <- flip_RtoL(subcortical, pet.structures)

# only flipped
subcortical.flipped <- flip_RtoL(subcortical.flipped, pet.structures)

# ----------------------------------------------------------------------
# T-test Ipsilateral hippocampus 
t.test(suvr.ipsi.hippocampus ~ group, data=subcortical)

t.test(suvr.ipsi.hippocampus ~ group, data=subcortical.flipped)

t.test(suvr.ipsi.hippocampus ~ group, data=mk6240)

t.test(suvr.ipsi.hippocampus ~ group, data=mk.df)

t.test(pet.L.Hippocampus ~ group, data=mk6240.2)

# ----------------------------------------------------------------------
# Tables 
# Define subcortical regions
sctx.str <- c("suvr.ipsi.thalamus", "suvr.ipsi.caudate", "suvr.ipsi.putamen", "suvr.ipsi.pallidus","suvr.ipsi.amygdala", "suvr.ipsi.hippocampus", "suvr.ipsi.accumbens", 
              "suvr.cntr.thalamus", "suvr.cntr.caudate", "suvr.cntr.putamen", "suvr.cntr.pallidus","suvr.cntr.amygdala", "suvr.cntr.hippocampus", "suvr.cntr.accumbens")

library(dplyr)
library(gtsummary)
library(knitr)
library(kableExtra)

# Mean and Std
subcortical %>%
  dplyr::select(
    group, all_of(sctx.str)) %>%
  tbl_summary(
    by = group,
    missing = "no",
    statistic = list(all_continuous() ~ "{mean}±{sd}"),
  ) %>%                                   # <<< closed tbl_summary() properly
  modify_header(label = "**zscore > flipping**") %>%
  add_p(test = all_continuous() ~ "t.test") %>%
  modify_header(statistic ~ "**Statistic**") %>%
  as_kable_extra(
    booktabs = TRUE,
    longtable = TRUE,
    linesep = ""
  ) %>%
  kableExtra::kable_styling(
    position = "left",
    latex_options = c("striped", "repeat_header"),
    stripe_color = "gray!15"
  )

# Mean and Std
subcortical.flipped %>%
  dplyr::select(
    group, all_of(sctx.str)) %>%
  tbl_summary(
    by = group,
    missing = "no",
    statistic = list(all_continuous() ~ "{mean}±{sd}"),
  ) %>%                                   # <<< closed tbl_summary() properly
  modify_header(label = "**Only Flipping**") %>%
  add_p(test = all_continuous() ~ "t.test") %>%
  modify_header(statistic ~ "**Statistic**") %>%
  as_kable_extra(
    booktabs = TRUE,
    longtable = TRUE,
    linesep = ""
  ) %>%
  kableExtra::kable_styling(
    position = "left",
    latex_options = c("striped", "repeat_header"),
    stripe_color = "gray!15"
  )


# Mean and Std
mk6240 %>%
  dplyr::select(
    group, all_of(sctx.str)) %>%
  tbl_summary(
    by = group,
    missing = "no",
    statistic = list(all_continuous() ~ "{mean}±{sd}"),
  ) %>%                                   # <<< closed tbl_summary() properly
  modify_header(label = "**Old Only Flipping**") %>%
  add_p(test = all_continuous() ~ "t.test") %>%
  modify_header(statistic ~ "**Statistic**") %>%
  as_kable_extra(
    booktabs = TRUE,
    longtable = TRUE,
    linesep = ""
  ) %>%
  kableExtra::kable_styling(
    position = "left",
    latex_options = c("striped", "repeat_header"),
    stripe_color = "gray!15"
  )


mk.df %>%
  dplyr::select(
    group, all_of(sctx.str)) %>%
  tbl_summary(
    by = group,
    missing = "no",
    statistic = list(all_continuous() ~ "{mean}±{sd}"),
  ) %>%                                   # <<< closed tbl_summary() properly
  modify_header(label = "**review flipping**") %>%
  add_p(test = all_continuous() ~ "t.test") %>%
  modify_header(statistic ~ "**Statistic**") %>%
  as_kable_extra(
    booktabs = TRUE,
    longtable = TRUE,
    linesep = ""
  ) %>%
  kableExtra::kable_styling(
    position = "left",
    latex_options = c("striped", "repeat_header"),
    stripe_color = "gray!15"
  )

