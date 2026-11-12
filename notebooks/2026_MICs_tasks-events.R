
# =========================================================================
# MICs dataset
# Functional memory task accuracy and bias-corrected scores
# =========================================================================
#
# OVERVIEW
#
# This script extracts behavioral performance measures from the functional
# task events.tsv files in the MICs BIDS dataset.
#
# The output contains one row per participant and session, including:
#
#   - Mean accuracy for four functional memory tasks:
#       MST2, Episodic, Semantic, and Spatial.
#   - Bias-corrected pattern separation (BPS) and pattern completion (BPC)
#       scores derived from the MST2 task.
#
# Task accuracy scores range from 0 to 1, where higher values indicate
# better performance. BPS and BPC are difference scores and can range
# from -1 to 1.
#
# -------------------------------------------------------------------------
# TASK DESCRIPTIONS
#
# MST2: Mnemonic Similarity / Discrimination Task
#
# A recognition memory task designed to assess participants' ability to
# distinguish previously encountered objects from similar objects and
# novel objects.
#
# Performance is reported as mean accuracy.
# Additionally, bias-corrected pattern separation (BPS) and pattern
# completion (BPC) scores are calculated from the MST2 responses.
#
# Reference:
# Li, Q., Tavakol, S., Royer, J., Larivière, S., Vos De Wael, R.,
# Park, B. Y., ... & Bernhardt, B. C. (2021).
# Atypical neural topographies underpin dysfunctional pattern separation
# in temporal lobe epilepsy.
# Brain, 144(8), 2486-2498.
#
# Stark, S. M., Yassa, M. A., Lacy, J. W., & Stark, C. E. (2013).
# A task to assess behavioral pattern separation (BPS) in humans:
# Data from healthy aging and mild cognitive impairment.
# Neuropsychologia, 51(12), 2442-2449.
#
#
# Episodic: Episodic Retrieval Task
#
# Assesses episodic memory retrieval. Participants are asked to identify
# pairs of objects that were previously presented during the encoding phase.
#
# Performance is reported as mean accuracy (0-1).
#
# Reference:
# Cabalo, D. G., DeKraker, J., Royer, J., Xie, K., Tavakol, S.,
# Rodríguez-Cruces, R., ... & Bernhardt, B. C. (2024).
# Differential reorganization of episodic and semantic memory systems
# in epilepsy-related mesiotemporal pathology.
# Brain, 147(11), 3918-3932.
#
#
# Semantic: Semantic Retrieval Task
#
# Assesses semantic memory retrieval. Participants are asked to identify
# the object that is most conceptually related to a given prime object.
#
# Performance is reported as mean accuracy (0-1).
#
# Reference:
# Cabalo, D. G., DeKraker, J., Royer, J., Xie, K., Tavakol, S.,
# Rodríguez-Cruces, R., ... & Bernhardt, B. C. (2024).
# Differential reorganization of episodic and semantic memory systems
# in epilepsy-related mesiotemporal pathology.
# Brain, 147(11), 3918-3932.
#
#
# Spatial: Spatial Discrimination Task
#
# Assesses memory for spatial configurations. Participants are presented
# with three objects, and the spatial configuration of these objects is
# manipulated (e.g., changes in the distances between objects or rotations).
#
# Performance is reported as mean accuracy (0-1).
#
# Reference:
# Tavakol, S., Li, Q., Royer, J., Vos de Wael, R., Larivière, S., Lowe, A.,
# ... & Bernhardt, B. (2021).
# A structure-function substrate of memory for spatial configurations
# in medial and lateral temporal cortices.
# Cerebral Cortex, 31(7), 3213-3225.
#
#
# BIAS-CORRECTED MST2 SCORES
# -------------------------------------------------------------------------
#
# Both BPS and BPC are calculated using the MST2 trial responses.
# These scores are difference measures that account for responses to
# novel (new) trials.
#
# BPS: Bias-Corrected Pattern Separation
# BPS = P(similar | similar) - P(similar | new)
#
# This score measures the difference in the probability of responding
# "similar" to similar versus new trials.
#
# Interpretation:
#   Positive BPS: "similar" responses are more frequent for similar
#                 trials than for new trials, indicating discrimination.
#
#   BPS = 0:      "similar" responses occur at the same rate for similar
#                 and new trials, indicating no discrimination based
#                 on this response.
#
#   Negative BPS: "similar" responses are more frequent for new trials
#                 than for similar trials, indicating a reversed
#                 response pattern.
#
#
# BPC: Bias-Corrected Pattern Completion
# BPC = P(old | similar) - P(old | new)
#
# This score measures the difference in the probability of responding
# "old" to similar versus new trials.
#
# Interpretation:
#   Positive BPC: "old" responses are more frequent for similar trials
#                 than for new trials.
#
#   BPC = 0:      "old" responses occur at the same rate for similar
#                 and new trials.
#
#   Negative BPC: "old" responses are more frequent for new trials
#                 than for similar trials, indicating a reversed
#                 response pattern.
#
# Note:
# BPS and BPC are difference scores, not overall task accuracy measures.
# Negative values should be retained and should not automatically be
# interpreted as impaired memory performance.
#
# -------------------------------------------------------------------------
# OUTPUT VARIABLES
#
# participant_id : Participant identifier (without the "sub-" prefix).
# session_id     : Session identifier (without the "ses-0" prefix).
# task.mst2      : Mean MST2 task accuracy.
# task.episodic  : Mean episodic retrieval task accuracy.
# task.semantic  : Mean semantic retrieval task accuracy.
# task.spatial   : Mean spatial discrimination task accuracy.
# task.BPS       : Bias-corrected pattern separation score.
# task.BPC       : Bias-corrected pattern completion score.
# ID             : Combined participant-session identifier.
#
# Missing task files or unavailable scores are represented as NA.
#
# =========================================================================

# ----------------------------------------------------------------------
# MICs dataset
# Functional task accuracy and group table generation
# ----------------------------------------------------------------------
# Working directory
rawdata <- "/data_/mica3/BIDS_MICs/rawdata"
setwd(rawdata)

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
####  Safe mean: return numeric NA if no valid observations  #### 
safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {return(NA_real_)}
  mean(x, na.rm = TRUE) }

# Read TSV file
read_events <- function(task, subid, sesid) {
  func <- file.path(rawdata, subid, sesid, "func")
  filename <- paste0(subid, "_", sesid, "_task-", task, "_events.tsv")
  filepath <- file.path(func, filename)
  if (!file.exists(filepath)) { message("[MISSING] ", subid, " ", sesid, " ", task) 
    return(NULL)}
  read.delim(filepath, sep = "\t")
}

# ----------------------------------------------------------------------
####  Calculate mean accuracy  #### 
get_accu <- function(data, task) {
  if (is.null(data)) return(NA_real_)
  if (task == "MST2") {
    if (!all(c("trial_type", "subject.response") %in% names(data))) {
      return(NA_real_)
    }
    accu <- safe_mean(as.character(data$trial_type) == as.character(data$subject.response))
  } else {
    if (!"accuracy" %in% names(data)) {return(NA_real_)}
    accu <- safe_mean(as.numeric(as.character(data$accuracy)))
  }
  return(accu)
}

# ----------------------------------------------------------------------
####  Calculate Bias Corrected Pattern Separation and Completion  #### 
get_bcp <- function(data) {
  if (is.null(data)) { return(c(BPS = NA_real_, BPC = NA_real_))}
  if (!all(c("trial_type", "subject.response") %in% names(data))) {
    return(c(BPS = NA_real_, BPC = NA_real_))}
  
  # BPS = pattern separation scores
  p.sim <- which(data$trial_type == "similar")
  p.new <- which(data$trial_type == "new")
  
  # BPS = P(similar|similar) - P(similar|new)
  BPS <- safe_mean(
    as.character(data$subject.response[p.sim]) ==
      as.character(data$trial_type[p.sim])) -
    safe_mean(data$subject.response[p.new] == "similar")
  
  # BPC = P(old|similar) - P(old|new)
  BPC <- safe_mean(data$subject.response[p.sim] == "old") -
    safe_mean(data$subject.response[p.new] == "old")
  bcp <- c(BPS = BPS, BPC = BPC)
  return(bcp)
}


# ----------------------------------------------------------------------
# Get subjects and sessions
subs <- list.dirs(
  rawdata,
  recursive = FALSE,
  full.names = FALSE
)
subs <- grep("^sub-", subs, value = TRUE)

# ----------------------------------------------------------------------
# Extract task data
results <- list()
for (i in seq_along(subs)) {
  subid <- subs[i]
  ses <- list.dirs(
    file.path(rawdata, subid),
    recursive = FALSE,
    full.names = FALSE
  )
  ses <- grep("^ses-", ses, value = TRUE)
  
  for (j in seq_along(ses)) {
    sesid <- ses[j]
    
    # Read each task TSV only once
    data.mst2      <- read_events("MST2", subid, sesid)
    data.retrieval <- read_events("retrieval", subid, sesid)
    data.semantic  <- read_events("semantic", subid, sesid)
    data.spatial   <- read_events("spatial", subid, sesid)
    
    # Calculate accuracy
    accu.mst2      <- get_accu(data.mst2, "MST2")
    accu.retrieval <- get_accu(data.retrieval, "retrieval")
    accu.semantic  <- get_accu(data.semantic, "semantic")
    accu.spatial   <- get_accu(data.spatial, "spatial")
    
    # Calculate BPS and BPC using the already loaded MST2 data
    bcp <- get_bcp(data.mst2)
    
    # Store results (no repeated rbind)
    results[[length(results) + 1L]] <- data.frame(
      participant_id = as.character(subid),
      session_id = as.character(sesid),
      task.mst2 = accu.mst2,
      task.episodic = accu.retrieval,
      task.semantic = accu.semantic,
      task.spatial = accu.spatial,
      task.BPS = unname(bcp["BPS"]),
      task.BPC = unname(bcp["BPC"]),
      stringsAsFactors = FALSE
    )
  }
}


# ----------------------------------------------------------------------
# Combine results
if (length(results) > 0) {
  task.data <- do.call(rbind, results)
  rownames(task.data) <- NULL
} else {
  task.data <- data.frame(
    participant_id = character(),
    session_id = character(),
    task.mst2 = numeric(),
    task.episodic = numeric(),
    task.semantic = numeric(),
    task.spatial = numeric(),
    task.BPS = numeric(),
    task.BPC = numeric()
  )
}


# ----------------------------------------------------------------------
# Clean IDs
task.data$participant_id <- sub("^sub-", "", task.data$participant_id)
task.data$session_id <- sub("^ses-0", "", task.data$session_id)

# Create ID
task.data$ID <- paste0(task.data$participant_id, "-",task.data$session_id)

# Save new database
write.csv(task.data, file = "/data_/mica3/BIDS_MICs/rawdata/participants_tasks-events.csv", quote = TRUE, row.names = FALSE)

# ----------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------
str(task.data)

# Missing values
colSums(is.na(task.data))

# Duplicated IDs
sum(duplicated(task.data$ID))

# Accuracy and BPS/BPC distributions
summary(task.data[, 3:8])

# Identify subjects with weird identifications (sub or ses)
task.data[!nchar(as.character(task.data$ID)) == 7,c("participant_id", "session_id")]
