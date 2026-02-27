library(reticulate)
library(lavaan)
library(stats)
library(RcppCNPy)

# Load npy matrices
mk6240.mat <- npyLoad("/Users/rcruces/Desktop/OSF_data_revision-2025/surf-fsLR-32k_desc-GroupDataFlippedMasked_smooth-10mm_pvc-probGM_ref-cerebellarGM_trc-18Fmk6240_pet.npy")
scneig.mat <- npyLoad("/Users/rcruces/Desktop/OSF_data_revision-2025/surf-fsLR-5k_desc-neighborsFlipped_SCw.npy")

# Assume n_vertices is number of columns in mk6240.mat
n_vertices <- ncol(mk6240.mat)

# Initialize result matrices
sem <- list(
  coef      = matrix(NA, nrow = n_vertices, ncol = 3),
  pval      = matrix(NA, nrow = n_vertices, ncol = 3),
  pval.corr = matrix(NA, nrow = n_vertices, ncol = 3)
)

# Iterate over the vertices
fit_vertex_sem <- function(v, outcome_vec) {
  predictor_vec <- scale(mk6240.mat[, v])
  mediator_vec  <- scale(scneig.mat[, v])
  
  # Complete cases
  keep <- complete.cases(predictor_vec, mediator_vec, outcome_vec)
  data <- data.frame(
    predictor = predictor_vec[keep],
    mediator  = mediator_vec[keep],
    outcome   = outcome_vec[keep]
  )
  
  # SEM model
  model <- '
    mediator ~ predictor
    outcome ~ predictor + mediator
  '
  
  fit <- sem(model, data = data)
  
  # Extract standardized coefficients and p-values
  coefs <- lavaan::standardizedSolution(fit)
  # Order: predictor->mediator, predictor->outcome, mediator->outcome
  coef_vec <- coefs$est.std[coefs$op == "~"]
  pval_vec <- coefs$pvalue[coefs$op == "~"]
  
  return(list(coef = coef_vec, pval = pval_vec))
}

# Apply over vertices
results <- lapply(1:n_vertices, function(x) fit_vertex_sem(x, mk.ses1$Epi.Acc.D.per))

# Coefficients matrix
sem$coef <- t(sapply(results, function(x) x$coef))

# P-values matrix
sem$pval <- t(sapply(results, function(x) x$pval))

# FDR correction across vertices for each path
sem$pval.corr <- apply(sem$pval, 2, p.adjust, method = "fdr")

write.table(sem$coef, "/Users/rcruces/Desktop/OSF_data_revision-2025/surf-fsLR-5k_desc-episodic_coef_SEM.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(sem$pval, "/Users/rcruces/Desktop/OSF_data_revision-2025/surf-fsLR-5k_desc-episodic_pval_SEM.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(sem$pval.corr, "/Users/rcruces/Desktop/OSF_data_revision-2025/surf-fsLR-5k_desc-episodic_pcor_SEM.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)