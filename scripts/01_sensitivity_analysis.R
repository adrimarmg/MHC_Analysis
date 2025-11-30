# --- Setup & options ----------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(MuMIn)
library(ggplot2)
library(glmnet)

set.seed(123)
options(na.action = "na.fail")


# ------------------------------------------------------------------------------ 
# Load predictors and response variables
# The script assumes that predictors and response variables are located inside a 
# 'data/' directory relative to the project root. If running the script 
# interactively, ensure that your working directory is set to the project folder
# ------------------------------------------------------------------------------

# Load predictors
coords1 <- read.csv("data/predictors.csv", header = TRUE, row.names = 1)

# Load response variable
coords <- read.csv("data/response_variables.csv", header = TRUE, row.names = 1)

# Merge the two datasets based on their row names
merged_df <- merge(coords, coords1, by = "row.names")

# Set the row names of the merged dataframe to the first column (Row.names)
rownames(merged_df) <- merged_df$Row.names

# Remove the "Row.names" column from the merged dataframe
merged_df <- merged_df[, !(names(merged_df) %in% c("Row.names"))]

# ------------------------------------------------------------------------------
# Select the response variable and the predictor variables used in 
# the multiple regression analysis
#-------------------------------------------------------------------------------

# Response variable
response_var <- "TD_WG"

# Predictor variables for the multiple regression
predictor_vars <- c(
  "Neolithic", "Urbanization", "Interval",
  "PC1Neolithic", "PC2Neolithic",
  "PC1urbanization", "PC2urbanization",
  "PC1Currenttime", "PC2Currenttime",
  "Total.Pathogen.Stress",
  "Alpha", "DistanceAfrica",
  "AUC_complex"
)

# Create final dataframe
final_df <- merged_df[, c(response_var, predictor_vars)]

# --- Columns to scale (predictor variables only) ------------------------------
# These are the predictor variables that will be standardized (e.g. scaled)
columns_to_scale <- c(
  "Neolithic", "Urbanization", "Interval",
  "PC1Neolithic", "PC2Neolithic",
  "PC1urbanization", "PC2urbanization",
  "PC1Currenttime", "PC2Currenttime",
  "Total.Pathogen.Stress",
  "Alpha",
  "DistanceAfrica",
  "AUC_complex"
)

#-------------------------------------------------------------------------------
# --- Helper: add uniform ±p% error to predictors
# Adds multiplicative noise to the predictor variables to perform a sensitivity
# analysis. For each predictor value x, the function applies:
#   x' = x * u   where u ~ Uniform(1 - level, 1 + level)
#
# Arguments:
#   data          : data.frame with response and predictors
#   level         : error level (e.g. 0.10 = ±10% error)
#   predictor_cols: character vector with the names of predictor columns
#
# Note: the response variable is not modified.

add_uniform_error <- function(data, level, predictor_cols) {
  if (level == 0) return(data)
  
  data_with_error <- data
  for (col in predictor_cols) {
    err_factor <- 1 + runif(nrow(data), -level, level)
    data_with_error[[col]] <- data[[col]] * err_factor
  }
  data_with_error
}

# --- Error levels & iterations ------------------------------------------------
# Sensitivity analysis settings:
#   error_levels     : proportion of uniform error added to predictors (0–40%)
#   num_iterations   : number of replicate datasets generated per error level

error_levels <- c(0.00, 0.10, 0.20, 0.30, 0.40)
num_iterations <- 100

# --- Init result tables -------------------------------------------------------
# List of predictor variables (excluding the response variable)
vars <- colnames(final_df)[colnames(final_df) != "TD_WG"]

# Column names for storing results across error levels
col_names <- c("0%", "10%", "20%", "30%", "40%")

# Table to store how many times each variable is selected by the LASSO model
lasso_selection_count <- data.frame(
  Variable = vars,
  `0%` = 0L, `10%` = 0L, `20%` = 0L, `30%` = 0L, `40%` = 0L,
  check.names = FALSE
)

# Table to store how many times each variable appears in the LASSO+AICc best model
best_model_selection_count <- data.frame(
  Variable = vars,
  `0%` = 0L, `10%` = 0L, `20%` = 0L, `30%` = 0L, `40%` = 0L,
  check.names = FALSE
)
# ------------------------------------------------------------------------------
# Main loop: run each error level in its own 100-iteration block
# ------------------------------------------------------------------------------
for (level in error_levels) {
  level_label <- paste0(round(level * 100), "%")
  message("Running error level: ", level_label)
  
  for (iteration in 1:num_iterations) {
    # 1) Add error for this scenario (uniform ±level across predictors)
    #    This perturbs only the predictor variables; the response is left unchanged.
    data_with_error <- add_uniform_error(
      data          = final_df,
      level         = level,
      predictor_cols = columns_to_scale
    )
    
    # 2) Scale predictors; keep response unscaled and as the last column
    #    - Standardize all predictors (mean = 0, sd = 1)
    #    - Remove the unscaled predictors
    #    - Restore original predictor names
    #    - Move the response variable to the last column
    scaled_data <- data_with_error %>%
      mutate(across(all_of(columns_to_scale), scale, .names = "scaled_{col}")) %>%
      select(-all_of(columns_to_scale)) %>%                              # drop unscaled
      rename_with(~ gsub("^scaled_", "", .), starts_with("scaled_")) %>% # clean names back
      relocate(TD_WG, .after = last_col())
    
    # 3) LASSO with cross-validation to select predictors
    #    Here, we use glmnet with alpha = 1 (LASSO) and choose lambda by CV.
    y <- scaled_data$TD_WG
    X <- as.matrix(scaled_data[, setdiff(names(scaled_data), "TD_WG")])
    
    cv_result   <- cv.glmnet(X, y, alpha = 1, lambda = seq(0.001, 0.1, length.out = 100))
    best_lambda <- cv_result$lambda.min
    
    lasso_fit   <- glmnet(X, y, alpha = 1, lambda = best_lambda)
    lasso_coefs <- coef(lasso_fit, s = best_lambda)
    selected_vars_lasso <- rownames(lasso_coefs)[as.numeric(lasso_coefs) != 0]
    selected_vars_lasso <- setdiff(selected_vars_lasso, "(Intercept)")
    
    # 4) Update LASSO selection counts
    if (length(selected_vars_lasso) > 0) {
      idx <- match(selected_vars_lasso, lasso_selection_count$Variable)
      lasso_selection_count[idx, level_label] <-
        lasso_selection_count[idx, level_label] + 1L
    }
    
    # 5) Fit linear model with LASSO-selected predictors and apply AICc dredge
    #    We fit an lm model using only the LASSO-selected variables, then use
    #    MuMIn::dredge to explore all possible submodels and identify the
    #    best model according to AICc. We then count which variables appear
    #    in the best AICc model.
    if (length(selected_vars_lasso) > 0) {
      form <- as.formula(
        paste("TD_WG ~", paste(selected_vars_lasso, collapse = " + "))
      )
      
      best_mod <- tryCatch(
        lm(form, data = scaled_data, na.action = na.fail),
        error = function(e) NULL
      )
      
      if (!is.null(best_mod)) {
        model_set <- tryCatch(dredge(best_mod), error = function(e) NULL)
        
        if (!is.null(model_set) && nrow(model_set) > 0) {
          best_model <- get.models(model_set, 1)[[1]]
          best_vars  <- names(coef(best_model))[-1]  # drop intercept
          
          if (length(best_vars) > 0) {
            idx2 <- match(best_vars, best_model_selection_count$Variable)
            best_model_selection_count[idx2, level_label] <-
              best_model_selection_count[idx2, level_label] + 1L
          }
        }
      }
    }
  }
}

# --- Results ------------------------------------------------------------------

# Print raw selection counts (LASSO and best AICc model)
print(lasso_selection_count)
print(best_model_selection_count)

# Save raw counts (recommended)
write.csv(lasso_selection_count, "lasso_selection_counts_by_error.csv", row.names = FALSE)
write.csv(best_model_selection_count, "best_model_selection_counts_by_error.csv", row.names = FALSE)