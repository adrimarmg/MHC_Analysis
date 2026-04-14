# ============================================================
# Area under the curve (AUC) for each  “Natural Geographic Areas” NGA from Turchin et al. 2017
# ============================================================

# Install packages if needed
# install.packages("tidyverse")
# install.packages("pracma")

library(tidyverse)
library(pracma)

# ------------------------------------------------------------
# 1. Read the CSV file
# ------------------------------------------------------------
df <- read.csv("sociopolitical_timeseries.csv", header = TRUE, stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 2. Make sure Time, V1 (PC1), and V2 (PC2) are numeric
# ------------------------------------------------------------
df <- df %>%
  mutate(
    Time = as.numeric(Time),
    V1   = as.numeric(V1),
    V2   = as.numeric(V2)
  )

# ------------------------------------------------------------
# 3. Rescale V1 and V2 from 0 to 10
# ------------------------------------------------------------
# This rescales each variable using its global minimum and maximum
# across all societies

df <- df %>%
  mutate(
    V1_scaled = (V1 - min(V1, na.rm = TRUE)) /
      (max(V1, na.rm = TRUE) - min(V1, na.rm = TRUE)) * 10,
    
    V2_scaled = (V2 - min(V2, na.rm = TRUE)) /
      (max(V2, na.rm = TRUE) - min(V2, na.rm = TRUE)) * 10
  )

# View data
print(df)

# ------------------------------------------------------------
# 4. Calculate AUC for each “Natural Geographic Areas” NGA
# ------------------------------------------------------------
auc_results <- df %>%
  group_by(NGA) %>%
  arrange(Time, .by_group = TRUE) %>%
  summarise(
    AUC_V1 = trapz(Time, V1_scaled),
    AUC_V2 = trapz(Time, V2_scaled),
    Time_span = max(Time, na.rm = TRUE) - min(Time, na.rm = TRUE),
    Mean_AUC_V1 = AUC_V1 / Time_span,
    Mean_AUC_V2 = AUC_V2 / Time_span,
    .groups = "drop"
  )

# View results
print(auc_results)

# ------------------------------------------------------------
# 5. Save results
# ------------------------------------------------------------
write.csv(auc_results, "AUC_results_V1_V2.csv", row.names = FALSE)
