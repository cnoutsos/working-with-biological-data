# =====================================================================
# Chapter 5 - Cleaning Messy Data (R)
# Runnable examples on a messy reef-fish field dataset (ships with the book).
#
#   install.packages("tidyverse")
#   source("clean_examples.R")
# =====================================================================
library(tidyverse)

# ---------------------------------------------------------------------
# 1. READ CORRECTLY  (turn placeholders into real missing values)
# ---------------------------------------------------------------------
raw <- read_csv("specimens_messy.csv")
cat("raw rows:", nrow(raw), "\n")                       # 123
cat("messy habitat spellings:", n_distinct(raw$habitat), "\n")  # 18

fish <- read_csv("specimens_messy.csv", na = c("", "NA", "-99"))

# ---------------------------------------------------------------------
# 2. MISSING VALUES
# ---------------------------------------------------------------------
print(colSums(is.na(fish)))                             # length 3, mass 2
cat("rows with any missing:", sum(!complete.cases(fish)), "\n")        # 5

# ---------------------------------------------------------------------
# 3. DUPLICATES
# ---------------------------------------------------------------------
cat("duplicate rows:", sum(duplicated(fish)), "\n")     # 3
dup_ids <- fish %>% group_by(specimen_id) %>% filter(n() > 1) %>%
  pull(specimen_id) %>% unique()
print(dup_ids)                                          # F010 F025 F060
fish <- distinct(fish)
cat("after distinct():", nrow(fish), "rows\n")          # 120

# ---------------------------------------------------------------------
# 4. STANDARDISE TEXT COLUMNS
# ---------------------------------------------------------------------
fish <- fish %>% mutate(
  habitat = str_to_title(str_trim(habitat)),
  sex     = str_sub(str_to_upper(sex), 1, 1))
print(count(fish, habitat))
print(count(fish, sex))                                 # F 59, M 61

# ---------------------------------------------------------------------
# 5. OUTLIERS  (IQR rule) and impute with the median
# ---------------------------------------------------------------------
q   <- quantile(fish$length_mm, c(.25, .75), na.rm = TRUE)
iqr <- q[2] - q[1]
out <- filter(fish, length_mm < q[1] - 1.5*iqr | length_mm > q[2] + 1.5*iqr)
print(select(out, specimen_id, length_mm))              # F026 F050 F099

fish <- fish %>% mutate(
  length_mm = if_else(length_mm < 50 | length_mm > 300, NA_real_, length_mm),
  mass_g    = if_else(mass_g < 5    | mass_g > 150,   NA_real_, mass_g),
  length_mm = replace_na(length_mm, median(length_mm, na.rm = TRUE)),
  mass_g    = replace_na(mass_g, median(mass_g, na.rm = TRUE)))

# ---------------------------------------------------------------------
# 6. DERIVED COLUMNS + VALIDATION  (Fulton's condition factor K)
# ---------------------------------------------------------------------
fish <- fish %>%
  rename(length = length_mm, mass = mass_g) %>%
  mutate(K = round(100000 * mass / length^3, 2),
         K_class = cut(K, breaks = c(0, 1.0, 1.5, 2.0, 100),
                       labels = c("poor", "fair", "good", "high")))
stopifnot(all(fish$length >= 50 & fish$length <= 300),
          all(fish$sex %in% c("M", "F")),
          !any(duplicated(fish$specimen_id)))
print(count(fish, K_class))

# ---------------------------------------------------------------------
# 7. RESHAPE  (wide tissue expression -> long)
# ---------------------------------------------------------------------
long <- read_csv("expression_wide.csv") %>%
  pivot_longer(c(gill, liver, muscle), names_to = "tissue", values_to = "expression")
cat("long dims:", dim(long), "\n")                      # 360 3

# ---------------------------------------------------------------------
# 8. JOINS  (specimens + sequencing metadata)
# ---------------------------------------------------------------------
gen <- read_csv("genetics.csv")
cat("inner:", nrow(inner_join(fish, gen, by = "specimen_id")), "\n")   # 100
left <- left_join(fish, gen, by = "specimen_id")
cat("left:", nrow(left), "| no sequencing record:", sum(is.na(left$read_count)), "\n")  # 120 | 20
cat("full:", nrow(full_join(fish, gen, by = "specimen_id")), "\n")     # 140
cat("anti:", nrow(anti_join(fish, gen, by = "specimen_id")), "\n")     # 20

# ---------------------------------------------------------------------
# 9. A COMBINED ANSWER  (reshape + join + group)
# ---------------------------------------------------------------------
long %>%
  inner_join(select(fish, specimen_id, habitat), by = "specimen_id") %>%
  group_by(habitat) %>%
  summarise(mean_expr = round(mean(expression), 1)) %>%
  print()

# ---------------------------------------------------------------------
# 10. SAVE THE CLEAN RESULT (never overwrite the raw file)
# ---------------------------------------------------------------------
clean <- left_join(fish, gen, by = "specimen_id")
write_csv(clean, "specimens_clean.csv")
cat("wrote specimens_clean.csv:", nrow(clean), "rows\n")
