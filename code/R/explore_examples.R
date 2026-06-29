# =====================================================================
# Chapter 6 - Exploring Your Data (R)
# Summary statistics, frequency tables, and quick checks (reef-fish data).
#
#   install.packages(c("tidyverse", "moments"))
#   source("explore_examples.R")
# =====================================================================
library(tidyverse)
library(moments)            # skewness(), kurtosis()

fish <- read_csv("reef_fish.csv")

# ---------------------------------------------------------------------
# A FIRST LOOK
# ---------------------------------------------------------------------
dim(fish)                                            # 300 7
summary(fish)

# ---------------------------------------------------------------------
# 1. MEASURES OF LOCATION
# ---------------------------------------------------------------------
cat("\n--- LOCATION ---\n")
cat("length mean/median:", mean(fish$length_mm), median(fish$length_mm), "\n")
cat("mass mean/median:", mean(fish$mass_g), median(fish$mass_g), "\n")
cat("mass geometric mean:", exp(mean(log(fish$mass_g))), "\n")        # 30.7
cat("mass 10% trimmed mean:", mean(fish$mass_g, trim = 0.1), "\n")    # 31.8
cat("parasite geometric mean:", exp(mean(log(fish$parasite_count))), "\n")  # 4.5
cat("most common habitat:", names(which.max(table(fish$habitat))), "\n")    # Reef
cat("weighted mean length (by mass):", weighted.mean(fish$length_mm, fish$mass_g), "\n")  # 123.9
print(fish %>% group_by(habitat) %>% summarise(mean_length = round(mean(length_mm), 1)))

# ---------------------------------------------------------------------
# 2. MEASURES OF SPREAD
# ---------------------------------------------------------------------
cat("\n--- SPREAD ---\n")
cat("range:", range(fish$length_mm), " width:", diff(range(fish$length_mm)), "\n")
cat("var/sd:", var(fish$length_mm), sd(fish$length_mm), "\n")
print(quantile(fish$length_mm, c(.25, .5, .75)))
cat("IQR:", IQR(fish$length_mm), "\n")
cat("MAD (scaled):", mad(fish$length_mm), "\n")
cat("SEM:", sd(fish$length_mm) / sqrt(nrow(fish)), "\n")
cat("CV length/mass %:", 100*sd(fish$length_mm)/mean(fish$length_mm),
    100*sd(fish$mass_g)/mean(fish$mass_g), "\n")
cat("skewness length/mass/parasite:", skewness(fish$length_mm),
    skewness(fish$mass_g), skewness(fish$parasite_count), "\n")
cat("excess kurtosis length/parasite:", kurtosis(fish$length_mm) - 3,
    kurtosis(fish$parasite_count) - 3, "\n")
cat("parasite skew raw vs log:", skewness(fish$parasite_count),
    skewness(log(fish$parasite_count)), "\n")        # 2.88 -> -0.15

# ---------------------------------------------------------------------
# 3. COUNTS AND FREQUENCIES
# ---------------------------------------------------------------------
cat("\n--- COUNTS ---\n")
print(table(fish$habitat))
print(round(prop.table(table(fish$habitat)) * 100, 1))
print(table(fish$sex, fish$habitat))
classes <- cut(fish$length_mm, c(0, 100, 120, 140, 1000),
               labels = c("small", "medium", "large", "xlarge"))
print(table(classes))
p <- prop.table(table(fish$habitat)); H <- -sum(p * log(p))
cat("Shannon H / evenness:", round(H, 3), round(H / log(length(p)), 3), "\n")   # 1.684 0.94
print(fish %>% group_by(habitat) %>%
  summarise(n = n(), mean_length = round(mean(length_mm), 1),
            mean_mass = round(mean(mass_g), 1),
            median_parasites = median(parasite_count)))

# ---------------------------------------------------------------------
# 4. QUICK CHECKS BEFORE ANALYSIS
# ---------------------------------------------------------------------
cat("\n--- QUICK CHECKS ---\n")
cat("missing:", sum(is.na(fish)), " duplicates:", sum(duplicated(fish)), "\n")
cat("group sizes:", range(table(fish$habitat)), "\n")
print(round(cor(fish[, c("length_mm","mass_g","parasite_count","expression")]), 2))
z <- scale(fish$length_mm)
cat("|z|>2 lengths:", sum(abs(z) > 2), "\n")
