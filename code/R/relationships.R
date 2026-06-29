# =====================================================================
# Chapter 12 - Relationships Between Variables (R)
# Correlation, linear/multiple regression, non-linear regression.
#
#   source("relationships.R")
#
# Data: reef_fish.csv, growth.csv, kinetics.csv, survival.csv.
# #> lines show the verified fitted values.
# =====================================================================
fish   <- read.csv("reef_fish.csv")
growth <- read.csv("growth.csv")
kin    <- read.csv("kinetics.csv")
surv   <- read.csv("survival.csv")

cat("\n== Correlation: length vs mass ==\n")
cat("Pearson :", round(cor(fish$length_mm, fish$mass_g), 3), "\n")                    #> 0.908
cat("Spearman:", round(cor(fish$length_mm, fish$mass_g, method="spearman"), 3), "\n") #> 0.917
cat("Kendall :", round(cor(fish$length_mm, fish$mass_g, method="kendall"), 3), "\n")  #> 0.761
print(round(cor(fish[c("length_mm","mass_g","parasite_count","expression")]), 2))

cat("\n== Simple linear regression: mass ~ length ==\n")
m <- lm(mass_g ~ length_mm, data = fish)
print(coef(m))                              #> intercept -63.7, slope 0.799
cat("R^2 =", round(summary(m)$r.squared, 3), "\n")   #> 0.824

cat("\n== Multiple regression: mass ~ length + parasites + sex ==\n")
m2 <- lm(mass_g ~ length_mm + parasite_count + sex, data = fish)
print(round(coef(m2), 3))
cat("adj R^2 =", round(summary(m2)$adj.r.squared, 3), "   AIC:", round(AIC(m2)), "\n")
cat("model comparison (does adding terms help?):\n"); print(anova(m, m2))

cat("\n== Polynomial regression: density ~ time ==\n")
print(round(coef(lm(density_OD ~ poly(time_h, 2, raw=TRUE), data = growth)), 4))

cat("\n== Allometry (log-log): exponent of mass ~ length ==\n")
cat(round(coef(lm(log(mass_g) ~ log(length_mm), data = fish))[2], 2), "\n")   #> 3.06

cat("\n== Non-linear: logistic growth (nls) ==\n")
fitL <- nls(density_OD ~ K/(1+exp(-r*(time_h-t0))), data = growth,
            start = c(K=1, r=0.5, t0=12))
print(round(coef(fitL), 3))                 #> K 1.06, r 0.46, t0 11.1

cat("\n== Non-linear: Michaelis-Menten (nls) ==\n")
fitM <- nls(velocity ~ Vmax*substrate_mM/(Km+substrate_mM), data = kin,
            start = c(Vmax=2, Km=1))
print(round(coef(fitM), 3))                 #> Vmax 2.18, Km 1.57

cat("\n== Logistic regression: survived ~ body_size ==\n")
fitG <- glm(survived ~ body_size_mm, data = surv, family = binomial)
print(round(coef(fitG), 3))
cat("odds ratio per mm =", round(exp(coef(fitG)[2]), 2), "\n")

cat("\nDone.\n")
