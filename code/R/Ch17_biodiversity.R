# =====================================================================
# Chapter 17 - Ecology and Biodiversity Data (R)
# Diversity indices, rarefaction, Bray-Curtis, PCoA and PERMANOVA.
#
# This version INSTALLS vegan itself (into your personal library, so it works
# even when the system library is read-only) and, if that is not possible,
# falls back to base-R computations so the script always runs.
#
# Data: community.csv (sites x species), sites_meta.csv (site, habitat).
# =====================================================================

# ---------------------------------------------------------------------
# 0. Make sure vegan is available -- install it into a writable user library
# ---------------------------------------------------------------------
ensure_packages <- function(pkgs, repo = "https://cloud.r-project.org") {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib)) user_lib <- file.path(path.expand("~"), "R-libs")
  if (!dir.exists(user_lib))
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))            # prefer the user library
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing '", p, "' into ", user_lib, " ...")
      try(install.packages(p, lib = user_lib, repos = repo), silent = TRUE)
    }
  }
}
ensure_packages("vegan")
HAVE_VEGAN <- requireNamespace("vegan", quietly = TRUE)
if (HAVE_VEGAN) suppressMessages(library(vegan)) else
  cat("(vegan could not be installed -- using the base-R fallback below)\n")

# ---------------------------------------------------------------------
# 1. Read the community matrix and the site metadata
# ---------------------------------------------------------------------
comm <- read.csv("community.csv", row.names = 1, check.names = FALSE)
comm <- as.matrix(comm)
meta <- read.csv("sites_meta.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta[[1]]                       # first column = site id
meta <- meta[rownames(comm), , drop = FALSE]      # align to the community rows

# pick the grouping column (habitat / land_use / first text column)
gcol <- intersect(c("habitat", "land_use", "group"), names(meta))
gcol <- if (length(gcol)) gcol[1] else
        names(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))][2]
group <- factor(meta[[gcol]])
cat("sites:", nrow(comm), " species:", ncol(comm),
    " groups:", paste(levels(group), collapse = ", "), "\n\n")

# ---------------------------------------------------------------------
# 2. Alpha diversity (richness + evenness)
# ---------------------------------------------------------------------
if (HAVE_VEGAN) {
  S    <- specnumber(comm)                 # richness
  H    <- diversity(comm, "shannon")
  D1   <- diversity(comm, "simpson")       # 1 - sum(p^2)
  invD <- diversity(comm, "invsimpson")    # 1 / sum(p^2)
} else {
  S    <- rowSums(comm > 0)
  shannon_row <- function(x){ p <- x / sum(x); p <- p[p > 0]; -sum(p * log(p)) }
  simpson_row <- function(x){ p <- x / sum(x); sum(p^2) }
  H    <- apply(comm, 1, shannon_row)
  D1   <- 1 - apply(comm, 1, simpson_row)
  invD <- 1 / apply(comm, 1, simpson_row)
}
J <- H / log(S)                            # Pielou's evenness

cat("== Alpha diversity by", gcol, "(means) ==\n")
alpha <- data.frame(richness = S, shannon = H, invsimpson = invD, evenness = J)
agg <- aggregate(alpha, list(group = group), mean)
agg[-1] <- round(agg[-1], 2)
print(agg, row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------
# 3. Beta diversity: Bray-Curtis dissimilarity
# ---------------------------------------------------------------------
if (HAVE_VEGAN) {
  bc <- vegdist(comm, method = "bray")
} else {
  bray <- function(m){
    n <- nrow(m); D <- matrix(0, n, n)
    for (i in 1:(n-1)) for (j in (i+1):n) {
      d <- sum(abs(m[i, ] - m[j, ])) / sum(m[i, ] + m[j, ])
      D[i, j] <- D[j, i] <- d
    }
    as.dist(D)
  }
  bc <- bray(comm)
}

# ---------------------------------------------------------------------
# 4. Ordination: PCoA (classical MDS -- base R, no vegan needed)
# ---------------------------------------------------------------------
pco <- cmdscale(bc, k = 2, eig = TRUE)
pos <- pco$eig[pco$eig > 0]
varexp <- round(100 * pco$eig[1:2] / sum(pos), 1)
cat("== PCoA ==\n")
cat("variance explained: axis1 =", varexp[1], "%  axis2 =", varexp[2], "%\n\n")

# ---------------------------------------------------------------------
# 5. PERMANOVA: do the communities differ by group?
# ---------------------------------------------------------------------
cat("== PERMANOVA (community ~", gcol, ") ==\n")
if (HAVE_VEGAN) {
  set.seed(1)
  ad <- adonis2(comm ~ group, data = data.frame(group = group),
                method = "bray", permutations = 999)
  print(ad)
} else {
  # pseudo-F from the distance matrix + label permutations
  permanova <- function(d, g, nperm = 999, seed = 1) {
    set.seed(seed)
    Dm <- as.matrix(d); n <- nrow(Dm); g <- factor(g); k <- nlevels(g)
    sst <- sum(Dm^2) / (2 * n)
    ssw <- function(gg) {
      s <- 0
      for (lv in levels(gg)) {
        idx <- which(gg == lv); ns <- length(idx)
        if (ns > 1) s <- s + sum(Dm[idx, idx]^2) / (2 * ns)
      }
      s
    }
    Fstat <- function(gg) {
      w <- ssw(gg); b <- sst - w
      (b / (k - 1)) / (w / (n - k))
    }
    F0 <- Fstat(g)
    Fp <- replicate(nperm, Fstat(factor(sample(as.character(g)))))
    p  <- (sum(Fp >= F0) + 1) / (nperm + 1)
    cat(sprintf("pseudo-F = %.2f   p = %.3f   (df = %d, %d)\n",
                F0, p, k - 1, n - k))
  }
  permanova(bc, group)
}

cat("\nDone.\n")
