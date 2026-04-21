# ============================================================
# Exercice 3 – Question 1 : Chaîne MCMC via Metropolis-Hastings
# ============================================================

mcmc_normal <- function(x0, L, s) {
  chain <- numeric(L)
  chain[1] <- x0
  
  for (t in 2:L) {
    # Proposition : uniforme sur [x_old - 2s, x_old + 2s]
    y <- runif(1, chain[t-1] - 2*s, chain[t-1] + 2*s)
    
    # Taux d'acceptation (ratio simplifié car loi instrumentale symétrique)
    rho <- min(1, exp(-0.5 * (y^2 - chain[t-1]^2)))
    
    # Acceptation / rejet
    if (runif(1) <= rho) {
      chain[t] <- y        # accepté
    } else {
      chain[t] <- chain[t-1]  # rejeté
    }
  }
  return(chain)
}

# ===========================================

# ============================================================
# Exercice 3 – Question 2 : Génération des 4 échantillons
# ============================================================
library(coda)

L <- 1600

# 1. Échantillon i.i.d. de la loi cible
iid       <- rnorm(L)

# 2. MCMC x0 = -10, s = 0.1 (sauts très petits)
mcmc_01   <- mcmc_normal(x0 = -10, L = L, s = 0.1)

# 3. MCMC x0 = -10, s = 0.8 (sauts intermédiaires)
mcmc_08   <- mcmc_normal(x0 = -10, L = L, s = 0.8)

# 4. MCMC x0 = -10, s = 2   (sauts larges)
mcmc_2    <- mcmc_normal(x0 = -10, L = L, s = 2)

# Assemblage en matrice mcmc (L x 4)
echantillons <- mcmc(matrix(
  c(iid, mcmc_01, mcmc_08, mcmc_2),
  nrow = L, ncol = 4,
  dimnames = list(NULL, c("iid", "s=0.1", "s=0.8", "s=2"))
))

# ===========================================
# ============================================================
# Exercice 3 – Question 3 : Historique des chaînes
# ============================================================

par(mfrow = c(2, 2))

# Tracé complet
traceplot(echantillons)

# Zoom sur les 200 premières valeurs (burn-in visible)
traceplot(window(echantillons, end = 100))

traceplot(window(echantillons, end = 1600,  start = 50))
# ============================================================
# Exercice 3 – Question 4 : Suppression du burn-in
# ============================================================

ech_burnin <- list(
  window(echantillons[, "iid"],   start = 1),
  window(echantillons[, "s=0.1"], start = 500),
  window(echantillons[, "s=0.8"], start = 100),
  window(echantillons[, "s=2"],   start = 50)
)

traceplot(ech_burnin[1])
traceplot(ech_burnin[2])
traceplot(ech_burnin[3])
traceplot(ech_burnin[4])

# ============================================================
# Exercice 3 – Question 5 : Quantiles ergodiques
# ============================================================
par(mfrow = c(2, 2))
cumuplot(window(echantillons,   start = 1))
cumuplot(window(echantillons, start = 500))
cumuplot(window(echantillons, start = 100))
cumuplot(window(echantillons,   start = 50))

# ============================================================
# Exercice 3 – Question 6 : Estimation de la densité
# ============================================================

densplot(window(echantillons, start = 1))
densplot(window(echantillons, start = 500))
densplot(window(echantillons, start = 100))
densplot(window(echantillons, start = 50))
# ============================================================
# Exercice 3 – Question 7 : Autocorrélations
# ============================================================

autocorr.plot(window(echantillons, start = 1))
autocorr.plot(window(echantillons, start = 500))
autocorr.plot(window(echantillons, start = 100))
autocorr.plot(window(echantillons, start = 50))
# ============================================================
# Exercice 3 – Question 8 : Erreur MCMC vs Monte Carlo
# ============================================================

# Erreur Monte Carlo classique (iid) : sigma/sqrt(n)
sigma2 <- 1  # Var(X) sous N(0,1)
n_iid  <- nrow(echantillons)
err_MC <- sqrt(sigma2 / n_iid)
cat("Erreur Monte Carlo classique :", round(err_MC, 6), "\n")
# ============================================================
# Exercice 3 – Question 8 : Erreur MCMC via summary
# ============================================================

summary(window(echantillons, start = 1))
summary(window(echantillons, start = 500))
summary(window(echantillons, start = 100))
summary(window(echantillons, start = 50))