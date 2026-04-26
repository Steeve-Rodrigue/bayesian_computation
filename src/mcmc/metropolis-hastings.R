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

#==================================================================================
#=======================================================================
# ============================================================
# Metropolis-Hastings – Loi a priori cos²(4πθ) sur [0,1]
# ============================================================

# --- Données (à adapter) ---
n <- 20   # nombre d'essais
x <- 14   # nombre de succès observés
x_obs <- x

# --- Log a posteriori (pour éviter les underflows) ---
log_post <- function(theta) {
  if (theta <= 0 || theta >= 1) return(-Inf)
  x * log(theta) + (n - x) * log(1 - theta) + 2 * log(abs(cos(4 * pi * theta)))
}

# --- Algorithme MH ---
mh_beta <- function(theta0, L, tau) {
  chain    <- numeric(L)
  chain[1] <- theta0
  n_accept <- 0
  
  for (t in 2:L) {
    # Proposition gaussienne
    theta_new <- rnorm(1, mean = chain[t-1], sd = 1/tau)
    
    # Rejet direct si hors support
    if (theta_new <= 0 || theta_new >= 1) {
      chain[t] <- chain[t-1]
      next
    }
    
    # Taux d'acceptation (sur log-échelle)
    log_rho <- log_post(theta_new) - log_post(chain[t-1])
    
    if (log(runif(1)) <= log_rho) {
      chain[t] <- theta_new
      n_accept <- n_accept + 1
    } else {
      chain[t] <- chain[t-1]
    }
  }
  
  cat("Taux d'acceptation :", round(n_accept / L, 3), "\n")
  return(chain)
}

# --- Simulation pour différentes valeurs de tau ---
set.seed(42)
L <- 5000

chain_tau1  <- mh_beta(theta0 = 0.5, L = L, tau = 1)
chain_tau5  <- mh_beta(theta0 = 0.5, L = L, tau = 5)
chain_tau50 <- mh_beta(theta0 = 0.5, L = L, tau = 50)

# --- Visualisation ---
par(mfrow = c(3, 2))

for (chain in list(chain_tau1, chain_tau5, chain_tau50)) {
  lab <- switch(which(sapply(list(chain_tau1, chain_tau5, chain_tau50),
                             identical, chain)),
                "tau=1", "tau=5", "tau=50")
  plot(chain, type = "l", col = "steelblue", 
       main = paste("Traceplot –", lab), 
       xlab = "Itération", ylab = expression(theta))
  hist(chain, breaks = 50, freq = FALSE, col = "lightblue",
       main = paste("Densité –", lab),
       xlab = expression(theta))
  curve(cos(4*pi*x)^2 * x^x_obs * (1-x)^(n-x_obs), 
        add = TRUE, col = "red", lwd = 2)
}

par(mfrow = c(1, 1))

# ===========================================

# ============================================================
# Étude de la convergence
# ============================================================
library(coda)

# Conversion en objets mcmc
mcmc_tau1  <- mcmc(chain_tau1)
mcmc_tau5  <- mcmc(chain_tau5)
mcmc_tau50 <- mcmc(chain_tau50)

# --- Traceplots ---
par(mfrow = c(2, 2))
traceplot(mcmc_tau1,  main = "tau=1")
traceplot(mcmc_tau5,  main = "tau=5")
traceplot(mcmc_tau50, main = "tau=50")

# --- Zoom début (burn-in) ---
par(mfrow = c(2, 2))
traceplot(window(mcmc_tau1,  end = 200), main = "tau=1  [début]")
traceplot(window(mcmc_tau5,  end = 200), main = "tau=5  [début]")
traceplot(window(mcmc_tau50, end = 200), main = "tau=50 [début]")

# --- Quantiles ergodiques ---
par(mfrow = c(2, 2))
cumuplot(window(mcmc_tau1,  start = 1), main = "tau=1")
cumuplot(window(mcmc_tau5,  start = 1), main = "tau=5")
cumuplot(window(mcmc_tau50, start = 1), main = "tau=50")

# --- Densités ---
par(mfrow = c(2, 2))
densplot(window(mcmc_tau1,  start = 500), main = "tau=1")
densplot(window(mcmc_tau5,  start = 500), main = "tau=5")
densplot(window(mcmc_tau50, start = 500), main = "tau=50")

# --- Autocorrélations ---
par(mfrow = c(2, 2))
autocorr.plot(window(mcmc_tau1,  start = 500), main = "tau=1")
autocorr.plot(window(mcmc_tau5,  start = 500), main = "tau=5")
autocorr.plot(window(mcmc_tau50, start = 500), main = "tau=50")

# --- Erreur MCMC ---
cat("\n=== Résumé tau=1 ===\n");  print(summary(window(mcmc_tau1,  start = 500)))
cat("\n=== Résumé tau=5 ===\n");  print(summary(window(mcmc_tau5,  start = 500)))
cat("\n=== Résumé tau=50 ===\n"); print(summary(window(mcmc_tau50, start = 500)))

# ===========================================

# ============================================================
# Exercice 2 – Question 2 : MH pour loi de Poisson
# ============================================================

mh_poisson <- function(lambda = 5, L = 5000, x0 = 0) {
  chain    <- integer(L)
  chain[1] <- x0
  
  for (t in 2:L) {
    x_old <- chain[t-1]
    Z     <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5))
    y     <- x_old + Z
    
    if (y < 0) { chain[t] <- x_old; next }
    
    rho <- min(1, exp(y*log(lambda) - lfactorial(y) -
                        x_old*log(lambda) + lfactorial(x_old)))
    
    chain[t] <- ifelse(runif(1) <= rho, y, x_old)
  }
  return(chain)
}

# ============================================================
# Exercice 2 – Question 3 : Comparaison MCMC vs rpois
# ============================================================

set.seed(42)
L      <- 5000
lambda <- 5

chain <- mh_poisson(lambda = lambda, L = L, x0 = 0)
iid   <- rpois(L, lambda = lambda)

# Comparaison graphique
par(mfrow = c(1, 1))
valeurs   <- 0:15
freq_mcmc <- table(factor(chain, levels = valeurs)) / L
freq_iid  <- table(factor(iid,   levels = valeurs)) / L
freq_theo <- dpois(valeurs, lambda = lambda)

barplot(rbind(freq_mcmc, freq_iid, freq_theo),
        beside = TRUE, col = c("steelblue","coral","gray"),
        main = "MCMC vs rpois vs théorique",
        xlab = "x", ylab = "Fréquence",
        legend.text = c("MCMC","rpois","théorique"))

# Comparaison numérique
cat("Moyenne  - MCMC:", round(mean(chain),4), 
    "| rpois:", round(mean(iid),4), 
    "| théo:", lambda, "\n")
cat("Variance - MCMC:", round(var(chain),4), 
    "| rpois:", round(var(iid),4),  
    "| théo:", lambda, "\n")

# ===========================================
# ============================================================
# Comparaison MCMC vs rpois
# ============================================================

par(mfrow = c(2, 2))

# Historiques
plot(chain, type = "l", col = "steelblue",
     main = "Historique MCMC",
     xlab = "Itération", ylab = "x")

plot(iid, type = "l", col = "coral",
     main = "Historique rpois (iid)",
     xlab = "Itération", ylab = "x")

# Distributions
valeurs   <- 0:15
freq_mcmc <- table(factor(chain, levels = valeurs)) / L
freq_iid  <- table(factor(iid,   levels = valeurs)) / L
freq_theo <- dpois(valeurs, lambda = lambda)

barplot(rbind(freq_mcmc, freq_theo),
        beside = TRUE, col = c("steelblue","gray"),
        main = "MCMC vs théorique",
        xlab = "x", ylab = "Fréquence",
        legend.text = c("MCMC","théorique"))

barplot(rbind(freq_iid, freq_theo),
        beside = TRUE, col = c("coral","gray"),
        main = "rpois vs théorique",
        xlab = "x", ylab = "Fréquence",
        legend.text = c("rpois","théorique"))

par(mfrow = c(1, 1))