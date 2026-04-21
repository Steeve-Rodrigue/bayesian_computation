# ============================================================
# Exercice 2 – Question 4 : Échantillonneur de Gibbs
# ============================================================
library(coda)
# --- Données ---
y <- c(1.80,1.85,1.87,1.77,2.02,2.27,2.15,2.26,2.47,2.19,
       2.26,2.40,2.39,2.41,2.50,2.32,2.32,2.43,2.47,2.56,
       2.65,2.47,2.64,2.56,2.70,2.72,2.57)
x <- c(1.0,1.5,1.5,1.5,2.5,4.0,5.0,5.0,7.0,8.0,8.5,9.0,
       9.5,9.5,10.0,12.0,12.0,13.0,13.0,14.5,15.5,15.5,
       16.5,17.0,22.5,29.0,31.5)

n    <- length(y)
logx <- log(x)
XtX  <- sum(logx^2)

# --- Algorithme de Gibbs ---
gibbs_sirenia <- function(L = 5000) {
  
  # Initialisation par MCO
  fit  <- lm(y ~ logx)
  beta <- matrix(0, nrow = L, ncol = 3,
                 dimnames = list(NULL, c("beta1","beta2","tau")))
  
  #inititalisation avec l'estimateur MCO
  #beta[1, ] <- c(coef(fit), 1/var(residuals(fit)))
  
  #inititalisation naive
  beta[1, ] <- c(0,0,1)
  
  for (t in 2:L) {
    b1  <- beta[t-1, "beta1"]
    b2  <- beta[t-1, "beta2"]
    tau <- beta[t-1, "tau"]
    
    # Étape 1 : beta1 | beta2, tau, y
    m1 <- mean(y - b2 * logx)
    b1 <- rnorm(1, mean = m1, sd = sqrt(1 / (n * tau)))
    
    # Étape 2 : beta2 | beta1, tau, y
    m2 <- sum(logx * (y - b1)) / XtX
    b2 <- rnorm(1, mean = m2, sd = sqrt(1 / (tau * XtX)))
    
    # Étape 3 : tau | beta1, beta2, y
    RSS    <- sum((y - b1 - b2 * logx)^2)
    tau    <- rgamma(1, shape = n/2 + 1, rate = RSS/2)
    
    beta[t, ] <- c(b1, b2, tau)
  }
  return(beta)
}

# --- Exécution ---
set.seed(42)
L      <- 10000
res    <- gibbs_sirenia(L = L)
chaine <- mcmc(res)

# ===========================================
# ============================================================
# Exercice 2 – Question 4 : Gibbs – Historique des chaînes
# ============================================================

par(mfrow = c(2, 2))
traceplot(chaine)

# Zoom début (burn-in)
par(mfrow = c(2, 2))
traceplot(window(chaine, end = 100))

# Zoom fin (régime stationnaire)
par(mfrow = c(2, 2))
traceplot(window(chaine, start = 800))

# ============================================================
# Exercice 2 – Question 4 : Gibbs – Quantiles ergodiques
# ============================================================

par(mfrow = c(2, 2))
cumuplot(window(chaine, start = 1))

# ============================================================
# Exercice 2 – Question 4 : Gibbs – Densités
# ============================================================

par(mfrow = c(2, 2))
densplot(window(chaine, start = 1))

# ============================================================
# Exercice 2 – Question 4 : Gibbs – Autocorrélations
# ============================================================

par(mfrow = c(2, 2))
autocorr.plot(window(chaine, start = 1))

# ============================================================
# Exercice 2 – Question 4 : Gibbs – Erreur MCMC
# ============================================================

summary(window(chaine, start = 1))

# ===========================================

# ============================================================
# Exercice 2 – Question 5 : Statistiques a posteriori
# ============================================================

chaine_clean <- window(chaine, start = 10)
res_clean    <- as.matrix(chaine_clean)

beta1       <- res_clean[, "beta1"]
beta2       <- res_clean[, "beta2"]
tau         <- res_clean[, "tau"]
sigma2_post <- 1 / tau

params      <- list(beta1 = beta1, beta2 = beta2, 
                    tau = tau, "1/tau" = sigma2_post)

par(mfrow = c(2, 2))
for (nom in names(params)) {
  p <- params[[nom]]
  cat("\n===", nom, "===\n")
  cat("Moyenne      :", round(mean(p), 4), "\n")
  cat("Écart-type   :", round(sd(p),   4), "\n")
  cat("IC 95%       :", round(quantile(p, c(0.025, 0.975)), 4), "\n")
  densplot(mcmc(p), main = paste("Densité a posteriori de", nom))
}

# ===========================================
# ============================================================
# Exercice 2 – Question 6 : Comparaison avec lm
# ============================================================


fit <- lm(y ~ logx)
summary(fit)

# Comparaison des estimations
cat("\n=== Comparaison Bayes vs fréquentiste ===\n")
cat("\nbeta1 :\n")
cat("  Bayes :", round(mean(beta1), 4), 
    "| IC 95% :", round(quantile(beta1, c(0.025, 0.975)), 4), "\n")
cat("  MCO   :", round(coef(fit)[1], 4),
    "| IC 95% :", round(confint(fit)[1,], 4), "\n")

cat("\nbeta2 :\n")
cat("  Bayes :", round(mean(beta2), 4),
    "| IC 95% :", round(quantile(beta2, c(0.025, 0.975)), 4), "\n")
cat("  MCO   :", round(coef(fit)[2], 4),
    "| IC 95% :", round(confint(fit)[2,], 4), "\n")

cat("\n1/tau :\n")
cat("  Bayes :", round(mean(sigma2_post), 4), "\n")
cat("  MCO   :", round(var(residuals(fit)), 4), "\n")

# ===========================================
# ============================================================
# Exercice 2 – Question 7 : Nuage de points beta1 vs beta2
# ============================================================

plot(beta1, beta2, pch = 1, cex = 0.2,
     col = rgb(0, 0, 1, 0.3),
     xlab = expression(beta[1]),
     ylab = expression(beta[2]),
     main = expression(paste("Nuage de points MCMC : ", 
                             beta[1], " vs ", beta[2])))

# Ajout des moyennes a posteriori
points(mean(beta1), mean(beta2), pch = 3, cex = 1, 
       col = "red", lwd = 2)

# ===========================================