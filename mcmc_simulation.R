# ============================================================
# Exercice 1 – Question 1 : Chaîne de Markov sur 3 produits
# ============================================================

# Loi initiale
pi0 <- c(P1 = 0.50, P2 = 0.30, P3 = 0.20)
cat("Loi initiale pi0 :\n")
print(pi0)
cat("Somme =", sum(pi0), "\n\n")

# Matrice de transition K
K <- matrix(
  c(0.60, 0.20, 0.20,
    0.35, 0.30, 0.35,
    0.05, 0.05, 0.90),
  nrow = 3, byrow = TRUE,
  dimnames = list(from = c("P1","P2","P3"),
                  to   = c("P1","P2","P3"))
)
cat("Matrice de transition K :\n")
print(K)

# Vérification lignes = 1
cat("\nSomme de chaque ligne :\n")
print(rowSums(K))

# ===========================================

# ============================================================
# Exercice 1 – Question 2 : Loi stationnaire
# ============================================================

# Vecteur propre de t(K) associé à la valeur propre 1
vp <- eigen(t(K))
pi_stat <- Re(vp$vectors[, 1])
pi_stat <- pi_stat / sum(pi_stat)       # normalisation
names(pi_stat) <- c("P1", "P2", "P3")

cat("Loi stationnaire :\n")
print(round(pi_stat, 4))

# Vérification : pi_stat %*% K doit redonner pi_stat
cat("\nVérification pi * K :\n")
print(round(pi_stat %*% K, 4))

# Comparaison avec la loi initiale
cat("\nLoi initiale    :", pi0, "\n")
cat("Loi stationnaire:", round(pi_stat, 4), "\n")

# ============================================================
# Exercice 1 – Question 3 : Évolution sur 20 mois
# ============================================================

n_mois <- 20
# Matrice pour stocker les distributions (1 ligne = 1 mois)
evol <- matrix(0, nrow = n_mois + 1, ncol = 3,
               dimnames = list(paste0("mois_", 0:n_mois), c("P1","P2","P3")))

evol[1, ] <- pi0  # distribution initiale

for (t in 1:n_mois) {
  evol[t+1, ] <- evol[t, ] %*% K  # multiplication par K à droite
}

print(round(evol, 4))

# Détection du mois de convergence (écart < seuil vs loi stationnaire)
seuil <- 1e-3
for (t in 1:(n_mois + 1)) {
  if (max(abs(evol[t, ] - pi_stat)) < seuil) {
    cat("\nRégime stationnaire atteint au mois", t - 1, "\n")
    break
  }
}

# Visualisation
matplot(0:n_mois, evol, type = "l", lty = 1, lwd = 2,
        col = c("steelblue","forestgreen","coral"),
        xlab = "Mois", ylab = "Proportion",
        main = "Évolution des parts de marché sur 20 mois")
abline(h = pi_stat, lty = 2, col = c("steelblue","forestgreen","coral"))
legend("right", legend = c("P1","P2","P3"),
       col = c("steelblue","forestgreen","coral"), lty = 1, lwd = 2)

# ============================================================
# Exercice 1 – Question 4 : Simulation d'une trajectoire
# ============================================================

# Fonction : tire un état (1,2,3) selon le vecteur p, via U[0,1]
tirage_discret <- function(p) {
  u <- runif(1)
  which(u < cumsum(p))[1]  # premier intervalle où u tombe
}

t = 200
# Simulation d'une trajectoire sur T pas
trajectoire <- function(T = t) {
  traj <- integer(T + 1)
  traj[1] <- tirage_discret(pi0)       # état initial selon pi0
  
  for (t in 1:T) {
    p_next <- K[traj[t], ]             # ligne de K = proba de transition
    traj[t+1] <- tirage_discret(p_next)
  }
  return(traj)
}

traj <- trajectoire(t)
cat("Trajectoire (1=P1, 2=P2, 3=P3) :\n")
print(traj)

# Visualisation
plot(0:t, traj, type = "s", lwd = 2, col = "steelblue",
     xlab = "Mois", ylab = "Produit", yaxt = "n",
     main = "Trajectoire d'un consommateur sur 20 mois")
axis(2, at = 1:3, labels = c("P1","P2","P3"))

# ===========================================
# ============================================================
# Exercice 1 – Question 5 : Trajectoire de taille n=500
# ============================================================

set.seed(42)
n <- 500
traj_long <- trajectoire(n)

# Fréquences empiriques (on exclut l'état initial)
freq_emp <- table(traj_long) / (n + 1)
names(freq_emp) <- c("P1","P2","P3")

cat("Fréquences empiriques :\n")
print(round(freq_emp, 4))

cat("\nLoi stationnaire :\n")
print(round(pi_stat, 4))

cat("\nÉcart :\n")
print(round(abs(freq_emp - pi_stat), 4))

# ===========================================
tailles <- 2:2000
ecarts  <- numeric(length(tailles))

set.seed(42)
for (i in seq_along(tailles)) {
  traj_i    <- trajectoire(tailles[i])
  freq_i    <- table(factor(traj_i, levels = 1:3)) / (tailles[i] + 1)
  ecarts[i] <- max(abs(freq_i - pi_stat))
}

plot(tailles, ecarts, type = "l",
     lwd = 1, col = "steelblue",
     xlab = "Taille n", ylab = "Écart max |freq - pi|",
     main = "Convergence vers la loi stationnaire")
abline(h = 0, lty = 2, col = "gray")
