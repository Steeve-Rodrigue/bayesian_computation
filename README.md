# Bayesian Computation Toolkit

## Overview

This repository provides practical implementations of core Bayesian computation techniques, with a focus on sampling-based methods such as Markov Chain Monte Carlo (MCMC).

The goal is to bridge theory and practice by offering clear, reproducible code for Bayesian inference.


---

## Features

* Markov Chain simulation and stationary distribution
* Metropolis-Hastings algorithm
* Gibbs sampler
* Posterior estimation and credibility intervals
* Convergence diagnostics (traceplots, ergodic quantiles, autocorrelations, MCMC error)
* Comparison with frequentist methods

---

## Project Structure

```
bayesian-computation-toolkit/
│
├── data/                        # Datasets
│
├── src/                         # Core implementations
├───├ mcmc/                         # Core implementations
│      ├── mcmc.R           # Markov chain simulation 
│      ├── metropolis.R              # Metropolis-Hastings 
│      ├── gibbs.R                          # Gibbs sampler – ANOVA
││    ├── ........                         
├── diagnostics/                 # Convergence diagnostics
│   └── coda_tools.R             # traceplot, cumuplot, autocorr, batchSE
│
├── README.md
└── LICENSE
```

---

## Installation

```bash
git clone https://github.com/your-username/bayesian-computation-toolkit.git
cd bayesian-computation-toolkit
```

Install required R packages :

```r
install.packages(c("coda"))
```

---

## Applications

* One-way ANOVA under Bayesian framework
* Bayesian linear regression with improper prior
* Posterior sampling for non-standard distributions
* Comparison of Bayesian vs frequentist estimates

---

## Roadmap

* [ ] Metropolis-Hastings for Poisson distribution
* [ ] Hamiltonian Monte Carlo (HMC)
* [ ] Variational Inference
* [ ] Advanced diagnostics (ESS, R-hat)
* [ ] Shiny dashboard for interactive diagnostics

---

## License

MIT License

---

## Author

ENSAI engineering student  
Bayesian Computation – 2026