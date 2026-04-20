# Bayesian Computation Toolkit

## Overview

This repository provides practical implementations of core Bayesian computation techniques, with a focus on sampling-based methods such as Markov Chain Monte Carlo (MCMC).

The goal is to bridge theory and practice by offering clear, reproducible code for Bayesian inference.

---

## Features

* Metropolis-Hastings algorithm
* Gibbs sampling
* Posterior estimation
* Convergence diagnostics
* Visualization tools

---

## Project Structure

```
bayesian-computation-toolkit/
│
├── data/                  # Datasets (if applicable)
├── notebooks/             # Jupyter notebooks (experiments & demos)
├── src/                   # Core implementation
│   ├── mh.py              # Metropolis-Hastings
│   ├── gibbs.py           # Gibbs sampler
│   ├── utils.py           # Helper functions
│
├── tests/                 # Unit tests
├── README.md
├── requirements.txt
├── LICENSE
└── .gitignore
```

---

## Installation

```bash
git clone https://github.com/your-username/bayesian-computation-toolkit.git
cd bayesian-computation-toolkit
pip install -r requirements.txt
```

---

## Example

```python
from src.mh import metropolis_hastings

samples = metropolis_hastings(
    target_density,
    initial_state=0,
    n_samples=10000,
    step_size=1.0
)
```

---

## Applications

* Bayesian inference
* Posterior sampling
* Probabilistic modeling
* Statistical learning

---

## Roadmap

* [ ] Hamiltonian Monte Carlo (HMC)
* [ ] Variational Inference
* [ ] Advanced diagnostics (ESS, R-hat)

---

## License

MIT License

---

## Author

Your Name
