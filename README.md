# thermodynamic_modeling_of_NP_endocytosis
This repository contains the Python implementation of a thermodynamic model to evaluate how nanoparticle (NP) size influences membrane wrapping during endocytosis, as described in our research.

## Overview
The model is built upon a free-energy framework incorporating:
- Mixing entropy of NP mixtures
- Ligand-receptor binding energy
- Membrane bending energy
- NP configurational entropy (simplified to retain only size-dependent energetic and entropic terms)

The free-energy functional used in the simulation is:  
$$W = M_f[\xi_f \ln(\xi_f) + (1-\xi_f)\ln(1-\xi_f)] - \mu L_b + 8\pi\kappa\eta N$$  

Where all energies are expressed in units of $k_BT$:
- $M_f$: Free membrane area (in units of receptor area $A_0$), $M_f = M - M_b$ (total membrane $M = 4\pi R^2/A_0$)
- $\xi_f$: Density of free receptors, $\xi_f = (\xi_0M - L_b)/M_f$
- $\mu$: Chemical energy gain from ligand-receptor binding
- $L_b = \eta NK$: Number of ligand-receptor bonds (scaled linearly with wrapped area)
- $\kappa$: Membrane bending rigidity (varied for clathrin-independent (CIE, $\kappa=10k_BT$) and clathrin-mediated (CME, $\kappa=60k_BT$) endocytosis)
- $\eta \in [0,1]$: Wrapping area fraction of NPs
- $N$: Number of NPs with wrapping fraction $\eta$ (contributes $8\pi\kappa\eta$ to total bending energy)
- $K = 4\pi R^2/A_0$: NP surface area in units of receptor area $A_0$
- $\xi_0$: Initial receptor density on the unbound membrane

For each NP size, the free energy $W$ is numerically minimized over $N$ (number of NPs) and $\eta$ (wrapping fraction). The resulting wrapping amount $\eta N$ is used as a proxy for endocytosis capacity.

## Requirements
- Python 3.12+ (consistent with the research implementation)
- numpy >= 1.21.0
- matplotlib >= 3.4.0
- pandas >= 1.3.0

Install dependencies via pip:
```bash
pip install numpy matplotlib pandas

## How to Run
1. Clone this repository:
```bash
git clone https://github.com/[Your-Username]/[Repo-Name].git
cd [Repo-Name]
```
2. Execute the main simulation script:
```bash
python np_endocytosis_simulation.py
```
3. The script will:
Run numerical minimization of free energy for NP sizes from 10 nm to 100 nm diameter
Generate separate plots for CIE and CME pathways (y-axis tick labels hidden as per visualization requirements)
Print optimal NP diameter (max endocytosis capacity) for both pathways in the terminal

## Key Parameters
| Parameter       | Description                                  | Unit          | CIE Value | CME Value |
|-----------------|----------------------------------------------|---------------|-----------|-----------|
| $\kappa$ (kappa)| Membrane bending rigidity                    | $k_BT$        | 10        | 60        |
| $\mu$ (mu)      | Ligand-receptor binding energy gain          | $k_BT$        | 20        | 20        |
| $\xi_0$ (xi0)   | Initial receptor density on unbound membrane | Dimensionless | 0.05      | 0.05      |
| $A_0$ (A0)      | Area per receptor                            | $m^2$         | $(15e-9)^2$ | $(15e-9)^2$ |
| $M$             | Total membrane lattice points                | Dimensionless | $3.14e6$  | $3.14e6$  |
| $c$             | NP surface concentration                     | $1/A_0$       | 0.003     | 0.003     |

## Implementation Details

- A grid-search approach is used to numerically minimize the free energy over N (number of NPs) and η (wrapping fraction)
- Physical constraints are enforced (e.g., non-negative free membrane area, valid receptor density range)
- Separate simulations for CIE and CME pathways (differentiated by bending rigidity κ)
- Output plots show NP diameter (x-axis) vs. endocytosis capacity (ηN, y-axis) with optimal diameter marked

## Output

1. Interactive Plots:

  - Two separate figures for CIE and CME endocytosis
  - X-axis: NP diameter (nm)
  - Y-axis: Wrapped/endocytosed amount (ηN)
  - Vertical dashed line: Optimal NP diameter for maximum endocytosis capacity


2. Terminal Summary:

  - Optimal diameter and maximum wrapped amount for CIE/CME pathways
