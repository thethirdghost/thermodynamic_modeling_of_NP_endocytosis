import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ==========================
# Core Parameters (Energy unit: kBT)
# ==========================
# Define two sets of parameters (CIE and CME)
PARAM_SETS = {
    "CIE": {
        "kappa": 10.0,       # Membrane bending rigidity (kBT)
        "mu": 20.0,          # Ligand-receptor binding energy (kBT)
        "xi0": 0.05,         # Initial receptor density (dimensionless)
        "A0": (15e-9)**2,    # Area per receptor (m^2)
        "M": 3.14e6,         # Total membrane lattice points (dimensionless, unit: A0)
        "c": 0.003           # NP surface concentration (unit: 1/A0)
    },
    "CME": {
        "kappa": 60.0,       # Membrane bending rigidity (kBT)
        "mu": 20.0,          # Ligand-receptor binding energy (kBT)
        "xi0": 0.05,         # Initial receptor density (dimensionless)
        "A0": (15e-9)**2,    # Area per receptor (m^2)
        "M": 3.14e6,         # Total membrane lattice points (dimensionless, unit: A0)
        "c": 0.003           # NP surface concentration (unit: 1/A0)
    }
}

# Simulation range settings
R_LIST = np.linspace(10e-9, 100e-9, 120)  # Radius range (m)
ETA_GRID = np.linspace(0.0, 1.0, 101)     # Eta grid
ETA_SUB = ETA_GRID[::5]                   # Subsampled eta grid
N_SAMPLES = 200                           # Number of N sampling points
EPS = 1e-12                               # Small value to avoid log(0)

# ==========================
# Core Functions
# ==========================
def mix_entropy(x):
    """
    Calculate mixing entropy term: x ln x + (1-x) ln(1-x)
    Args:
        x: Receptor density (dimensionless)
    Returns:
        Entropy value (dimensionless)
    """
    x = np.clip(x, EPS, 1 - EPS)  # Avoid log(0) or log(1)
    return x * np.log(x) + (1 - x) * np.log(1 - x)

def W_vectorized(R, N_grid, eta, params):
    """
    Vectorized calculation of free energy W(N; R, eta) (unit: kBT)
    Args:
        R: NP radius (m)
        N_grid: Array of N values to sample
        eta: Wrapping fraction (0 ≤ eta ≤ 1)
        params: Dictionary of simulation parameters
    Returns:
        Array of W values corresponding to N_grid
    """
    # Calculate surface area factor K (dimensionless)
    K = (4 * np.pi * R**2) / params["A0"]
    
    N = N_grid.astype(float)
    Mb = eta * N * K       # Membrane area bound to NP
    Mf = params["M"] - Mb  # Free membrane area
    Lb = eta * N * K       # Number of bound ligands/receptors
    
    total_receptors = params["xi0"] * params["M"]
    W = np.full_like(N, np.inf, dtype=float)  # Initialize with infinity
    
    # Filter valid states (physical constraints)
    valid = (Mf > 0) & (Lb >= 0) & (Lb < total_receptors)
    if not np.any(valid):
        return W
    
    # Calculate free receptor density
    xi_f = (total_receptors - Lb[valid]) / Mf[valid]
    valid2 = (xi_f > 0) & (xi_f < 1)
    if not np.any(valid2):
        return W
    
    # Extract valid subsets for calculation
    Mf2 = Mf[valid][valid2]
    Lb2 = Lb[valid][valid2]
    N2 = N[valid][valid2]
    xi2 = xi_f[valid2]
    
    # Calculate energy components
    W_entropy = Mf2 * mix_entropy(xi2)       # Entropy term
    W_binding = -params["mu"] * Lb2          # Binding energy term
    W_bending = 8 * np.pi * params["kappa"] * eta * N2  # Bending energy term
    
    # Total free energy
    W_valid = W_entropy + W_binding + W_bending
    
    # Assign valid values back to W array
    idx_valid = np.where(valid)[0]
    idx_valid2 = idx_valid[valid2]
    W[idx_valid2] = W_valid
    
    return W

def run_simulation(param_set_name):
    """
    Run simulation for a given parameter set
    Args:
        param_set_name: Name of parameter set ("CIE" or "CME")
    Returns:
        DataFrame containing simulation results
    """
    params = PARAM_SETS[param_set_name]
    records = []
    
    for R in R_LIST:
        K = (4 * np.pi * R**2) / params["A0"]
        
        # Calculate maximum possible N (geometric and receptor constraints)
        N_geo_max = int(params["c"] * params["M"])
        N_rec_max = int((params["xi0"] * params["M"]) / (K + 1e-30))
        N_max = max(1, min(N_geo_max, N_rec_max))
        
        # Create N sampling grid (unique integer values)
        N_grid = np.linspace(0, N_max, N_SAMPLES).astype(int)
        N_grid = np.unique(N_grid)
        
        # Find optimal N and eta that minimize W
        best_W = np.inf
        best_N = 0
        best_eta = 0.0
        
        for eta in ETA_SUB:
            W_vals = W_vectorized(R, N_grid, eta, params)
            idx = np.argmin(W_vals).astype(int)
            if W_vals[idx] < best_W:
                best_W = float(W_vals[idx])
                best_N = int(N_grid[idx])
                best_eta = float(eta)
        
        # Record results
        records.append({
            "param_set": param_set_name,
            "R_nm": R * 1e9,
            "D_nm": 2.0 * R * 1e9,
            "N_opt": best_N,
            "eta_opt": best_eta,
            "wrapped_etaN": best_eta * best_N,
            "W_min": best_W
        })
    
    return pd.DataFrame(records)

def plot_simulation_results(df, param_set_name, D_opt, max_wrapped):
    """
    Plot simulation results for a single parameter set (separate figure)
    Args:
        df: DataFrame containing simulation results
        param_set_name: Name of parameter set ("CIE" or "CME")
        D_opt: Optimal diameter (nm)
        max_wrapped: Maximum wrapped amount (ηN)
    """
    # Create a new figure for each parameter set
    plt.figure(figsize=(12, 8))
    
    # Plot wrapped amount vs diameter
    plt.plot(df["D_nm"], df["wrapped_etaN"], linewidth=2, 
             label=r"Wrapped amount $\eta N$")
    # Add vertical line for optimal diameter
    plt.axvline(D_opt, linestyle="--", alpha=0.7,
                label=f"$D_{{opt}}$ ≈ {D_opt:.1f} nm")
    
    # Plot formatting
    plt.xlabel("NP diameter (nm)")
    plt.ylabel(r"Wrapped / endocytosed amount ($\eta N$)")
    plt.title(f"Endocytosis with Partial Wrapping - {param_set_name} (Minimize W over N and η)")
    plt.legend()
    plt.grid(alpha=0.3)
    
    # Hide y-axis tick labels (keep y-axis label)
    plt.gca().set_yticklabels([])
    
    plt.tight_layout()
    # Save each plot as separate PDF (uncomment if needed)
    # plt.savefig(f"wrapped_amount_vs_diameter_{param_set_name}.pdf")
    plt.show()

# ==========================
# Run Simulations and Plot Results
# ==========================
# Run for both CIE and CME parameter sets
df_cie = run_simulation("CIE")
df_cme = run_simulation("CME")

# Calculate Optimal Diameters and plot separately
results_summary = {}
for param_set in ["CIE", "CME"]:
    # Get corresponding DataFrame
    df = df_cie if param_set == "CIE" else df_cme
    D_nm = df["D_nm"].values
    wrapped_etaN = df["wrapped_etaN"].values
    
    # Calculate optimal diameter
    valid = wrapped_etaN > 0
    D_opt = D_nm[np.argmax(wrapped_etaN)] if valid.any() else 0.0
    max_wrapped = np.max(wrapped_etaN) if valid.any() else 0.0
    
    results_summary[param_set] = {
        "D_opt": D_opt,
        "max_wrapped": max_wrapped
    }
    
    # Plot single figure for this parameter set
    plot_simulation_results(df, param_set, D_opt, max_wrapped)

# Print summary of optimal diameters
print("\n=== Simulation Summary ===")
for param_set, stats in results_summary.items():
    print(f"{param_set}: Optimal diameter = {stats['D_opt']:.1f} nm, Max wrapped amount = {stats['max_wrapped']:.2f}")
