import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.integrate import cumulative_trapezoid, trapezoid

def solve_analytical(D_age, gamma, F_values, t_array):
    """
    Re-implementation of age-dependent degradation model for validation.
    
    Uses convolution with survival probability:
    m(T) = gamma * Integral( F(t) * S(T-t) dt )
    where S(age) = exp(-Integral D(tau) dtau)
    """
    # Ensure arrays
    D_age = np.asarray(D_age, dtype=np.float64)
    gamma = float(gamma)
    F_values = np.asarray(F_values, dtype=np.float64)
    
    # Calculate survival curve
    cumulative_hazard = cumulative_trapezoid(D_age, t_array, initial=0)
    survival_prob = np.exp(-cumulative_hazard)
    
    # Reverse to align with transcription times
    # mRNA made at t=0 has age=T, mRNA made at t=T has age=0
    survival_profile_reversed = survival_prob[::-1]
    
    # Convolve
    integrand = F_values * survival_profile_reversed
    integral = trapezoid(integrand, t_array)
    
    return gamma * integral

def main():
    parser = argparse.ArgumentParser(
        description='Validate age-dependent degradation MCMC results with posterior predictive check'
    )
    parser.add_argument('--chain', required=True, help='MCMC chain CSV')
    parser.add_argument('--transcription', required=True, help='Original transcription CSV')
    parser.add_argument('--mrna', required=True, help='Original mRNA CSV')
    parser.add_argument('--output', default='posterior_predictive_check.png')
    args = parser.parse_args()

    # 1. Load Data
    df_samples = pd.read_csv(args.chain)
    F_data = pd.read_csv(args.transcription, header=None).values
    m_obs = pd.read_csv(args.mrna, header=None).values.flatten()
    t_array = np.arange(0, 1201, 20) / 60.0

    # 2. Extract D_age posterior means
    # For age-dependent model, D columns represent degradation at different molecular ages
    D_columns = [col for col in df_samples.columns if col.startswith('D[') and col.endswith(']')]
    D_columns_sorted = sorted(D_columns, key=lambda x: int(x.split('[')[1].split(']')[0]))
    
    n_ages = len(D_columns_sorted)
    print(f"Found {n_ages} age bins in chain")
    
    # Get posterior mean for each age
    D_age_means = np.array([df_samples[col].mean() for col in D_columns_sorted])
    gamma_mean = df_samples['gamma'].mean()
    
    print(f"Mean gamma: {gamma_mean:.2f}")
    print(f"Mean D(age=0): {D_age_means[0]:.3f} min⁻¹")
    print(f"Mean D(age=max): {D_age_means[-1]:.3f} min⁻¹")

    # 3. Predict mRNA for every trace
    # All traces use the same age-dependent degradation curve
    m_pred = []
    for i in range(5):  # 5 spatial bins
        for j in range(5):  # 5 traces per bin
            idx = i * 5 + j
            F_trace = F_data[idx, :]
            # Use the SAME D_age curve for all traces
            val = solve_analytical(D_age_means, gamma_mean, F_trace, t_array)
            m_pred.append(val)
    
    m_pred = np.array(m_pred)

    # 4. Calculate fit statistics
    residuals = m_obs - m_pred
    rmse = np.sqrt(np.mean(residuals**2))
    r_squared = 1 - (np.sum(residuals**2) / np.sum((m_obs - np.mean(m_obs))**2))
    
    print(f"\nModel fit statistics:")
    print(f"  RMSE: {rmse:.2f}")
    print(f"  R²: {r_squared:.3f}")

    # 5. Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left panel: Posterior Predictive Check
    max_val = max(max(m_obs), max(m_pred))
    ax1.plot([0, max_val], [0, max_val], color='red', linestyle='--', alpha=0.5, 
             linewidth=2, label='Perfect Fit')
    
    # Plot all bins (they share the same D_age, so color by spatial bin for reference)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for i in range(5):
        start, end = i*5, (i+1)*5
        ax1.scatter(m_obs[start:end], m_pred[start:end], 
                    color=colors[i], label=f'Spatial Bin {i+1}', 
                    s=100, edgecolors='k', alpha=0.7)

    ax1.set_xlabel('Observed mRNA (Data)', fontsize=12)
    ax1.set_ylabel('Predicted mRNA (Model)', fontsize=12)
    ax1.set_title(f'Posterior Predictive Check\nR² = {r_squared:.3f}, RMSE = {rmse:.2f}', 
                  fontsize=13)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # Right panel: Residuals
    bin_indices = np.repeat(np.arange(5), 5)
    ax2.scatter(m_pred, residuals, c=bin_indices, cmap='viridis', 
                s=100, edgecolors='k', alpha=0.7)
    ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=2)
    ax2.set_xlabel('Predicted mRNA', fontsize=12)
    ax2.set_ylabel('Residuals (Observed - Predicted)', fontsize=12)
    ax2.set_title('Residual Plot', fontsize=13)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\nValidation plot saved to {args.output}")

if __name__ == "__main__":
    main()