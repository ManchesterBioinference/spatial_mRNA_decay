import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

def solve_analytical(D, gamma, F_values, t_array):
    """Re-implementation of the model logic for validation."""
    t_final = t_array[-1]
    dt = t_array[1] - t_array[0]
    weights = np.ones_like(t_array)
    weights[0], weights[-1] = 0.5, 0.5
    
    exp_term = np.exp(D * (t_array - t_final))
    integrand = F_values * exp_term
    integral = np.sum(integrand * weights * dt)
    return gamma * integral

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chain', required=True, help='MCMC chain CSV')
    parser.add_argument('--transcription', required=True, help='Original transcription CSV')
    parser.add_argument('--mrna', required=True, help='Original mRNA CSV')
    parser.add_argument('--max_time', type=int, required=True, help='Maximum time in seconds')
    parser.add_argument('--output', default='posterior_predictive_check.png')
    args = parser.parse_args()

    # 1. Load Data
    df_samples = pd.read_csv(args.chain)
    F_data = pd.read_csv(args.transcription, header=None).values
    m_obs = pd.read_csv(args.mrna, header=None).values.flatten()
    t_array = np.arange(0, args.max_time + 1, 20) / 60.0

    # 2. Get Posterior Means
    # We use the mean of the samples as our "best fit" estimate
    D_means = [df_samples[f'D[{i}]'].mean() for i in range(5)]
    gamma_mean = df_samples['gamma'].mean()

    # 3. Predict mRNA for every trace
    m_pred = []
    for i in range(5): # 5 AP Bins
        for j in range(5): # 5 traces per bin
            idx = i * 5 + j
            F_trace = F_data[idx, :]
            val = solve_analytical(D_means[i], gamma_mean, F_trace, t_array)
            m_pred.append(val)

    # 4. Plotting
    plt.figure(figsize=(10, 6))
    
    # Plot Identity Line (Perfect Fit)
    max_val = max(max(m_obs), max(m_pred))
    plt.plot([0, max_val], [0, max_val], color='red', linestyle='--', alpha=0.5, label='Perfect Fit')
    
    # Plot Bins with different colors
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for i in range(5):
        start, end = i*5, (i+1)*5
        plt.scatter(m_obs[start:end], m_pred[start:end], 
                    color=colors[i], label=f'Bin {i+1} (D={D_means[i]:.2f})', s=100, edgecolors='k')

    plt.xlabel('Observed mRNA (Data)')
    plt.ylabel('Predicted mRNA (Model)')
    plt.title('Posterior Predictive Check: Model vs. Data')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(args.output, dpi=300)
    print(f"Validation plot saved to {args.output}")

if __name__ == "__main__":
    main()