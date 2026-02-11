#!/usr/bin/env python
"""
Infer age-dependent mRNA degradation rates using Bayesian inference.

This script implements a Bayesian model to infer degradation rates as a function
of MOLECULAR AGE (not spatial position or global time). The model calculates how
the probability of an individual mRNA molecule surviving depends on how long ago
it was transcribed.

Mathematical formulation:
    - D(age) = degradation rate as function of molecular age
    - S(age) = exp(-∫[0 to age] D(τ)dτ) = survival probability to age τ
    - m(T) = γ∫[0 to T] F(t)·S(T-t)dt = convolution of transcription and survival

Key biological insight: This models intrinsic molecular processes like poly-A tail
shortening, where each mRNA has an internal "timer" determining its stability,
independent of what's happening globally in the cell.

The model is fit using PyMC for MCMC sampling with NUTS.

Based on infer_D_across_stripe2.jl, updated for age-dependent degradation.
"""

import os
import sys
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from scipy.integrate import odeint, trapezoid
from scipy.interpolate import interp1d
from scipy.stats import invgamma
import pytensor
import pytensor.tensor as pt
from pytensor.graph import Op, Apply
import argparse


def smfish_ode(m, t, F_interp, D, gamma):
    """
    ODE model for mRNA dynamics with transcription input and degradation.
    
    dm/dt = γ*F(t) - D*m
    
    Args:
        m: mRNA concentration
        t: time
        F_interp: interpolated transcription function
        D: degradation rate
        gamma: transcription scaling factor
    
    Returns:
        dm/dt
    """
    F_t = F_interp(t)
    dmdt = gamma * F_t - D * m
    return dmdt


def integrate_ode_solution_numerical(D, gamma, F_interp, t_array):
    """
    Solve ODE numerically and return final mRNA concentration.
    
    Args:
        D: degradation rate (scalar)
        gamma: transcription scaling factor (scalar)
        F_interp: interpolated transcription function
        t_array: time points
    
    Returns:
        mRNA concentration at final time point
    """
    m0 = 0.0  # initial condition
    sol = odeint(smfish_ode, m0, t_array, args=(F_interp, D, gamma))
    return sol[-1, 0]


def solve_ode_analytical(D_age, gamma, F_values, t_array):
    """
    Analytical solution for AGE-DEPENDENT degradation.
    
    Calculates m(T) = gamma * Integral( F(t) * S(T-t) dt )
    
    Where S(age) is the survival probability for an mRNA molecule of given age:
    S(age) = exp(-∫[0 to age] D(τ) dτ)
    
    This is a CONVOLUTION between transcription history and survival probability,
    not an ODE. Each mRNA molecule "remembers" when it was made and degrades
    according to its individual age, not the global clock time.
    
    Args:
        D_age: Degradation rate as a function of AGE [n_timepoints]
        gamma: transcription scaling factor (scalar)
        F_values: transcription values at time points (array)
        t_array: time points (array)
    
    Returns:
        mRNA concentration at final time point
    """
    # Ensure inputs are numpy arrays for numerical computation
    D_age = np.asarray(D_age, dtype=np.float64)
    gamma = np.asarray(gamma, dtype=np.float64)
    F_values = np.asarray(F_values, dtype=np.float64)
    t_array = np.asarray(t_array, dtype=np.float64)
    
    # Clip D to reasonable range to avoid numerical issues
    # Lower upper bound since we want more stable mRNA
    D_age = np.clip(D_age, 1e-6, 1.0)
    
    # 1. Calculate Survival Curve S(tau)
    # Cumulative Hazard = Integral of D(age) from 0 to tau
    from scipy.integrate import cumulative_trapezoid
    cumulative_hazard = cumulative_trapezoid(D_age, t_array, initial=0)
    survival_prob = np.exp(-cumulative_hazard)
    
    # 2. Match Ages to Timepoints
    # At the final time T (end of experiment):
    # - The mRNA made at t=0 has age = T (needs survival_prob[-1])
    # - The mRNA made at t=T has age = 0 (needs survival_prob[0])
    # So we reverse the survival probability array to line up with F(t)
    survival_profile_reversed = survival_prob[::-1]
    
    # 3. Convolve (Integrate product)
    integrand = F_values * survival_profile_reversed
    integral = trapezoid(integrand, t_array)
    
    # Final solution: m(T) = γ * integral
    m_final = gamma * integral
    
    # Ensure output is valid
    if not np.isfinite(m_final):
        m_final = 0.0
    
    return float(m_final)


# ODESolverOp class removed - using vectorized PyMC operations instead
# This enables NUTS sampling with automatic differentiation


def calculate_expected_mRNA(D_age_array, gamma, F_data_arrays, t_array):
    """
    Calculate expected mRNA for all spatial bins using AGE-DEPENDENT degradation.
    
    Unlike spatial or temporal models, D here represents the degradation rate
    as a function of molecular age (how long since transcription), not position
    or global time. This models biological processes like poly-A tail shortening
    that act as molecular timers.
    
    Args:
        D_age_array: Array of degradation rates indexed by molecular age [n_timepoints]
        gamma: transcription scaling factor
        F_data_arrays: List of transcription data arrays for each bin [n_bins][n_traces_per_bin, n_timepoints]
        t_array: time points
    
    Returns:
        Array of expected mRNA concentrations [n_bins * n_traces_per_bin]
    """
    n_bins = len(F_data_arrays)
    expected_m_all = []
    
    for i in range(n_bins):
        # For each transcription trace, calculate expected mRNA
        # using the age-dependent degradation rate
        n_traces = F_data_arrays[i].shape[0]
        for j in range(n_traces):
            F_trace = F_data_arrays[i][j, :]
            # Pass the age-dependent rate vector to the convolution solver
            m_expected = solve_ode_analytical(D_age_array, gamma, F_trace, t_array)
            expected_m_all.append(m_expected)
    
    return np.array(expected_m_all)


def build_pymc_model(F_data_arrays, m_data, t_array, n_bins=5):
    """
    Build PyMC model for AGE-DEPENDENT degradation.
    
    D depends on the AGE of individual mRNA molecules (τ), not on spatial location
    or global clock time. This models biological processes like poly-A tail shortening
    where each molecule has an internal "timer" that determines its stability.
    
    Mathematical formulation:
        - D(age) = degradation rate as function of molecular age
        - S(age) = exp(-∫[0 to age] D(τ)dτ) = survival probability
        - m(T) = γ∫[0 to T] F(t)·S(T-t)dt = convolution of transcription and survival
    
    Key insight: An mRNA transcribed at t=5 min that is measured at t=15 min
    has age=10 min and uses D(10 min), regardless of what's happening globally at t=15.
    
    Prior distributions:
        log_D_age ~ GaussianRandomWalk with init_dist Normal(-2, 1)
            Enforces smooth age-dependency (e.g., gradual poly-A shortening)
        γ ~ HalfNormal(100.0) - transcription scaling
        σ ~ InverseGamma(2, 3) - observation noise
    
    Likelihood:
        m_obs ~ Normal(expected_m_all, σ²)
    
    Args:
        F_data_arrays: List of transcription data arrays [n_bins][n_traces_per_bin, n_timepoints]
        m_data: Observed mRNA data [n_bins * n_traces_per_bin] = 25 values
        t_array: time points
        n_bins: number of spatial bins
    
    Returns:
        PyMC model
    """
    print("Building PyMC model with Age-Dependent Degradation (Convolution)...")
    
    # Get structure info
    n_traces_per_bin = F_data_arrays[0].shape[0]
    n_timepoints = len(t_array)
    total_observations = n_bins * n_traces_per_bin
    
    print(f"  {n_bins} spatial bins (all share same D(age))")
    print(f"  {n_traces_per_bin} traces per bin")
    print(f"  {n_timepoints} age bins for D(age)")
    print(f"  {total_observations} total expected mRNA values")
    
    # Pre-calculate trapezoidal weights for integration
    dt = t_array[1] - t_array[0]
    weights = np.ones_like(t_array)
    weights[0] = 0.5
    weights[-1] = 0.5
    weights_scaled = pt.as_tensor_variable(weights * dt)
    
    with pm.Model() as model:
        # --- Priors ---
        
        # D_age represents the degradation rate at different AGES (0 min old, 1 min old...)
        # We use a Gaussian Random Walk to enforce that age-dependency is smooth.
        # (e.g., degradation might ramp up slowly as poly-A tails shorten)
        # Prior: log_D ~ N(-3.5, 0.5) → D ~ 0.03 min⁻¹ (half-life ~ 23 min)
        # This is more biologically realistic for mRNA stability
        log_D_age = pm.GaussianRandomWalk(
            'log_D_age', 
            sigma=0.1, 
            init_dist=pm.Normal.dist(-3.5, 0.5),
            shape=n_timepoints
        )
        D_age = pm.math.exp(log_D_age)
        pm.Deterministic('D', D_age)  # Track actual rates for saving and diagnostics
        
        # Scalar parameters
        # Increased gamma prior to allow higher transcription scaling
        gamma = pm.HalfNormal('gamma', sigma=500.0)
        sigma = pm.InverseGamma('sigma', alpha=2, beta=3)
        
        # --- Calculate Survival Curve S(tau) ---
        # 1. Cumulative Hazard = Integral of D(age) from 0 to tau
        cumulative_hazard = pm.math.cumsum(D_age) * dt
        
        # 2. Survival Probability S(tau) = exp( - Cumulative Hazard )
        # This vector describes: [Prob surviving 0 min, Prob surviving 1 min, ...]
        survival_prob = pm.math.exp(-cumulative_hazard)
        
        # --- The Convolution (History Integral) ---
        # For an observation at time T:
        # - mRNA produced at t=0 has age = T (needs survival_prob[-1])
        # - mRNA produced at t=T has age = 0 (needs survival_prob[0])
        # So we reverse the survival vector to match the timepoints of F.
        
        # Reverse the survival profile to align ages with transcription times
        # PyTensor doesn't support [::-1], so we use explicit indexing
        survival_profile_reversed = survival_prob[::-1]
        
        # --- Vectorized Convolution ---
        expected_m_list = []
        
        for i in range(n_bins):
            # Loop over traces in this bin
            n_traces = F_data_arrays[i].shape[0]
            for j in range(n_traces):
                F_trace = F_data_arrays[i][j, :]
                
                # Contribution = Transcription * Probability of Surviving until End
                # This is the convolution: each F(t) is weighted by S(T-t)
                integrand = F_trace * survival_profile_reversed
                
                # Integrate using trapezoidal rule
                integral = pm.math.sum(integrand * weights_scaled)
                
                # Final result: gamma * integral
                m_ij = gamma * integral
                expected_m_list.append(m_ij)
        
        # Stack results into a single tensor
        expected_m = pm.Deterministic(
            'expected_m',
            pt.stack(expected_m_list)
        )
        
        # Likelihood
        m_obs = pm.Normal(
            'm_obs', 
            mu=expected_m, 
            sigma=sigma, 
            observed=m_data
        )
    
    return model


def build_pymc_model_with_polya_protection(F_data_arrays, m_data, t_array, n_bins=5):
    """
    Build PyMC model with mechanistic poly-A tail protection.
    Details pulled from: https://pmc.ncbi.nlm.nih.gov/articles/PMC11649921/
    
    This model explicitly represents the biological mechanism of deadenylation
    and Pab1-mediated protection, creating a biphasic degradation pattern:
    1. Slow deadenylation phase while poly-A tail is long (>20 As, Pab1-protected)
    2. Rapid decay phase after decapping when poly-A tail is short (<20 As)
    
    Biological model:
        - NA(age) = NA_0 - deadenylation_rate × age  (poly-A tail shortens linearly)
        - protection(NA) = tanh(β × NA)  (Pab1 binding strength, saturates at ~20 As)
        - D(age) = D_protected + D_unprotected × (1 - protection(NA(age)))
    
    This creates the observed biphasic pattern from the literature:
    - Low D while protected (slow deadenylation)
    - Sharp transition as protection drops
    - High D after decapping (fast Xrn1-mediated decay)
    
    Prior distributions:
        NA_0 ~ Normal(50, 15) - Initial poly-A tail length (adjusted for short-lived transcripts)
        deadenylation_rate ~ HalfNormal(10.0) - Adenosines removed per minute (faster for eve)
        β ~ Normal(0.096, 0.03) - Protection sharpness (from literature)
        D_protected ~ HalfNormal(0.1) - Slow decay rate while protected
        D_unprotected ~ HalfNormal(1.0) - Fast decay rate after decapping
        γ ~ HalfNormal(500.0) - Transcription scaling
        σ ~ InverseGamma(2, 3) - Observation noise
    
    Args:
        F_data_arrays: List of transcription data arrays [n_bins][n_traces_per_bin, n_timepoints]
        m_data: Observed mRNA data [n_bins * n_traces_per_bin]
        t_array: time points (molecular ages)
        n_bins: number of spatial bins
    
    Returns:
        PyMC model
    """
    print("Building PyMC model with Mechanistic Poly-A Protection...")
    
    n_traces_per_bin = F_data_arrays[0].shape[0]
    n_timepoints = len(t_array)
    total_observations = n_bins * n_traces_per_bin
    
    print(f"  {n_bins} spatial bins")
    print(f"  {n_traces_per_bin} traces per bin")
    print(f"  {n_timepoints} time points")
    print(f"  {total_observations} total observations")
    print("  Model: Biphasic degradation via poly-A tail dynamics")
    
    dt = t_array[1] - t_array[0]
    weights = np.ones_like(t_array)
    weights[0] = 0.5
    weights[-1] = 0.5
    weights_scaled = pt.as_tensor_variable(weights * dt)
    
    # Convert time array to PyTensor
    t_tensor = pt.as_tensor_variable(t_array)
    
    with pm.Model() as model:
        # --- Priors for Poly-A Tail Dynamics ---
        
        # Initial poly-A tail length
        # For short-lived transcripts like eve (half-life ~7 min), may start with shorter tails
        # or have faster deadenylation to reach critical threshold quickly
        NA_0 = pm.Normal('NA_0', mu=60, sigma=10)
        
        # Deadenylation rate (adenosines removed per minute)
        # For transition at ~4 min: need to remove ~30 As in 4 min → ~7.5 As/min
        # Using broader prior to let data inform the rate
        deadenylation_rate = pm.HalfNormal('deadenylation_rate', sigma=10.0)
        
        # Protection parameter (from literature: β = 0.096)
        # This controls how sharply protection drops below ~20 As
        # TIGHTENED PRIOR: Previous loose prior (sigma=0.03) led to bimodality
        # and non-convergence. We trust the biochemical literature here.
        beta = pm.Normal('beta', mu=0.096, sigma=0.005)
        
        # Base degradation rate (slow decay while Pab1-protected)
        # For eve transcript with 7-min half-life, even "slow" decay needs to be substantial
        D_protected = pm.HalfNormal('D_protected', sigma=0.1)
        
        # Fast degradation rate (after decapping when poly-A is short)
        # Should be much higher - represents Xrn1-mediated decay
        # For 7-min half-life, total decay needs: ln(2)/7 ≈ 0.1 min⁻¹
        D_unprotected = pm.HalfNormal('D_unprotected', sigma=1.0)
        
        # Transcription and noise parameters
        gamma = pm.HalfNormal('gamma', sigma=500.0)
        sigma = pm.InverseGamma('sigma', alpha=2, beta=3)
        
        # --- Calculate Age-Dependent Degradation via Poly-A Mechanism ---
        
        # 1. Poly-A tail length as function of age
        # NA(age) = NA_0 - deadenylation_rate × age
        # Clip to minimum of 0 adenosines
        NA_age = pm.math.maximum(0, NA_0 - deadenylation_rate * t_tensor)
        
        # 2. Protection factor (Pab1 binding strength)
        # tanh(β × NA) saturates at ~20 As, drops sharply below
        # This models the modified gamma distribution from the literature
        protection_factor = pm.math.tanh(beta * NA_age)
        
        # 3. Age-dependent degradation rate
        # When protection is high (1.0): D ≈ D_protected (slow deadenylation)
        # When protection is low (0.0): D ≈ D_protected + D_unprotected (fast decay)
        D_age = D_protected + D_unprotected * (1 - protection_factor)
        
        # Store for diagnostics and visualization
        pm.Deterministic('D', D_age)
        pm.Deterministic('NA', NA_age)
        pm.Deterministic('protection', protection_factor)
        
        # --- Calculate Survival Probability ---
        cumulative_hazard = pm.math.cumsum(D_age) * dt
        survival_prob = pm.math.exp(-cumulative_hazard)
        survival_profile_reversed = survival_prob[::-1]
        
        # --- Vectorized Convolution ---
        expected_m_list = []
        
        for i in range(n_bins):
            n_traces = F_data_arrays[i].shape[0]
            for j in range(n_traces):
                F_trace = F_data_arrays[i][j, :]
                integrand = F_trace * survival_profile_reversed
                integral = pm.math.sum(integrand * weights_scaled)
                m_ij = gamma * integral
                expected_m_list.append(m_ij)
        
        expected_m = pm.Deterministic('expected_m', pt.stack(expected_m_list))
        
        # Likelihood
        m_obs = pm.Normal('m_obs', mu=expected_m, sigma=sigma, observed=m_data)
    
    return model


def run_mcmc_inference(model, n_samples=10000, n_chains=4, target_accept=0.99):
    """
    Run MCMC sampling using NUTS algorithm.
    
    Args:
        model: PyMC model
        n_samples: number of samples per chain
        n_chains: number of parallel chains
        target_accept: target acceptance rate for NUTS
    
    Returns:
        InferenceData object with samples
    """
    print(f"Running MCMC inference: {n_samples} samples × {n_chains} chains...")
    
    with model:
        # Sample using NUTS (No-U-Turn Sampler)
        trace = pm.sample(
            draws=n_samples,
            chains=n_chains,
            cores=n_chains,
            target_accept=target_accept,
            return_inferencedata=True,
            random_seed=14
        )
    
    return trace


def save_chain_results(trace, output_path):
    """
    Save MCMC chain samples to CSV.
    
    Now handles temporal degradation rates D(t) with shape (chains, draws, n_timepoints)
    instead of spatial degradation rates.
    
    Args:
        trace: InferenceData object
        output_path: path to save CSV
    """
    print(f"Saving chain results to: {output_path}")
    
    posterior = trace.posterior
    samples_dict = {}
    
    # Handle temporal degradation - save 'D' (the deterministic with actual rates)
    # If 'D' is not available, try 'log_D' and transform
    if 'D' in posterior:
        var_data = posterior['D'].values
        # D has shape (chains, draws, n_timepoints)
        n_chains, n_draws, n_timepoints = var_data.shape
        flattened = var_data.reshape(-1, n_timepoints)
        for i in range(n_timepoints):
            samples_dict[f'D[{i}]'] = flattened[:, i]
    elif 'log_D' in posterior:
        var_data = posterior['log_D'].values
        n_chains, n_draws, n_timepoints = var_data.shape
        flattened = var_data.reshape(-1, n_timepoints)
        # Transform to actual scale
        flattened_D = np.exp(flattened)
        for i in range(n_timepoints):
            samples_dict[f'D[{i}]'] = flattened_D[:, i]
    
    # Handle poly-A tail dynamics (if present)
    for array_var in ['NA', 'protection']:
        if array_var in posterior:
            var_data = posterior[array_var].values
            n_chains, n_draws, n_timepoints = var_data.shape
            flattened = var_data.reshape(-1, n_timepoints)
            for i in range(n_timepoints):
                samples_dict[f'{array_var}[{i}]'] = flattened[:, i]
    
    # Handle scalar parameters (standard and poly-A specific)
    scalar_params = ['gamma', 'sigma', 'NA_0', 'deadenylation_rate', 'beta', 'D_protected', 'D_unprotected']
    for var_name in scalar_params:
        if var_name in posterior:
            var_data = posterior[var_name].values
            samples_dict[var_name] = var_data.flatten()
    
    df = pd.DataFrame(samples_dict)
    
    # Save to CSV
    df.to_csv(output_path, index=False)
    
    print(f"Saved {len(df)} samples")


def plot_trace(trace, output_path):
    """Plot MCMC trace plots for diagnostics."""
    print(f"Generating trace plot: {output_path}")
    
    # Check which variables are present
    posterior = trace.posterior
    var_names = ['gamma', 'sigma']
    
    # Add poly-A specific parameters if present (but not D, NA, protection - too many dimensions)
    polya_params = ['NA_0', 'deadenylation_rate', 'beta', 'D_protected', 'D_unprotected']
    for param in polya_params:
        if param in posterior:
            var_names.append(param)
    
    az.plot_trace(trace, var_names=var_names)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def print_summary_statistics(trace):
    """Print summary statistics of posterior distributions with age-dependent degradation rates."""
    print("\n=== Posterior Summary Statistics ===")
    
    # Check which variables are present
    posterior = trace.posterior
    var_names = ['D', 'gamma', 'sigma']
    
    # Add poly-A specific parameters if present
    polya_params = ['NA_0', 'deadenylation_rate', 'beta', 'D_protected', 'D_unprotected']
    for param in polya_params:
        if param in posterior:
            var_names.append(param)
    
    # Add log_D_age if present (smooth model)
    if 'log_D_age' in posterior:
        var_names.insert(0, 'log_D_age')
    
    summary = az.summary(trace, var_names=var_names)
    print(summary)
    
    # Print poly-A specific summary if parameters are present
    if 'NA_0' in posterior:
        print("\n=== Poly-A Tail Model Parameters ===")
        NA_0_mean = posterior['NA_0'].values.mean()
        deadenyl_mean = posterior['deadenylation_rate'].values.mean()
        beta_mean = posterior['beta'].values.mean()
        D_prot_mean = posterior['D_protected'].values.mean()
        D_unprot_mean = posterior['D_unprotected'].values.mean()
        
        print(f"Initial poly-A length (NA_0):        {NA_0_mean:.1f} adenosines")
        print(f"Deadenylation rate:                  {deadenyl_mean:.3f} As/min")
        print(f"Protection sharpness (β):            {beta_mean:.4f}")
        print(f"Protected degradation rate:          {D_prot_mean:.4f} min⁻¹ (t_1/2 = {np.log(2)/D_prot_mean:.1f} min)")
        print(f"Unprotected degradation rate:        {D_unprot_mean:.4f} min⁻¹ (contributes to fast decay)")
        print(f"Total fast decay rate:               ~{D_prot_mean + D_unprot_mean:.4f} min⁻¹ (t_1/2 = {np.log(2)/(D_prot_mean + D_unprot_mean):.1f} min)")
        
        # Calculate when poly-A drops below 20 As
        if deadenyl_mean > 0:
            time_to_20As = max(0, (NA_0_mean - 20) / deadenyl_mean)
            print(f"\nTime until poly-A < 20 As:           ~{time_to_20As:.1f} min")
            print("(This is when protection drops and fast decay begins)")
    
    # Get D values to compute half-lives as function of age
    posterior = trace.posterior
    if 'D' in posterior:
        D_data = posterior['D'].values  # shape: (chains, draws, n_ages)
    elif 'log_D_age' in posterior:
        log_D_data = posterior['log_D_age'].values
        D_data = np.exp(log_D_data)
    else:
        print("No degradation rate variable found in trace")
        return
    
    # Average across chains and draws
    D_mean = np.mean(D_data, axis=(0, 1))
    D_std = np.std(np.mean(D_data, axis=1), axis=0)  # std of chain means
    
    print("\n=== Age-Dependent Degradation Rate D(age) ===")
    print("Degradation rate depends on molecular age, not spatial position or global time:")
    print(f"{'Age (min)':<12} {'D (min⁻¹)':<15} {'t_1/2 (min)':<15} {'Std Dev':<12}")
    print("-" * 54)
    for age_idx, D_val in enumerate(D_mean):
        halflife = np.log(2) / D_val
        # age_idx * time_step gives the actual age in minutes
        print(f"{age_idx:<12} {D_val:<15.3f} {halflife:<15.2f} {D_std[age_idx]:<12.3f}")


def load_data(transcription_path, mrna_path):
    """
    Load transcription and mRNA data.
    
    Args:
        transcription_path: Path to transcription data CSV (no headers)
        mrna_path: Path to mRNA data CSV (no headers)
    
    Returns:
        F_data: transcription data array [n_bins * n_traces_per_bin, n_timepoints]
        m_data: mRNA data array [n_bins]
    """
    print(f"Loading transcription data: {transcription_path}")
    F_data = pd.read_csv(transcription_path, header=None).values
    
    print(f"Loading mRNA data: {mrna_path}")
    m_data = pd.read_csv(mrna_path, header=None).values.flatten()
    
    print(f"Transcription data shape: {F_data.shape}")
    print(f"mRNA data shape: {m_data.shape}")
    
    return F_data, m_data


def organize_transcription_data(F_data, n_ap_bins=5, n_dv_bins=5):
    """
    Organize transcription data by AP position.
    
    The preprocessing creates a 5×5 grid where rows are ordered by spatial position.
    We group by AP position to infer AP-specific degradation rates.
    
    Assuming row order from preprocessing is:
    - Rows 0-4: AP position 1, DV positions 1-5
    - Rows 5-9: AP position 2, DV positions 1-5
    - Rows 10-14: AP position 3, DV positions 1-5
    - Rows 15-19: AP position 4, DV positions 1-5
    - Rows 20-24: AP position 5, DV positions 1-5
    
    Args:
        F_data: Transcription data [25 spatial bins, n_timepoints]
        n_ap_bins: Number of AP bins (default: 5)
        n_dv_bins: Number of DV bins per AP position (default: 5)
    
    Returns:
        List of arrays [n_ap_bins][n_dv_bins, n_timepoints]
        Each array contains the 5 DV traces for one AP position
    """
    total_bins = n_ap_bins * n_dv_bins
    
    if len(F_data) != total_bins:
        raise ValueError(
            f"Expected {total_bins} bins ({n_ap_bins} AP × {n_dv_bins} DV), "
            f"got {len(F_data)} rows"
        )
    
    F_data_arrays = []
    for ap_idx in range(n_ap_bins):
        start_idx = ap_idx * n_dv_bins
        end_idx = (ap_idx + 1) * n_dv_bins
        F_data_arrays.append(F_data[start_idx:end_idx, :])
    
    print(f"Organized transcription data:")
    print(f"  {n_ap_bins} AP positions")
    print(f"  {n_dv_bins} DV bins per AP position")
    print(f"  Total: {total_bins} spatial bins")
    print(f"  Each bin has one averaged transcription trace")
    
    return F_data_arrays


def organize_mrna_data(m_data, n_ap_bins=5, n_dv_bins=5):
    """
    Validate mRNA data matches expected structure from preprocessing.
    
    The smFISH preprocessing produces exactly 25 values:
    - 5 AP bins × 5 DV bins = 25 spatial bins
    - Each value = average mRNA/nucleus in that spatial bin
    - Ordered same as transcription data (by AP then DV)
    
    Args:
        m_data: mRNA data [25] - one value per spatial bin
        n_ap_bins: Number of AP bins (default: 5)
        n_dv_bins: Number of DV bins per AP position (default: 5)
    
    Returns:
        mRNA data unchanged [25]
    
    Raises:
        ValueError if data size doesn't match 25 bins
    """
    total_bins = n_ap_bins * n_dv_bins
    
    if len(m_data) != total_bins:
        raise ValueError(
            f"Expected {total_bins} mRNA values ({n_ap_bins} AP × {n_dv_bins} DV), "
            f"got {len(m_data)}. Check that smFISH preprocessing uses same 5×5 binning "
            f"as transcription preprocessing."
        )
    
    print(f"mRNA data validated:")
    print(f"  {len(m_data)} values (one per spatial bin)")
    print(f"  {n_ap_bins} AP positions × {n_dv_bins} DV bins")
    print(f"  Each value = average mRNA/nucleus in that bin")
    
    return m_data


def main():
    parser = argparse.ArgumentParser(
        description='Infer spatially-varying mRNA degradation rates'
    )
    parser.add_argument( '--transcription', required=True, help='Path to transcription data CSV (no headers)')
    parser.add_argument( '--mrna', required=True, help='Path to mRNA count data CSV (no headers)')
    parser.add_argument( '--output-chain', required=True, help='Output path for MCMC chain CSV')
    parser.add_argument( '--output-trace', required=True, help='Output path for trace plot')
    parser.add_argument( '--n-samples', type=int, default=10000, help='Number of MCMC samples per chain (default: 10000)')
    parser.add_argument( '--n-chains', type=int, default=4, help='Number of parallel MCMC chains (default: 4)')
    parser.add_argument( '--n-ap-bins', type=int, default=5, help='Number of anterior-posterior bins (default: 5)')
    parser.add_argument( '--n-dv-bins', type=int, default=5, help='Number of dorsal-ventral bins (default: 5)')

    args = parser.parse_args()
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Load data
    F_data, m_data_raw = load_data(args.transcription, args.mrna)
    
    # Organize data by AP position
    # Structure: 5 AP bins, each containing 5 DV transcription traces and 5 DV mRNA values
    # Model infers 5 degradation rates (one per AP position)
    # Each D[i] applies to all 5 DV bins at that AP position
    F_data_arrays = organize_transcription_data(F_data, args.n_ap_bins, args.n_dv_bins)
    m_data = organize_mrna_data(m_data_raw, args.n_ap_bins, args.n_dv_bins)
    
    # Time array (0 to 20 minutes)
    n_timepoints = F_data.shape[1]
    t_array = np.arange(0, n_timepoints * 20, 20) / 60.0
    
    # Build Bayesian model with mechanistic poly-A protection
    # This creates a biphasic degradation pattern:
    # - Slow decay while Pab1-protected (poly-A > 20 As)
    # - Fast decay after decapping (poly-A < 20 As)
    model = build_pymc_model_with_polya_protection(F_data_arrays, m_data, t_array, args.n_ap_bins)
    
    # Run MCMC inference
    trace = run_mcmc_inference(model, args.n_samples, args.n_chains)
    
    # Save results
    save_chain_results(trace, args.output_chain)
    plot_trace(trace, args.output_trace)
    
    # Print summary
    print_summary_statistics(trace)
    
    print("\n=== Inference complete ===")


if __name__ == '__main__':
    main()
