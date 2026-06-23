#!/usr/bin/env python
"""
Tests for null constant degradation rate inference script.

Test strategy:
1. Test with synthetic data where true degradation rate is known
2. Verify output format (CSV columns, trace plot)
3. Verify model can recover reasonable parameter estimates
"""

import os
import sys
import numpy as np
import pandas as pd
import pytest
import tempfile
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))


def create_synthetic_data(n_bins=5, n_traces_per_bin=5, n_timepoints=61, 
                         true_D=0.05, true_gamma=50.0, noise_level=2.0):
    """
    Create synthetic transcription and mRNA data with known constant degradation.
    
    Args:
        n_bins: Number of spatial bins
        n_traces_per_bin: Traces per bin
        n_timepoints: Number of time points
        true_D: True constant degradation rate
        true_gamma: True gamma scaling factor
        noise_level: Observation noise std
    
    Returns:
        F_data: transcription data array
        m_data: mRNA data with noise
        true_params: dict with ground truth
    """
    np.random.seed(42)
    
    # Time array in minutes
    t_array = np.arange(n_timepoints) * 20 / 60.0
    dt = t_array[1] - t_array[0]
    
    # Simple Gaussian transcription pulse
    t_peak = 10.0  # minutes
    F_data = np.exp(-((t_array - t_peak) ** 2) / (2 * 2.0**2))
    
    # Broadcast to all spatial positions
    total_traces = n_bins * n_traces_per_bin
    F_data = np.tile(F_data, (total_traces, 1))
    
    # Calculate expected mRNA using convolution with constant D
    # m(T) = gamma * integral(F(t) * exp(-D*(T-t)) dt)
    t_final = t_array[-1]
    survival_profile = np.exp(-true_D * (t_final - t_array))
    
    # Trapezoidal integration weights
    weights = np.ones_like(t_array)
    weights[0] = 0.5
    weights[-1] = 0.5
    weights_scaled = weights * dt
    
    # Expected mRNA for each trace (all same since F is same)
    integrand = F_data[0] * survival_profile
    m_expected = true_gamma * np.sum(integrand * weights_scaled)
    
    # Add noise to mRNA observations
    m_data = m_expected + np.random.normal(0, noise_level, total_traces)
    m_data = np.maximum(m_data, 0)  # Ensure non-negative
    
    return F_data, m_data, {'D': true_D, 'gamma': true_gamma, 'sigma': noise_level}


class TestNullConstantInference:
    """Test suite for null constant degradation inference script."""
    
    def test_synthetic_data_generation(self):
        """Test that synthetic data generation works."""
        F_data, m_data, true_params = create_synthetic_data()
        
        assert F_data.shape == (25, 61)
        assert m_data.shape == (25,)
        assert all(m_data >= 0)
        assert 'D' in true_params
        assert 'gamma' in true_params
        assert 'sigma' in true_params
    
    def test_script_cli_interface(self):
        """Test that script accepts required CLI arguments."""
        import subprocess
        
        script_path = Path(__file__).parent.parent / "scripts" / "02_infer_degradation_rates_null_constant.py"
        
        # Test help message
        result = subprocess.run(
            ["python", str(script_path), "--help"],
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0
        assert "--transcription" in result.stdout
        assert "--mrna" in result.stdout
        assert "--output-chain" in result.stdout
        assert "--output-trace" in result.stdout
    
    def test_output_csv_format(self):
        """Test that output CSV has correct format: D[0], gamma, sigma."""
        # Create temporary directory for outputs
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Create synthetic data files
            F_data, m_data, true_params = create_synthetic_data()
            
            trans_file = tmpdir / "transcription.csv"
            mrna_file = tmpdir / "mrna.csv"
            
            pd.DataFrame(F_data).to_csv(trans_file, index=False, header=False)
            pd.DataFrame(m_data).to_csv(mrna_file, index=False, header=False)
            
            # Output paths
            output_chain = tmpdir / "chains" / "degradation_chain.csv"
            output_trace = tmpdir / "trace.pdf"
            
            # Run inference with minimal samples for speed
            import subprocess
            script_path = Path(__file__).parent.parent / "scripts" / "02_infer_degradation_rates_null_constant.py"
            
            result = subprocess.run([
                "python", str(script_path),
                "--transcription", str(trans_file),
                "--mrna", str(mrna_file),
                "--output-chain", str(output_chain),
                "--output-trace", str(output_trace),
                "--n-samples", "100",  # Minimal for testing
                "--n-chains", "2"
            ], capture_output=True, text=True)
            
            # Check execution
            if result.returncode != 0:
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                pytest.fail(f"Script failed with return code {result.returncode}")
            
            # Check output files exist
            assert output_chain.exists(), f"Chain CSV not created at {output_chain}"
            assert output_trace.exists(), f"Trace plot not created at {output_trace}"
            
            # Check CSV format
            df = pd.read_csv(output_chain)
            
            # Must have columns: D[0], gamma, sigma
            assert "D[0]" in df.columns, f"Missing D[0] column. Found: {df.columns.tolist()}"
            assert "gamma" in df.columns, f"Missing gamma column. Found: {df.columns.tolist()}"
            assert "sigma" in df.columns, f"Missing sigma column. Found: {df.columns.tolist()}"
            
            # Check that we have samples (chains × draws rows)
            assert len(df) > 0, "No samples in output"
            assert len(df) == 200, f"Expected 200 samples (2 chains × 100 draws), got {len(df)}"
            
            # Check that D[0] values are reasonable (positive)
            assert (df["D[0]"] > 0).all(), "Some D[0] values are non-positive"
    
    def test_parameter_recovery(self):
        """Test that model can recover parameters from synthetic data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Create synthetic data with known parameters
            true_D = 0.05  # min^-1
            true_gamma = 50.0
            F_data, m_data, true_params = create_synthetic_data(
                true_D=true_D, 
                true_gamma=true_gamma,
                noise_level=1.0  # Low noise for better recovery
            )
            
            trans_file = tmpdir / "transcription.csv"
            mrna_file = tmpdir / "mrna.csv"
            
            pd.DataFrame(F_data).to_csv(trans_file, index=False, header=False)
            pd.DataFrame(m_data).to_csv(mrna_file, index=False, header=False)
            
            output_chain = tmpdir / "chains" / "degradation_chain.csv"
            output_trace = tmpdir / "trace.pdf"
            
            # Run inference
            import subprocess
            script_path = Path(__file__).parent.parent / "scripts" / "02_infer_degradation_rates_null_constant.py"
            
            result = subprocess.run([
                "python", str(script_path),
                "--transcription", str(trans_file),
                "--mrna", str(mrna_file),
                "--output-chain", str(output_chain),
                "--output-trace", str(output_trace),
                "--n-samples", "500",  # More samples for better estimates
                "--n-chains", "2",
                "--target-accept", "0.9"
            ], capture_output=True, text=True, timeout=120)
            
            if result.returncode != 0:
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                pytest.fail(f"Script failed with return code {result.returncode}")
            
            # Load results
            df = pd.read_csv(output_chain)
            
            # Check posterior means are in reasonable range
            D_mean = df["D[0]"].mean()
            gamma_mean = df["gamma"].mean()
            
            # Very loose bounds - just check ballpark
            assert 0.01 < D_mean < 0.5, f"D mean {D_mean} outside reasonable range"
            assert 10 < gamma_mean < 200, f"gamma mean {gamma_mean} outside reasonable range"
            
            print(f"\nTrue D: {true_D:.3f}, Estimated D: {D_mean:.3f}")
            print(f"True gamma: {true_gamma:.1f}, Estimated gamma: {gamma_mean:.1f}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
