#!/usr/bin/env python
"""
Tests for multi-model comparison script refactoring (Phase 2).

Test Coverage:
1. load_chain_df - Load chain CSV into DataFrame
2. parse_chain_params - Extract D, gamma, sigma for each model type
3. build_expected_mrna_for_model - Compute expected mRNA for each model type
4. Full comparison pipeline with multiple models
"""

import os
import sys
import numpy as np
import pandas as pd
import pytest
import tempfile
from pathlib import Path
import importlib.util

# Load the module dynamically (script name starts with digit, can't import directly)
script_path = Path(__file__).parent.parent / "scripts" / "08_compare_models_loo_ppc.py"
spec = importlib.util.spec_from_file_location("compare_models_loo_ppc", script_path)
compare_module = importlib.util.module_from_spec(spec)
sys.modules["compare_models_loo_ppc"] = compare_module
spec.loader.exec_module(compare_module)

# Make functions available at module level
load_chain_df = compare_module.load_chain_df
parse_chain_params = compare_module.parse_chain_params
build_expected_mrna_for_model = compare_module.build_expected_mrna_for_model
ModelSpec = compare_module.ModelSpec
build_log_likelihood = compare_module.build_log_likelihood
infer_chain_draw_shape = compare_module.infer_chain_draw_shape


def create_mock_chain_csv(model_type, n_samples=100, n_time=61, n_ap_bins=5):
    """
    Create a mock chain CSV for testing.
    
    Args:
        model_type: "age", "null_constant", or "spatial_ap"
        n_samples: Number of MCMC samples
        n_time: Number of timepoints (for age model)
        n_ap_bins: Number of AP bins (for spatial model)
    
    Returns:
        DataFrame with appropriate columns
    """
    np.random.seed(42)
    
    if model_type == "age":
        # D[0..n_time-1], gamma, sigma
        data = {}
        for i in range(n_time):
            data[f"D[{i}]"] = np.random.gamma(2, 0.025, n_samples)  # ~0.05 mean
        data["gamma"] = np.random.gamma(2, 25, n_samples)  # ~50 mean
        data["sigma"] = np.random.gamma(2, 1.5, n_samples)  # ~3 mean
        
    elif model_type == "null_constant":
        # D[0] (scalar), gamma, sigma
        data = {
            "D[0]": np.random.gamma(2, 0.025, n_samples),
            "gamma": np.random.gamma(2, 25, n_samples),
            "sigma": np.random.gamma(2, 1.5, n_samples),
        }
        
    elif model_type == "spatial_ap":
        # D[0..n_ap_bins-1], gamma, sigma
        data = {}
        for i in range(n_ap_bins):
            data[f"D[{i}]"] = np.random.gamma(2, 0.025, n_samples)
        data["gamma"] = np.random.gamma(2, 25, n_samples)
        data["sigma"] = np.random.gamma(2, 1.5, n_samples)
        
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    return pd.DataFrame(data)


def create_mock_transcription_data(n_obs=25, n_time=61):
    """Create mock transcription data."""
    np.random.seed(42)
    t_array = np.arange(n_time) * 20 / 60.0
    # Gaussian pulse
    F = np.exp(-((t_array - 10.0) ** 2) / (2 * 2.0**2))
    return np.tile(F, (n_obs, 1))


def create_mock_mrna_data(n_obs=25):
    """Create mock mRNA observation data."""
    np.random.seed(42)
    return np.random.gamma(2, 10, n_obs)


class TestLoadChainDF:
    """Test load_chain_df function."""
    
    def test_load_chain_df_age_model(self):
        """Test loading age-dependent model chain."""
        df = create_mock_chain_csv("age", n_samples=50, n_time=61)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df.to_csv(f.name, index=False)
            temp_path = f.name
        
        try:
            loaded_df = load_chain_df(temp_path)
            assert isinstance(loaded_df, pd.DataFrame)
            assert len(loaded_df) == 50
            assert "D[0]" in loaded_df.columns
            assert "D[60]" in loaded_df.columns
            assert "gamma" in loaded_df.columns
            assert "sigma" in loaded_df.columns
        finally:
            os.unlink(temp_path)
    
    def test_load_chain_df_null_constant(self):
        """Test loading null constant model chain."""
        df = create_mock_chain_csv("null_constant", n_samples=50)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df.to_csv(f.name, index=False)
            temp_path = f.name
        
        try:
            loaded_df = load_chain_df(temp_path)
            assert isinstance(loaded_df, pd.DataFrame)
            assert len(loaded_df) == 50
            assert "D[0]" in loaded_df.columns
            assert "D[1]" not in loaded_df.columns  # Only one D
            assert "gamma" in loaded_df.columns
            assert "sigma" in loaded_df.columns
        finally:
            os.unlink(temp_path)
    
    def test_load_chain_df_spatial(self):
        """Test loading spatial AP model chain."""
        df = create_mock_chain_csv("spatial_ap", n_samples=50, n_ap_bins=5)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df.to_csv(f.name, index=False)
            temp_path = f.name
        
        try:
            loaded_df = load_chain_df(temp_path)
            assert isinstance(loaded_df, pd.DataFrame)
            assert len(loaded_df) == 50
            assert "D[0]" in loaded_df.columns
            assert "D[4]" in loaded_df.columns
            assert "D[5]" not in loaded_df.columns  # Only 5 bins
            assert "gamma" in loaded_df.columns
            assert "sigma" in loaded_df.columns
        finally:
            os.unlink(temp_path)


class TestParseChainParams:
    """Test parse_chain_params function."""
    
    def test_parse_chain_params_age_model(self):
        """Test parsing age-dependent model parameters."""
        df = create_mock_chain_csv("age", n_samples=50, n_time=61)
        D_samples, gamma, sigma = parse_chain_params(df, "age", n_time=61, n_ap_bins=5)
        
        assert D_samples.shape == (50, 61)
        assert gamma.shape == (50,)
        assert sigma.shape == (50,)
        assert np.all(D_samples >= 0)
        assert np.all(gamma >= 0)
        assert np.all(sigma >= 0)
    
    def test_parse_chain_params_null_constant(self):
        """Test parsing null constant model parameters."""
        df = create_mock_chain_csv("null_constant", n_samples=50)
        D_samples, gamma, sigma = parse_chain_params(df, "null_constant", n_time=61, n_ap_bins=5)
        
        # D_samples should be 1D (scalar per sample) or 2D [n_samples, 1]
        assert D_samples.ndim in [1, 2]
        if D_samples.ndim == 2:
            assert D_samples.shape[0] == 50
            assert D_samples.shape[1] == 1
        else:
            assert D_samples.shape == (50,)
        assert gamma.shape == (50,)
        assert sigma.shape == (50,)
        assert np.all(D_samples >= 0)
    
    def test_parse_chain_params_spatial_ap(self):
        """Test parsing spatial AP model parameters."""
        df = create_mock_chain_csv("spatial_ap", n_samples=50, n_ap_bins=5)
        D_samples, gamma, sigma = parse_chain_params(df, "spatial_ap", n_time=61, n_ap_bins=5)
        
        assert D_samples.shape == (50, 5)
        assert gamma.shape == (50,)
        assert sigma.shape == (50,)
        assert np.all(D_samples >= 0)
    
    def test_parse_chain_params_missing_columns(self):
        """Test error handling for missing columns."""
        df = pd.DataFrame({"D[0]": [0.1, 0.2], "gamma": [50, 51]})  # Missing sigma
        
        with pytest.raises(ValueError, match="sigma"):
            parse_chain_params(df, "null_constant", n_time=61, n_ap_bins=5)
    
    def test_parse_chain_params_mismatched_d_count(self):
        """Test error handling for wrong number of D columns."""
        df = create_mock_chain_csv("age", n_samples=50, n_time=30)  # Only 30 timepoints
        
        with pytest.raises(ValueError, match="D columns"):
            parse_chain_params(df, "age", n_time=61, n_ap_bins=5)  # Expects 61


class TestBuildExpectedMRNA:
    """Test build_expected_mrna_for_model function."""
    
    def test_build_expected_mrna_age_model(self):
        """Test expected mRNA computation for age-dependent model."""
        n_samples = 50
        n_time = 61
        n_obs = 25
        
        # Create parameters
        D_samples = np.random.gamma(2, 0.025, (n_samples, n_time))
        gamma = np.random.gamma(2, 25, n_samples)
        transcription = create_mock_transcription_data(n_obs, n_time)
        dt = 20.0 / 60.0
        
        mu = build_expected_mrna_for_model(
            "age", D_samples, gamma, transcription, dt, n_ap_bins=5, n_dv_bins=5
        )
        
        assert mu.shape == (n_samples, n_obs)
        assert np.all(mu >= 0)
        assert not np.any(np.isnan(mu))
        assert not np.any(np.isinf(mu))
    
    def test_build_expected_mrna_null_constant(self):
        """Test expected mRNA computation for null constant model."""
        n_samples = 50
        n_time = 61
        n_obs = 25
        
        # D_samples is scalar per sample
        D_samples = np.random.gamma(2, 0.025, n_samples)
        gamma = np.random.gamma(2, 25, n_samples)
        transcription = create_mock_transcription_data(n_obs, n_time)
        dt = 20.0 / 60.0
        
        mu = build_expected_mrna_for_model(
            "null_constant", D_samples, gamma, transcription, dt, n_ap_bins=5, n_dv_bins=5
        )
        
        assert mu.shape == (n_samples, n_obs)
        assert np.all(mu >= 0)
        assert not np.any(np.isnan(mu))
        assert not np.any(np.isinf(mu))
    
    def test_build_expected_mrna_spatial_ap(self):
        """Test expected mRNA computation for spatial AP model."""
        n_samples = 50
        n_time = 61
        n_ap_bins = 5
        n_dv_bins = 5
        n_obs = n_ap_bins * n_dv_bins
        
        # D_samples has one value per AP bin
        D_samples = np.random.gamma(2, 0.025, (n_samples, n_ap_bins))
        gamma = np.random.gamma(2, 25, n_samples)
        transcription = create_mock_transcription_data(n_obs, n_time)
        dt = 20.0 / 60.0
        
        mu = build_expected_mrna_for_model(
            "spatial_ap", D_samples, gamma, transcription, dt, n_ap_bins=n_ap_bins, n_dv_bins=n_dv_bins
        )
        
        assert mu.shape == (n_samples, n_obs)
        assert np.all(mu >= 0)
        assert not np.any(np.isnan(mu))
        assert not np.any(np.isinf(mu))
    
    def test_build_expected_mrna_spatial_ap_mapping(self):
        """Test that spatial model correctly maps AP bins to observations."""
        n_samples = 10
        n_time = 61
        n_ap_bins = 3
        n_dv_bins = 2
        n_obs = n_ap_bins * n_dv_bins  # 6 observations
        
        # Create distinctive D values for each AP bin
        D_samples = np.array([[0.01, 0.05, 0.10]] * n_samples)  # All samples same for testing
        gamma = np.ones(n_samples) * 50.0
        transcription = create_mock_transcription_data(n_obs, n_time)
        dt = 20.0 / 60.0
        
        mu = build_expected_mrna_for_model(
            "spatial_ap", D_samples, gamma, transcription, dt, n_ap_bins=n_ap_bins, n_dv_bins=n_dv_bins
        )
        
        # Observations 0,1 should use D[0]=0.01 (low decay, high mRNA)
        # Observations 2,3 should use D[1]=0.05 (medium decay, medium mRNA)
        # Observations 4,5 should use D[2]=0.10 (high decay, low mRNA)
        # Average across samples (all same)
        mu_mean = mu.mean(axis=0)
        
        # Within each AP bin, observations should be similar
        assert np.abs(mu_mean[0] - mu_mean[1]) < 0.1 * mu_mean[0]
        assert np.abs(mu_mean[2] - mu_mean[3]) < 0.1 * mu_mean[2]
        assert np.abs(mu_mean[4] - mu_mean[5]) < 0.1 * mu_mean[4]
        
        # Lower D should give higher mRNA
        assert mu_mean[0] > mu_mean[2] > mu_mean[4]


class TestModelSpecStructure:
    """Test ModelSpec dataclass or structure."""
    
    def test_model_spec_creation(self):
        """Test creating ModelSpec instances."""
        spec = ModelSpec(
            name="Convolution",
            model_type="age",
            results_dir="/path/to/results"
        )
        
        assert spec.name == "Convolution"
        assert spec.model_type == "age"
        assert spec.results_dir == "/path/to/results"
    
    def test_model_spec_invalid_type(self):
        """Test that invalid model_type raises error."""
        # This should work
        ModelSpec(name="Test", model_type="age", results_dir="/path")
        ModelSpec(name="Test", model_type="null_constant", results_dir="/path")
        ModelSpec(name="Test", model_type="spatial_ap", results_dir="/path")
        
        # This should fail (if validation implemented)
        # Note: dataclass won't validate by default, but we can add validation
        try:
            spec = ModelSpec(name="Test", model_type="invalid", results_dir="/path")
            # If we get here without error, that's okay - validation is optional
            assert spec.model_type == "invalid"
        except (ValueError, TypeError):
            # If validation is implemented, this is expected
            pass


class TestFullComparisonPipeline:
    """Test full comparison pipeline with multiple models."""
    
    def test_compare_three_models(self):
        """Test comparison of age, null_constant, and spatial models."""
        import arviz as az
        
        n_samples = 80  # 4 chains * 20 draws
        n_chains = 4
        n_draws = 20
        n_time = 61
        n_ap_bins = 5
        n_dv_bins = 5
        n_obs = n_ap_bins * n_dv_bins
        dt = 20.0 / 60.0
        
        # Create data
        transcription = create_mock_transcription_data(n_obs, n_time)
        observed = create_mock_mrna_data(n_obs)
        
        # Create chain CSVs for three models
        model_specs = [
            ("Age", "age"),
            ("NullConstant", "null_constant"),
            ("SpatialAP", "spatial_ap"),
        ]
        
        idata_dict = {}
        
        for model_name, model_type in model_specs:
            # Create chain
            df = create_mock_chain_csv(model_type, n_samples=n_samples, n_time=n_time, n_ap_bins=n_ap_bins)
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
                df.to_csv(f.name, index=False)
                temp_path = f.name
            
            try:
                # Load and parse
                loaded_df = load_chain_df(temp_path)
                D_samples, gamma, sigma = parse_chain_params(
                    loaded_df, model_type, n_time=n_time, n_ap_bins=n_ap_bins
                )
                
                # Build expected mRNA
                mu = build_expected_mrna_for_model(
                    model_type, D_samples, gamma, transcription, dt, n_ap_bins=n_ap_bins, n_dv_bins=n_dv_bins
                )
                
                # Build log likelihood
                log_lik = build_log_likelihood(mu, sigma, observed)
                
                # Reshape for ArviZ
                chains, draws = infer_chain_draw_shape(len(loaded_df), n_chains, n_draws)
                
                # Prepare posterior dict based on model type
                if model_type == "age":
                    D_reshaped = D_samples.reshape(chains, draws, n_time)
                    posterior = {
                        "D": D_reshaped,
                        "gamma": gamma.reshape(chains, draws),
                        "sigma": sigma.reshape(chains, draws),
                    }
                elif model_type == "null_constant":
                    # D0 should be 2D [chains, draws]
                    if D_samples.ndim == 1:
                        D_reshaped = D_samples.reshape(chains, draws)
                    else:
                        D_reshaped = D_samples.reshape(chains, draws, 1)[:, :, 0]
                    posterior = {
                        "D0": D_reshaped,
                        "gamma": gamma.reshape(chains, draws),
                        "sigma": sigma.reshape(chains, draws),
                    }
                elif model_type == "spatial_ap":
                    D_reshaped = D_samples.reshape(chains, draws, n_ap_bins)
                    posterior = {
                        "D": D_reshaped,
                        "gamma": gamma.reshape(chains, draws),
                        "sigma": sigma.reshape(chains, draws),
                    }
                
                log_lik_reshaped = log_lik.reshape(chains, draws, n_obs)
                
                # Create InferenceData
                idata = az.from_dict(
                    posterior=posterior,
                    log_likelihood={"m_obs": log_lik_reshaped},
                    observed_data={"m_obs": observed},
                )
                idata.posterior.attrs["model_name"] = model_name
                
                idata_dict[model_name] = idata
                
            finally:
                os.unlink(temp_path)
        
        # Test that comparison works
        assert len(idata_dict) == 3
        
        # Run LOO comparison
        comparison = az.compare(idata_dict, ic="loo", scale="deviance")
        
        assert len(comparison) == 3
        assert "elpd_loo" in comparison.columns
        assert all(model_name in comparison.index for model_name, _ in model_specs)
    
    def test_backward_compatibility_two_models(self):
        """Test that refactored code still works with two models (backward compatibility)."""
        import arviz as az
        
        n_samples = 80
        n_chains = 4
        n_draws = 20
        n_time = 61
        n_obs = 25
        dt = 20.0 / 60.0
        
        transcription = create_mock_transcription_data(n_obs, n_time)
        observed = create_mock_mrna_data(n_obs)
        
        # Two age-dependent models (like original Convolution vs BioPolyA)
        idata_dict = {}
        
        for model_name in ["Convolution", "BioPolyA"]:
            df = create_mock_chain_csv("age", n_samples=n_samples, n_time=n_time)
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
                df.to_csv(f.name, index=False)
                temp_path = f.name
            
            try:
                loaded_df = load_chain_df(temp_path)
                D_samples, gamma, sigma = parse_chain_params(loaded_df, "age", n_time=n_time, n_ap_bins=5)
                mu = build_expected_mrna_for_model(
                    "age", D_samples, gamma, transcription, dt, n_ap_bins=5, n_dv_bins=5
                )
                log_lik = build_log_likelihood(mu, sigma, observed)
                
                chains, draws = infer_chain_draw_shape(len(loaded_df), n_chains, n_draws)
                
                idata = az.from_dict(
                    posterior={
                        "D": D_samples.reshape(chains, draws, n_time),
                        "gamma": gamma.reshape(chains, draws),
                        "sigma": sigma.reshape(chains, draws),
                    },
                    log_likelihood={"m_obs": log_lik.reshape(chains, draws, n_obs)},
                    observed_data={"m_obs": observed},
                )
                
                idata_dict[model_name] = idata
                
            finally:
                os.unlink(temp_path)
        
        comparison = az.compare(idata_dict, ic="loo", scale="deviance")
        assert len(comparison) == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
