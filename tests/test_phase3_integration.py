#!/usr/bin/env python
"""
Integration test for Phase 3 model comparison - backward compatibility and new features.
"""
import os
import sys
import tempfile
import json
import shutil
import numpy as np
import pandas as pd
from pathlib import Path

# Get the script directory
script_path = Path(__file__).parent.parent / "scripts" / "08_compare_models_loo_ppc.py"

def create_mock_chain_csv(output_path, model_type="age", n_samples=100, n_time=61, n_ap_bins=5):
    """Create a mock chain CSV file for testing."""
    np.random.seed(42)
    
    if model_type == "age":
        data = {}
        for i in range(n_time):
            data[f"D[{i}]"] = np.random.gamma(2, 0.025, n_samples)
        data["gamma"] = np.random.gamma(2, 25, n_samples)
        data["sigma"] = np.random.gamma(2, 1.5, n_samples)
    elif model_type == "null_constant":
        data = {
            "D[0]": np.random.gamma(2, 0.025, n_samples),
            "gamma": np.random.gamma(2, 25, n_samples),
            "sigma": np.random.gamma(2, 1.5, n_samples),
        }
    elif model_type == "spatial_ap":
        data = {}
        for i in range(n_ap_bins):
            data[f"D[{i}]"] = np.random.gamma(2, 0.025, n_samples)
        data["gamma"] = np.random.gamma(2, 25, n_samples)
        data["sigma"] = np.random.gamma(2, 1.5, n_samples)
    
    df = pd.DataFrame(data)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)

def create_mock_data_files(temp_dir):
    """Create mock transcription and mRNA data files."""
    n_obs = 25
    n_time = 61
    
    # Transcription data
    np.random.seed(42)
    t_array = np.arange(n_time) * 20 / 60.0
    F = np.exp(-((t_array - 10.0) ** 2) / (2 * 2.0**2))
    transcription = np.tile(F, (n_obs, 1))
    trans_path = os.path.join(temp_dir, "transcription.csv")
    pd.DataFrame(transcription).to_csv(trans_path, header=False, index=False)
    
    # mRNA data
    mrna = np.random.gamma(2, 10, n_obs)
    mrna_path = os.path.join(temp_dir, "mrna.csv")
    pd.DataFrame(mrna).to_csv(mrna_path, header=False, index=False)
    
    return trans_path, mrna_path

def test_backward_compatibility():
    """Test backward compatible 2-model comparison via CLI."""
    print("Testing backward compatibility (2 models via CLI)...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create mock model directories
        conv_dir = os.path.join(temp_dir, "conv_model")
        bio_dir = os.path.join(temp_dir, "bio_model")
        output_dir = os.path.join(temp_dir, "output")
        
        # Create mock chain files
        create_mock_chain_csv(os.path.join(conv_dir, "chains", "degradation_chain.csv"), "age")
        create_mock_chain_csv(os.path.join(bio_dir, "chains", "degradation_chain.csv"), "age")
        
        # Create mock data
        trans_path, mrna_path = create_mock_data_files(temp_dir)
        
        # Run script
        cmd = f"""
        python {script_path} \\
            --convolution-dir {conv_dir} \\
            --biopolya-dir {bio_dir} \\
            --transcription {trans_path} \\
            --mrna {mrna_path} \\
            --output-dir {output_dir} \\
            --n-chains 4 \\
            --n-draws 25
        """
        
        result = os.system(cmd)
        
        if result == 0 and os.path.exists(os.path.join(output_dir, "loo_comparison.csv")):
            print("✓ Backward compatibility test PASSED")
            return True
        else:
            print("✗ Backward compatibility test FAILED")
            return False

def test_four_model_cli():
    """Test 4-model comparison via CLI args."""
    print("\nTesting 4-model comparison via CLI...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create mock model directories
        conv_dir = os.path.join(temp_dir, "conv_model")
        bio_dir = os.path.join(temp_dir, "bio_model")
        null_dir = os.path.join(temp_dir, "null_model")
        spatial_dir = os.path.join(temp_dir, "spatial_model")
        output_dir = os.path.join(temp_dir, "output")
        
        # Create mock chain files
        create_mock_chain_csv(os.path.join(conv_dir, "chains", "degradation_chain.csv"), "age")
        create_mock_chain_csv(os.path.join(bio_dir, "chains", "degradation_chain.csv"), "age")
        create_mock_chain_csv(os.path.join(null_dir, "chains", "degradation_chain.csv"), "null_constant")
        create_mock_chain_csv(os.path.join(spatial_dir, "chains", "degradation_chain.csv"), "spatial_ap")
        
        # Create mock data
        trans_path, mrna_path = create_mock_data_files(temp_dir)
        
        # Run script
        cmd = f"""
        python {script_path} \\
            --convolution-dir {conv_dir} \\
            --biopolya-dir {bio_dir} \\
            --null-dir {null_dir} \\
            --spatial-dir {spatial_dir} \\
            --transcription {trans_path} \\
            --mrna {mrna_path} \\
            --output-dir {output_dir} \\
            --n-chains 4 \\
            --n-draws 25
        """
        
        result = os.system(cmd)
        
        # Check outputs
        comparison_exists = os.path.exists(os.path.join(output_dir, "loo_comparison.csv"))
        
        if result == 0 and comparison_exists:
            # Check that all 4 models are in the comparison
            df = pd.read_csv(os.path.join(output_dir, "loo_comparison.csv"), index_col=0)
            if len(df) == 4:
                print("✓ 4-model CLI test PASSED")
                return True
        
        print("✗ 4-model CLI test FAILED")
        return False

def test_json_config():
    """Test JSON configuration file."""
    print("\nTesting JSON config file...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create mock model directories
        conv_dir = os.path.join(temp_dir, "conv_model")
        bio_dir = os.path.join(temp_dir, "bio_model")
        output_dir = os.path.join(temp_dir, "output")
        
        # Create mock chain files
        create_mock_chain_csv(os.path.join(conv_dir, "chains", "degradation_chain.csv"), "age")
        create_mock_chain_csv(os.path.join(bio_dir, "chains", "degradation_chain.csv"), "age")
        
        # Create mock data
        trans_path, mrna_path = create_mock_data_files(temp_dir)
        
        # Create JSON config
        config = {
            "chain_relative_path": "chains/degradation_chain.csv",
            "models": [
                {"name": "Convolution", "type": "age", "dir": conv_dir},
                {"name": "BioPolyA", "type": "age", "dir": bio_dir},
            ]
        }
        config_path = os.path.join(temp_dir, "models.json")
        with open(config_path, 'w') as f:
            json.dump(config, f)
        
        # Run script
        cmd = f"""
        python {script_path} \\
            --models-config {config_path} \\
            --transcription {trans_path} \\
            --mrna {mrna_path} \\
            --output-dir {output_dir} \\
            --n-chains 4 \\
            --n-draws 25
        """
        
        result = os.system(cmd)
        
        if result == 0 and os.path.exists(os.path.join(output_dir, "loo_comparison.csv")):
            print("✓ JSON config test PASSED")
            return True
        else:
            print("✗ JSON config test FAILED")
            return False

def main():
    print("="*60)
    print("Phase 3 Integration Tests")
    print("="*60)
    
    results = []
    
    # Run tests
    results.append(("Backward Compatibility", test_backward_compatibility()))
    results.append(("4-Model CLI", test_four_model_cli()))
    results.append(("JSON Config", test_json_config()))
    
    # Summary
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{name:30s} {status}")
    
    all_passed = all(r[1] for r in results)
    print("="*60)
    if all_passed:
        print("All integration tests PASSED ✓")
        return 0
    else:
        print("Some integration tests FAILED ✗")
        return 1

if __name__ == "__main__":
    sys.exit(main())
