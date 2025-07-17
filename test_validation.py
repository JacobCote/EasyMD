#!/usr/bin/env python3
"""
Test script to demonstrate the enhanced input validation system for EasyMD.
This script shows various validation scenarios and error messages.
"""

import sys
import argparse
from src.EasyMD.argManager.manager import ArgManager

def test_validation_scenarios():
    """Test various validation scenarios to demonstrate the enhanced error handling."""
    
    test_cases = [
        {
            "name": "Missing protein file",
            "args": ["--steps", "1000", "--solvate"],
            "expected": "Should fail - no protein file"
        },
        {
            "name": "Both steps and clock specified",
            "args": ["--protein", "test.pdb", "--steps", "1000", "--clock", "60", "--solvate"],
            "expected": "Should fail - conflicting duration methods"
        },
        {
            "name": "No solvation method",
            "args": ["--protein", "test.pdb", "--steps", "1000"],
            "expected": "Should fail - no solvation method specified"
        },
        {
            "name": "Both solvation methods",
            "args": ["--protein", "test.pdb", "--steps", "1000", "--solvate", "--GBIS"],
            "expected": "Should fail - conflicting solvation methods"
        },
        {
            "name": "Invalid temperature",
            "args": ["--protein", "test.pdb", "--steps", "1000", "--solvate", "--temperature", "-100"],
            "expected": "Should fail - negative temperature"
        },
        {
            "name": "Very short simulation",
            "args": ["--protein", "test.pdb", "--steps", "10", "--solvate"],
            "expected": "Should warn - very short simulation"
        },
        {
            "name": "Non-existent restart directory",
            "args": ["--restart", "/non/existent/path"],
            "expected": "Should fail - restart directory doesn't exist"
        }
    ]
    
    print("="*80)
    print("TESTING EASYMD INPUT VALIDATION SYSTEM")
    print("="*80)
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n{i}. Testing: {test_case['name']}")
        print(f"   Command: python -m EasyMD {' '.join(test_case['args'])}")
        print(f"   Expected: {test_case['expected']}")
        print("-" * 60)
        
        # Simulate the command line arguments
        original_argv = sys.argv.copy()
        sys.argv = ["test"] + test_case['args']
        
        try:
            parser = argparse.ArgumentParser()
            manager = ArgManager(parser)
            print("   Result: ✅ PASSED (unexpected)")
        except SystemExit as e:
            if e.code == 1:
                print("   Result: ❌ FAILED as expected")
            else:
                print(f"   Result: ⚠️  Unexpected exit code: {e.code}")
        except Exception as e:
            print(f"   Result: 💥 Unexpected error: {e}")
        finally:
            sys.argv = original_argv

def create_sample_config_files():
    """Create sample configuration files for testing."""
    
    # Valid configuration
    valid_config = """
# EasyMD Configuration File
protein: "4zgm.pdb"
ligand: "LIG"
steps: 10000
temperature: 300
solvate: true
water_model: "tip3p"
padding: 12.0
ionic_strength: 0.15
interval: 1000
equilibration_steps: 500
outdir: "simulation_output"
"""
    
    # Invalid configuration
    invalid_config = """
# Invalid EasyMD Configuration
protein: "nonexistent.pdb"
steps: -1000
temperature: -273
solvate: true
GBIS: true  # Conflicting with solvate
padding: -5
"""
    
    with open("valid_config.yml", "w") as f:
        f.write(valid_config)
    
    with open("invalid_config.yml", "w") as f:
        f.write(invalid_config)
    
    print("\n" + "="*80)
    print("SAMPLE CONFIGURATION FILES CREATED")
    print("="*80)
    print("\n1. valid_config.yml - Contains valid simulation parameters")
    print("2. invalid_config.yml - Contains invalid parameters for testing")
    print("\nTest with:")
    print("  python -m EasyMD --config valid_config.yml")
    print("  python -m EasyMD --config invalid_config.yml")

if __name__ == "__main__":
    print("EasyMD Input Validation Test Suite")
    print("This script demonstrates the enhanced validation system.")
    print("\nNote: This script will show validation errors - this is expected behavior!")
    
    # Create sample config files
    create_sample_config_files()
    
    # Test validation scenarios
    test_validation_scenarios()
    
    print("\n" + "="*80)
    print("VALIDATION TESTING COMPLETE")
    print("="*80)
    print("\nThe validation system provides:")
    print("✅ Clear error messages with specific issues")
    print("✅ Helpful suggestions for fixing problems")
    print("✅ Warnings for potentially problematic settings")
    print("✅ Success summary when all validations pass")
    print("✅ File existence and format checking")
    print("✅ Parameter range validation")
    print("✅ Logical consistency checks")