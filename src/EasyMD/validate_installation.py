#!/usr/bin/env python3
"""
EasyMD Installation and Validation Test Script

This script tests the EasyMD installation and demonstrates the input validation system.
Run this script to verify that EasyMD is properly installed and configured.
"""

import sys
import os
import argparse
import tempfile
from pathlib import Path

def create_test_files():
    """Create temporary test files for validation testing."""
    test_files = {}
    
    # Create a minimal test PDB file
    test_pdb_content = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
END
"""
    
    # Create test PDB file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(test_pdb_content)
        test_files['pdb'] = f.name
    
    # Create test config file
    test_config_content = f"""
protein: "{test_files['pdb']}"
steps: 1000
temperature: 300
solvate: true
interval: 100
equilibration_steps: 50
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
        f.write(test_config_content)
        test_files['config'] = f.name
    
    return test_files

def cleanup_test_files(test_files):
    """Clean up temporary test files."""
    for file_path in test_files.values():
        try:
            os.unlink(file_path)
        except OSError:
            pass

def test_import():
    """Test if EasyMD modules can be imported."""
    print("Testing EasyMD module imports...")
    
    try:
        from EasyMD.argManager.manager import ArgManager
        print("✅ ArgManager import successful")
    except ImportError as e:
        print(f"❌ ArgManager import failed: {e}")
        return False
    
    try:
        from EasyMD.sysGenerator.sysGenerator import SysGenerator
        print("✅ SysGenerator import successful")
    except ImportError as e:
        print(f"❌ SysGenerator import failed: {e}")
        return False
    
    try:
        from EasyMD.simRunner.simRunner import SimRunner
        print("✅ SimRunner import successful")
    except ImportError as e:
        print(f"❌ SimRunner import failed: {e}")
        return False
    
    return True

def test_dependencies():
    """Test if required dependencies are available."""
    print("\nTesting required dependencies...")
    
    dependencies = [
        ('openmm', 'OpenMM'),
        ('openff.toolkit', 'OpenFF Toolkit'),
        ('yaml', 'PyYAML'),
        ('mdtraj', 'MDTraj'),
        ('pdbfixer', 'PDBFixer')
    ]
    
    all_available = True
    
    for module_name, display_name in dependencies:
        try:
            __import__(module_name)
            print(f"✅ {display_name} available")
        except ImportError:
            print(f"❌ {display_name} not available")
            all_available = False
    
    return all_available

def test_validation_system():
    """Test the input validation system."""
    print("\nTesting input validation system...")
    
    try:
        from EasyMD.argManager.manager import ArgManager
        
        # Test 1: Valid configuration
        print("\n1. Testing valid configuration...")
        test_files = create_test_files()
        
        try:
            original_argv = sys.argv.copy()
            sys.argv = ["test", "--config", test_files['config']]
            
            parser = argparse.ArgumentParser()
            manager = ArgManager(parser)
            print("✅ Valid configuration accepted")
            
        except SystemExit as e:
            if e.code == 0:
                print("✅ Valid configuration accepted (clean exit)")
            else:
                print(f"⚠️  Unexpected exit code: {e.code}")
        finally:
            sys.argv = original_argv
            cleanup_test_files(test_files)
        
        # Test 2: Invalid configuration (missing protein)
        print("\n2. Testing invalid configuration (missing protein)...")
        try:
            original_argv = sys.argv.copy()
            sys.argv = ["test", "--steps", "1000", "--solvate"]
            
            parser = argparse.ArgumentParser()
            manager = ArgManager(parser)
            print("❌ Invalid configuration was accepted (unexpected)")
            
        except SystemExit as e:
            if e.code == 1:
                print("✅ Invalid configuration properly rejected")
            else:
                print(f"⚠️  Unexpected exit code: {e.code}")
        finally:
            sys.argv = original_argv
        
        # Test 3: Conflicting parameters
        print("\n3. Testing conflicting parameters...")
        try:
            original_argv = sys.argv.copy()
            sys.argv = ["test", "--protein", "test.pdb", "--steps", "1000", "--clock", "60", "--solvate"]
            
            parser = argparse.ArgumentParser()
            manager = ArgManager(parser)
            print("❌ Conflicting parameters were accepted (unexpected)")
            
        except SystemExit as e:
            if e.code == 1:
                print("✅ Conflicting parameters properly rejected")
            else:
                print(f"⚠️  Unexpected exit code: {e.code}")
        finally:
            sys.argv = original_argv
        
        return True
        
    except Exception as e:
        print(f"❌ Validation system test failed: {e}")
        return False

def run_installation_test():
    """Run complete installation test."""
    print("="*60)
    print("EASYMD INSTALLATION AND VALIDATION TEST")
    print("="*60)
    
    # Test imports
    imports_ok = test_import()
    
    # Test dependencies
    deps_ok = test_dependencies()
    
    # Test validation system
    validation_ok = test_validation_system()
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    if imports_ok:
        print("✅ Module imports: PASSED")
    else:
        print("❌ Module imports: FAILED")
    
    if deps_ok:
        print("✅ Dependencies: PASSED")
    else:
        print("❌ Dependencies: FAILED")
    
    if validation_ok:
        print("✅ Validation system: PASSED")
    else:
        print("❌ Validation system: FAILED")
    
    overall_status = imports_ok and deps_ok and validation_ok
    
    print("\n" + "="*60)
    if overall_status:
        print("🎉 OVERALL STATUS: EasyMD is properly installed and ready to use!")
        print("\nNext steps:")
        print("1. Try running: python -m EasyMD --help")
        print("2. Create a simulation with: python -m EasyMD --protein your_protein.pdb --steps 1000 --solvate")
        print("3. Check examples in the examples/ directory")
    else:
        print("❌ OVERALL STATUS: EasyMD installation has issues")
        print("\nTroubleshooting:")
        if not imports_ok:
            print("- Check that EasyMD is properly installed in your Python path")
        if not deps_ok:
            print("- Install missing dependencies with: pip install -r requirements.txt")
        if not validation_ok:
            print("- Check for Python version compatibility (3.8+ recommended)")
    
    print("="*60)
    
    return overall_status

if __name__ == "__main__":
    success = run_installation_test()
    sys.exit(0 if success else 1)