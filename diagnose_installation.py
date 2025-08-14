#!/usr/bin/env python3
"""
Diagnostic script to identify why EasyMD works locally but not after git clone.
"""

import sys
import os
from pathlib import Path

def print_section(title):
    print(f"\n{'='*60}")
    print(f"🔍 {title}")
    print('='*60)

def check_python_path():
    print_section("PYTHON PATH ANALYSIS")
    print("Current Python path:")
    for i, path in enumerate(sys.path):
        print(f"  {i}: {path}")
    
    # Check if src directory is in path
    current_dir = Path.cwd()
    src_dir = current_dir / "src"
    
    print(f"\nCurrent directory: {current_dir}")
    print(f"Source directory: {src_dir}")
    print(f"Source directory exists: {src_dir.exists()}")
    print(f"Source directory in Python path: {str(src_dir) in sys.path}")

def check_file_structure():
    print_section("FILE STRUCTURE ANALYSIS")
    
    # Check if all necessary files exist
    files_to_check = [
        "setup.py",
        "requirements.yml",
        "src/EasyMD/__init__.py",
        "src/EasyMD/__main__.py",
        "src/EasyMD/sysGenerator/__init__.py",
        "src/EasyMD/sysGenerator/sysGenerator.py",
    ]
    
    for file_path in files_to_check:
        path = Path(file_path)
        exists = path.exists()
        print(f"  {'✅' if exists else '❌'} {file_path}")
        if exists and file_path.endswith('.py'):
            try:
                with open(path, 'r') as f:
                    content = f.read()
                    lines = len(content.splitlines())
                    print(f"      ({lines} lines)")
            except:
                print("      (could not read)")

def check_imports():
    print_section("IMPORT TESTING")
    
    # Test different import methods
    import_tests = [
        ("sys.path.insert(0, 'src')", lambda: sys.path.insert(0, 'src')),
        ("import EasyMD", lambda: __import__('EasyMD')),
        ("from EasyMD.sysGenerator.sysGenerator import SysGenerator", 
         lambda: __import__('EasyMD.sysGenerator.sysGenerator', fromlist=['SysGenerator'])),
        ("from EasyMD.sysGenerator import SysGenerator",
         lambda: __import__('EasyMD.sysGenerator', fromlist=['SysGenerator'])),
    ]
    
    for description, test_func in import_tests:
        try:
            test_func()
            print(f"  ✅ {description}")
        except Exception as e:
            print(f"  ❌ {description}")
            print(f"      Error: {e}")

def check_installation_state():
    print_section("INSTALLATION STATE")
    
    try:
        import pkg_resources
        try:
            dist = pkg_resources.get_distribution('EasyMD')
            print(f"  ✅ EasyMD is installed")
            print(f"      Version: {dist.version}")
            print(f"      Location: {dist.location}")
            print(f"      Egg info: {dist.egg_info}")
        except pkg_resources.DistributionNotFound:
            print(f"  ❌ EasyMD is not installed as a package")
    except ImportError:
        print(f"  ⚠️  pkg_resources not available")
    
    # Check if easymd command exists
    import shutil
    easymd_path = shutil.which('easymd')
    if easymd_path:
        print(f"  ✅ easymd command found at: {easymd_path}")
    else:
        print(f"  ❌ easymd command not found")

def check_environment():
    print_section("ENVIRONMENT INFORMATION")
    
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    
    # Check conda environment
    conda_env = os.environ.get('CONDA_DEFAULT_ENV')
    if conda_env:
        print(f"Conda environment: {conda_env}")
    else:
        print("Not in a conda environment")
    
    # Check working directory
    print(f"Working directory: {os.getcwd()}")

def main():
    print("EasyMD Installation Diagnostic Tool")
    print("This will help identify why EasyMD works locally but not after git clone")
    
    check_environment()
    check_python_path()
    check_file_structure()
    check_installation_state()
    check_imports()
    
    print_section("RECOMMENDATIONS")
    print("Based on the analysis above:")
    print("1. If 'Source directory in Python path' is False, that's likely the issue")
    print("2. If 'EasyMD is not installed as a package', run: pip install -e .")
    print("3. If imports fail, there might be circular import issues")
    print("4. Compare this output between your working local and fresh clone")

if __name__ == "__main__":
    main()