#!/usr/bin/env python3
"""
Fix script for EasyMD installation issues after git clone.
"""

import sys
import os
import subprocess
from pathlib import Path

def run_command(cmd, description, check=True):
    """Run a command and return success status."""
    print(f"🔄 {description}...")
    try:
        result = subprocess.run(cmd, shell=True, check=check, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ {description} - Success")
            if result.stdout.strip():
                print(f"   Output: {result.stdout.strip()}")
            return True
        else:
            print(f"❌ {description} - Failed")
            if result.stderr.strip():
                print(f"   Error: {result.stderr.strip()}")
            return False
    except Exception as e:
        print(f"❌ {description} - Exception: {e}")
        return False

def fix_method_1_clean_install():
    """Method 1: Clean pip installation."""
    print("\n" + "="*60)
    print("🔧 METHOD 1: Clean pip installation")
    print("="*60)
    
    # Uninstall any existing installation
    run_command("pip uninstall EasyMD -y", "Removing existing EasyMD installation", check=False)
    
    # Clean any cached files
    print("🧹 Cleaning cached files...")
    for root, dirs, files in os.walk("src"):
        for d in dirs[:]:  # Use slice to modify list during iteration
            if d == "__pycache__":
                import shutil
                shutil.rmtree(os.path.join(root, d))
                print(f"   Removed: {os.path.join(root, d)}")
                dirs.remove(d)
    
    # Install in development mode
    success = run_command("pip install -e .", "Installing EasyMD in development mode")
    
    if success:
        # Test the installation
        test_success = run_command(
            'python -c "import EasyMD; print(\'EasyMD imported successfully\')"',
            "Testing EasyMD import"
        )
        if test_success:
            run_command("easymd --help", "Testing easymd command", check=False)
            return True
    
    return False

def fix_method_2_python_path():
    """Method 2: Add src to Python path."""
    print("\n" + "="*60)
    print("🔧 METHOD 2: Python path setup")
    print("="*60)
    
    current_dir = Path.cwd()
    src_dir = current_dir / "src"
    
    if not src_dir.exists():
        print("❌ src directory not found")
        return False
    
    # Test with manual path addition
    test_code = f'''
import sys
sys.path.insert(0, "{src_dir}")
try:
    import EasyMD
    from EasyMD.sysGenerator.sysGenerator import SysGenerator
    print("✅ EasyMD imports work with manual path")
except Exception as e:
    print(f"❌ Import failed: {{e}}")
    sys.exit(1)
'''
    
    success = run_command(f'python -c "{test_code}"', "Testing with manual Python path")
    
    if success:
        # Create a startup script
        startup_script = '''#!/bin/bash
# EasyMD startup script
export PYTHONPATH="$PWD/src:$PYTHONPATH"
python -m EasyMD "$@"
'''
        with open("easymd_dev.sh", "w") as f:
            f.write(startup_script)
        os.chmod("easymd_dev.sh", 0o755)
        print("✅ Created easymd_dev.sh script")
        print("   You can use: ./easymd_dev.sh --help")
        return True
    
    return False

def fix_method_3_direct_execution():
    """Method 3: Direct module execution."""
    print("\n" + "="*60)
    print("🔧 METHOD 3: Direct module execution")
    print("="*60)
    
    # Test direct module execution
    success = run_command(
        'python -m EasyMD --help',
        "Testing direct module execution"
    )
    
    if success:
        print("✅ Direct module execution works!")
        print("   Use: python -m EasyMD --help")
        print("   Use: python -m EasyMD info --help")
        print("   Use: python -m EasyMD analyze --help")
        return True
    
    return False

def main():
    print("EasyMD Installation Fix Tool")
    print("=" * 50)
    
    # Check if we're in the right directory
    if not Path("setup.py").exists():
        print("❌ Error: setup.py not found. Please run this from the EasyMD root directory.")
        sys.exit(1)
    
    # Check conda environment
    conda_env = os.environ.get('CONDA_DEFAULT_ENV')
    if conda_env:
        print(f"📦 Conda environment: {conda_env}")
    else:
        print("⚠️  Warning: Not in a conda environment")
    
    # Try different fix methods
    methods = [
        fix_method_1_clean_install,
        fix_method_2_python_path,
        fix_method_3_direct_execution,
    ]
    
    for i, method in enumerate(methods, 1):
        try:
            if method():
                print(f"\n🎉 SUCCESS! Method {i} worked.")
                print("\nEasyMD should now be working. Try:")
                if i == 1:
                    print("  easymd --help")
                elif i == 2:
                    print("  ./easymd_dev.sh --help")
                else:
                    print("  python -m EasyMD --help")
                return
        except Exception as e:
            print(f"❌ Method {i} failed with exception: {e}")
    
    print("\n❌ All methods failed. Please check the diagnostic output above.")
    print("Consider running: python diagnose_installation.py")

if __name__ == "__main__":
    main()