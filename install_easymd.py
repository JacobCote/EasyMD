#!/usr/bin/env python3
"""
Alternative installation script for EasyMD.
Use this if the conda environment installation is having issues.
"""

import sys
import os
import subprocess
from pathlib import Path

def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"🔄 {description}...")
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(f"✅ {description} completed successfully")
        if result.stdout:
            print(f"   Output: {result.stdout.strip()}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed")
        print(f"   Error: {e.stderr.strip()}")
        return False

def main():
    """Main installation function."""
    print("EasyMD Alternative Installation Script")
    print("=" * 50)
    
    # Check if we're in the right directory
    if not Path("setup.py").exists():
        print("❌ Error: setup.py not found. Please run this script from the EasyMD root directory.")
        sys.exit(1)
    
    # Check if we're in a conda environment
    conda_env = os.environ.get('CONDA_DEFAULT_ENV')
    if not conda_env:
        print("⚠️  Warning: No conda environment detected. Make sure you've activated your environment.")
    else:
        print(f"📦 Installing in conda environment: {conda_env}")
    
    # Method 1: Try pip install in development mode
    print("\n🔧 Method 1: Installing with pip in development mode...")
    if run_command("pip install -e .", "Installing EasyMD in development mode"):
        print("\n🎉 Installation successful!")
        print("\nTesting installation...")
        if run_command("python -c 'import EasyMD; print(\"EasyMD imported successfully\")'", "Testing EasyMD import"):
            print("\n✅ EasyMD is ready to use!")
            print("\nYou can now run:")
            print("  easymd --help")
            print("  easymd info --help")
            print("  easymd analyze --help")
            return
    
    # Method 2: Try direct Python path setup
    print("\n🔧 Method 2: Setting up Python path...")
    current_dir = Path.cwd()
    src_dir = current_dir / "src"
    
    if src_dir.exists():
        print(f"Adding {src_dir} to Python path...")
        # Create a .pth file in site-packages
        try:
            import site
            site_packages = site.getsitepackages()[0]
            pth_file = Path(site_packages) / "easymd.pth"
            with open(pth_file, 'w') as f:
                f.write(str(src_dir) + '\n')
            print(f"✅ Created {pth_file}")
            
            # Test the import
            if run_command("python -c 'import EasyMD; print(\"EasyMD imported successfully\")'", "Testing EasyMD import"):
                print("\n✅ EasyMD is ready to use!")
                return
        except Exception as e:
            print(f"❌ Method 2 failed: {e}")
    
    # Method 3: Manual instructions
    print("\n🔧 Method 3: Manual setup instructions")
    print("If the automatic installation failed, you can set up EasyMD manually:")
    print(f"1. Add this to your Python path: {src_dir}")
    print("2. Or run Python scripts from the EasyMD directory using:")
    print("   python -m EasyMD --help")
    print("   python -m EasyMD info --help")
    print("   python -m EasyMD analyze --help")

if __name__ == "__main__":
    main()