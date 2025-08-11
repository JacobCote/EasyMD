# EasyMD Installation Troubleshooting Guide

If you're experiencing import errors when installing EasyMD, try these solutions in order:

## Method 1: Standard Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/JacobCote/EasyMD.git
cd EasyMD

# Create and activate conda environment
conda env create -f requirements.yml
conda activate mdEnv

# Test the installation
easymd --help
```

## Method 2: Manual pip installation

If Method 1 fails, try installing manually:

```bash
# After activating the conda environment
pip install -e .

# Test the installation
python -c "import EasyMD; print('EasyMD imported successfully')"
easymd --help
```

## Method 3: Alternative Installation Script

Run the provided installation script:

```bash
python install_easymd.py
```

## Method 4: Development Mode

If you're still having issues, run EasyMD in development mode:

```bash
# From the EasyMD directory
python -m EasyMD --help
python -m EasyMD info --help
python -m EasyMD analyze --help
```

## Common Issues and Solutions

### Issue: "ModuleNotFoundError: No module named 'EasyMD.sysGenerator.sysGenerator'"

**Solution 1**: Reinstall in development mode
```bash
pip uninstall EasyMD
pip install -e .
```

**Solution 2**: Check Python path
```bash
python -c "import sys; print('\\n'.join(sys.path))"
```
Make sure the EasyMD/src directory is in the path.

**Solution 3**: Use absolute imports
If you're still having issues, the package has been updated to use absolute imports throughout.

### Issue: "easymd command not found"

**Solution**: The console script wasn't installed properly
```bash
pip install -e .
# or
python -m EasyMD --help  # Use module mode instead
```

### Issue: Import errors in conda environment

**Solution**: Ensure all dependencies are installed
```bash
conda install -c conda-forge setuptools wheel pip
pip install -e .
```

## Verification Steps

After installation, verify everything works:

```bash
# Test basic import
python -c "import EasyMD; print('✅ EasyMD imported successfully')"

# Test command line tools
easymd --help
easymd info --help
easymd analyze --help

# Test specific functionality
python -c "from EasyMD.sysGenerator.sysGenerator import SysGenerator; print('✅ SysGenerator imported')"
python -c "from EasyMD.argManager.manager import ArgManager; print('✅ ArgManager imported')"
```

## Environment Information

If you're still having issues, please provide:

```bash
# Python and environment info
python --version
conda --version
conda list | grep -E "(numpy|pandas|openmm|mdtraj)"

# Package installation info
pip show EasyMD
python -c "import EasyMD; print(EasyMD.__file__)"
```

## Contact

If none of these solutions work, please open an issue with:
1. Your operating system
2. Python version
3. Conda version
4. The exact error message
5. Output from the verification steps above