# EasyMD Main Entry Points

## Overview
The EasyMD program provides multiple entry points for different functionalities, all accessible through the main `__main__.py` module. The program supports three primary modes: simulation, analysis, and structure information.

## Entry Point Architecture

### 1. Main Entry Point (`__main__.py`)
The primary entry point that handles command routing and execution.

#### Key Functions:
- **`main()`**: Primary entry point for simulation mode
- **`run_analysis()`**: Entry point for trajectory analysis
- **`run_info()`**: Entry point for structure information analysis

#### Command Routing Logic:
```python
# Analysis mode
if len(sys.argv) > 1 and sys.argv[1] == 'analyze':
    run_analysis()

# Info mode  
if len(sys.argv) > 1 and sys.argv[1] == 'info':
    run_info()

# Default simulation mode
main()
```

### 2. Package Initialization (`__init__.py`)
Defines the core module imports and public API:
```python
from EasyMD.simRunner.simRunner import SimRunner
from EasyMD.sysGenerator import SysGenerator
from EasyMD.argManager.manager import ArgManager
```

### 3. Installation Validator (`validate_installation.py`)
Comprehensive installation testing and validation script.

## Command Line Interface

### Simulation Mode (Default)
```bash
python -m EasyMD --protein protein.pdb --steps 5000000 --solvate
```

**Features:**
- Protein-only and protein-ligand simulations
- Explicit and implicit solvation
- Time-based or step-based duration
- Configuration file support
- Restart capabilities

**Example Commands:**
```bash
# Basic protein simulation
python -m EasyMD --protein protein.pdb --steps 5000000 --solvate

# Protein-ligand complex
python -m EasyMD --protein complex.pdb --ligand ATP --steps 10000000 --solvate

# Implicit solvent simulation
python -m EasyMD --protein protein.pdb --steps 1000000 --GBIS

# Time-based simulation
python -m EasyMD --protein protein.pdb --clock 50ns --solvate

# Restart simulation
python -m EasyMD --restart out_0/
```

### Analysis Mode
```bash
python -m EasyMD analyze out_0/ --rmsd --rmsf
```

**Features:**
- RMSD and RMSF calculations
- Inter-chain distance analysis
- Radius of gyration
- Secondary structure analysis
- Custom analysis scripts

**Example Commands:**
```bash
# Basic analysis
python -m EasyMD analyze out_0/ --rmsd --rmsf

# All analyses with data export
python -m EasyMD analyze out_0/ --all --save-data

# Custom parameters
python -m EasyMD analyze out_0/ --rmsd --atom-selection backbone --reference-frame 10
```

### Info Mode
```bash
python -m EasyMD info protein.pdb
```

**Features:**
- Chain composition analysis
- Missing residue detection (PDBFixer-based)
- Ligand identification
- Water molecule analysis
- Disulfide bond prediction
- Metal ion detection

**Example Commands:**
```bash
# Basic structure analysis
python -m EasyMD info protein.pdb

# Custom output format
python -m EasyMD info protein.pdb --format summary

# JSON output
python -m EasyMD info protein.pdb --format json

# Terminal only (no file output)
python -m EasyMD info protein.pdb --no-file
```

## Error Handling and User Experience

### Graceful Error Handling
All entry points implement comprehensive error handling:
```python
try:
    # Main functionality
    runner.run()
except KeyboardInterrupt:
    print("\n\n⚠️  Analysis interrupted by user")
    sys.exit(1)
except Exception as e:
    print(f"\n❌ Analysis failed: {str(e)}")
    sys.exit(1)
```

### User-Friendly Help System
Each mode provides detailed help with:
- Command descriptions
- Usage examples
- Parameter explanations
- Output file descriptions
- Common use cases

## Installation Validation

### Validation Script Features
The `validate_installation.py` script provides:

1. **Module Import Testing**
   - Tests all core EasyMD modules
   - Verifies proper installation

2. **Dependency Checking**
   - OpenMM availability
   - OpenFF Toolkit
   - PyYAML, MDTraj, PDBFixer
   - Reports missing dependencies

3. **Validation System Testing**
   - Tests argument parsing
   - Validates error handling
   - Checks parameter validation

4. **Comprehensive Reporting**
   - Clear pass/fail indicators
   - Troubleshooting guidance
   - Next steps recommendations

### Running Validation
```bash
python -m EasyMD.validate_installation
```

## Output Files and Structure

### Simulation Output
- `output_traj_*.dcd`: Trajectory files
- `log.txt`: Energy and temperature logs
- `last_state_*.pdb`: Final structures
- `topology.pkl`: System topology
- `restart_setup.yml`: Restart configuration

### Analysis Output
- Analysis plots: `rmsd.png`, `rmsf.png`, etc.
- Data files: `rmsd_data.csv`, `rmsf_data.csv`
- Terminal summaries

### Info Output
- Detailed reports: `<pdb_name>.info`
- Custom output files
- Terminal display with formatting

## Integration Points

The main entry points integrate with all EasyMD components:
- **ArgManager**: Command-line parsing and validation
- **SysGenerator**: System preparation and validation
- **SimRunner**: Simulation execution
- **AnalysisManager/Runner**: Trajectory analysis
- **InfoManager/Runner**: Structure information
- **Utils**: Platform detection and utilities

This architecture provides a unified interface while maintaining modular functionality and clear separation of concerns.