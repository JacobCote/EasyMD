# Utils Component

## Overview
The Utils component provides essential utility functions and helper tools that support the core EasyMD functionality. It includes platform detection, molecular manipulation tools, charge calculation utilities, and specialized algorithms like simulated annealing. These utilities ensure optimal performance and provide common functionality across all EasyMD components.

## Architecture

### Component Structure
```
utils/
├── __init__.py              # Module initialization and exports
├── utils.py                 # Core utility functions
├── openbabel_charge.py      # Ligand charge calculation
└── simulated_annealing.py   # Simulated annealing implementation
```

### Key Modules
- **utils.py**: Platform detection, molecular manipulation, file operations
- **openbabel_charge.py**: OpenBabel-based charge calculation for ligands
- **simulated_annealing.py**: Specialized annealing simulation protocol

## Core Utilities (utils.py)

### 1. Platform Detection and Optimization
```python
def get_platform():
    """
    Automatically detect and configure the optimal OpenMM platform.
    
    Returns:
        Platform: The fastest available OpenMM platform with optimized settings
    """
    os_platform = os.getenv('PLATFORM')
    if os_platform:
        platform = Platform.getPlatformByName(os_platform)
    else:
        # Automatically select fastest platform
        speed = 0
        for i in range(Platform.getNumPlatforms()):
            p = Platform.getPlatform(i)
            if p.getSpeed() > speed:
                platform = p
                speed = p.getSpeed()
    
    # Optimize GPU platforms
    if platform.getName() in ['CUDA', 'OpenCL']:
        platform.setPropertyDefaultValue('Precision', 'mixed')
        print('Set precision for platform', platform.getName(), 'to mixed')
    
    return platform
```

**Features:**
- **Automatic Detection**: Finds fastest available platform (CPU, CUDA, OpenCL)
- **Environment Override**: Respects PLATFORM environment variable
- **GPU Optimization**: Automatically sets mixed precision for GPU platforms
- **Performance Reporting**: Displays selected platform and optimizations

**Platform Hierarchy:**
1. **CUDA**: NVIDIA GPU acceleration (fastest)
2. **OpenCL**: Cross-platform GPU acceleration
3. **CPU**: Multi-threaded CPU execution (fallback)

### 2. Molecular Manipulation Tools
```python
def insert_molecule_and_remove_clashes(
    topology: Topology,
    insert: Molecule,
    radius: Quantity = 1.5 * unit.angstrom,
    keep: list[Molecule] = [],
) -> Topology:
    """
    Add a molecule to a topology while removing clashing molecules.
    
    Parameters:
        topology: The topology to insert a molecule into
        insert: The molecule to insert
        radius: Clash detection radius
        keep: Molecules to preserve even if clashing
    
    Returns:
        Topology: New topology with inserted molecule and clashes removed
    """
```

**Functionality:**
- **Clash Detection**: Identifies overlapping molecules within specified radius
- **Selective Removal**: Removes clashing molecules while preserving important ones
- **Topology Preservation**: Maintains topology integrity and box vectors
- **Coordinate Analysis**: Uses 3D distance calculations for clash detection

**Applications:**
- **Ligand Insertion**: Add ligands to protein structures
- **Solvent Optimization**: Remove problematic solvent molecules
- **Structure Cleaning**: Eliminate overlapping artifacts

### 3. File Operations and Formatting
```python
def writeFooter(topology, file):
    """Write out the footer for a PDB file."""
    # Handles proper PDB file termination
    # Ensures format compliance
```

**Additional Utilities:**
- **PDB Writing**: Specialized PDB output functions
- **Format Validation**: Ensure proper file formats
- **Index Management**: Handle atom and residue indexing

## OpenBabel Charge Calculation (openbabel_charge.py)

### Charge Calculation Engine
```python
def get_charges(input_file, output_file, file_type_input, file_type_output):
    """
    Calculate partial charges for ligands using OpenBabel and MMFF94.
    
    Parameters:
        input_file (str): Path to input ligand file
        output_file (str): Path to output file with charges
        file_type_input (str): Input file format (pdb, mol2, etc.)
        file_type_output (str): Output file format (sdf, mol2, etc.)
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(file_type_input, file_type_output)
    
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_file)
    
    # Calculate MMFF94 charges
    charge_model = openbabel.OBChargeModel.FindType("MMFF94")
    charge_model.ComputeCharges(mol)
    
    obConversion.WriteFile(mol, output_file)
```

**Features:**
- **MMFF94 Force Field**: Industry-standard charge calculation method
- **Format Flexibility**: Supports multiple input/output formats
- **Automatic Processing**: Seamless integration with ligand preparation
- **Quality Assurance**: Reliable charge assignment for MD simulations

**Supported Formats:**
- **Input**: PDB, MOL2, SDF, MOL
- **Output**: SDF, MOL2, PDB (with charges)

**Integration Points:**
- **SysGenerator**: Used during ligand parameterization
- **Complex Preparation**: Essential for protein-ligand simulations
- **Force Field Assignment**: Provides charges for OpenFF/GAFF

### Command-Line Interface
```bash
# Calculate charges for ligand
python -m EasyMD.utils.openbabel_charge -i ligand.pdb -f pdb -o ligand_charged.sdf
```

## Simulated Annealing Implementation (simulated_annealing.py)

### Annealing Protocol
```python
def simulated_annealing(modeller, system, temperature, out_dir, step_size, 
                       friction_coeff, reporting_interval, equilibration_steps):
    """
    Perform simulated annealing molecular dynamics simulation.
    
    This function implements a temperature-ramping protocol to help systems
    escape local energy minima and explore conformational space more effectively.
    
    Parameters:
        modeller: OpenMM Modeller with system topology and positions
        system: OpenMM System with force field parameters
        temperature: Starting temperature in Kelvin
        out_dir: Output directory for results
        step_size: Integration time step in picoseconds
        friction_coeff: Langevin friction coefficient
        reporting_interval: Frequency of output reporting
        equilibration_steps: Initial equilibration steps
    """
```

### Annealing Algorithm
```python
# Annealing protocol: 2000 cycles of temperature ramping
for i in range(2000):
    # Increase temperature by 0.1K per cycle
    integrator.setTemperature(temperature + (0.1 * i) * unit.kelvin)
    simulation.step(100)  # 100 steps per temperature point
```

**Protocol Details:**
- **Cycles**: 2000 annealing cycles
- **Steps per Cycle**: 100 MD steps
- **Temperature Increment**: 0.1K per cycle
- **Total Temperature Range**: Base temperature + 200K
- **Total Steps**: 200,000 steps

**Applications:**
- **Protein Folding**: Help proteins find native conformations
- **Conformational Sampling**: Explore multiple stable states
- **Structure Optimization**: Escape local energy minima
- **Refinement**: Improve initial structure quality

### Error Handling
```python
try:
    simulation.minimizeEnergy(maxIterations=10000)
    print('✓ Energy minimization completed successfully')
except Exception as e:
    if "NaN" in str(e) or "coordinate is NaN" in str(e):
        print("❌ Energy minimization failed due to NaN coordinates")
        print("Troubleshooting suggestions:")
        print("  1. Check input PDB structure for overlapping atoms")
        print("  2. Verify ligand coordinates are reasonable")
        print("  3. Try using --keep-water flag")
        raise RuntimeError(f"Energy minimization failed: {e}")
```

**Error Detection:**
- **NaN Coordinates**: Detects coordinate failures
- **Minimization Issues**: Identifies energy minimization problems
- **User Guidance**: Provides specific troubleshooting steps
- **Graceful Failure**: Clean error reporting and exit

## Integration Across EasyMD

### 1. Platform Optimization
```python
# Used in main entry point
platform = utils.get_platform()
```
- **Automatic Selection**: Chooses optimal platform for each system
- **Performance Optimization**: Maximizes simulation speed
- **Hardware Utilization**: Efficiently uses available resources

### 2. Ligand Processing
```python
# Used in SysGenerator for ligand parameterization
from EasyMD.utils.openbabel_charge import get_charges
get_charges(f'{outdir}/ligand.pdb', f'{outdir}/ligand.sdf', 'pdb', 'sdf')
```
- **Charge Assignment**: Essential for ligand force field parameterization
- **Format Conversion**: Converts between molecular file formats
- **Quality Control**: Ensures proper ligand preparation

### 3. Specialized Protocols
```python
# Used in AnnealingRunner
from EasyMD.utils.simulated_annealing import simulated_annealing
```
- **Alternative Protocols**: Provides specialized simulation methods
- **Enhanced Sampling**: Improves conformational exploration
- **Research Applications**: Supports advanced simulation techniques

## Performance Considerations

### 1. Platform Optimization
- **GPU Acceleration**: Automatic detection and configuration of CUDA/OpenCL
- **Mixed Precision**: Optimal balance of speed and accuracy
- **Memory Management**: Efficient resource utilization

### 2. Computational Efficiency
- **Vectorized Operations**: NumPy-based calculations where applicable
- **Optimized Algorithms**: Efficient molecular manipulation
- **Memory Conservation**: Minimal memory footprint for utilities

### 3. I/O Optimization
- **Format Handling**: Efficient file format conversions
- **Batch Processing**: Optimized for multiple operations
- **Error Recovery**: Robust handling of file I/O issues

## Quality Assurance

### 1. Input Validation
- **Format Verification**: Ensures proper file formats
- **Parameter Validation**: Checks reasonable parameter ranges
- **Dependency Checking**: Verifies required libraries

### 2. Error Handling
- **Graceful Degradation**: Continues operation when possible
- **Clear Error Messages**: Specific guidance for common issues
- **Recovery Mechanisms**: Fallback options for failures

### 3. Testing Support
- **Mock Objects**: Handles test environments gracefully
- **Validation Functions**: Supports comprehensive testing
- **Debug Information**: Detailed logging for troubleshooting

## Dependencies and Requirements

### 1. Core Dependencies
- **OpenMM**: Platform detection and simulation engine
- **OpenBabel**: Chemical structure processing and charge calculation
- **OpenFF**: Molecular topology and force field tools
- **NumPy**: Numerical computations

### 2. Optional Dependencies
- **CUDA**: GPU acceleration (if available)
- **OpenCL**: Cross-platform GPU support
- **MDTraj**: Trajectory analysis support

## Usage Examples

### Platform Detection
```python
from EasyMD.utils import get_platform
platform = get_platform()
# Automatically selects fastest available platform
```

### Charge Calculation
```bash
# Command line usage
python -m EasyMD.utils.openbabel_charge -i ligand.pdb -f pdb -o ligand.sdf

# Programmatic usage
from EasyMD.utils import get_charges
get_charges('ligand.pdb', 'ligand.sdf', 'pdb', 'sdf')
```

### Molecular Manipulation
```python
from EasyMD.utils import insert_molecule_and_remove_clashes
new_topology = insert_molecule_and_remove_clashes(
    topology, ligand_molecule, radius=2.0*unit.angstrom
)
```

## Benefits

1. **Performance Optimization**: Automatic platform detection and configuration
2. **Chemical Accuracy**: Reliable charge calculation for ligands
3. **Specialized Methods**: Advanced simulation protocols like annealing
4. **Code Reusability**: Common functions shared across components
5. **Error Resilience**: Robust error handling and recovery
6. **Hardware Utilization**: Optimal use of available computational resources
7. **Format Flexibility**: Support for multiple molecular file formats

The Utils component provides the essential foundation that enables EasyMD's high performance, reliability, and flexibility, ensuring that all components can operate efficiently and handle diverse molecular systems with appropriate computational resources.