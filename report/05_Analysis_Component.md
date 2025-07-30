# Analysis Component

## Overview
The Analysis component provides comprehensive trajectory analysis tools for molecular dynamics simulations. It supports multiple analysis types including RMSD, RMSF, inter-chain distances, radius of gyration, and secondary structure analysis, with flexible output options and robust data handling.

## Architecture

### Component Structure
```
analysis/
├── __init__.py           # Module initialization
├── analysisManager.py    # Argument parsing and validation
└── analysisRunner.py     # Analysis execution and data processing
```

### Key Classes
- **AnalysisManager**: Handles command-line arguments and validation
- **AnalysisRunner**: Performs trajectory analysis and generates outputs

## AnalysisManager - Argument Management

### Argument Categories

#### 1. Input/Output Parameters
```python
io_group.add_argument("trajectory_directory", type=str,
                     help="Path to directory containing trajectory files and topology")
io_group.add_argument("-o", "--output-dir", type=str, default=None,
                     help="Output directory for analysis results")
io_group.add_argument("--output-format", choices=["png", "pdf", "svg"], default="png",
                     help="Output format for plots")
```

#### 2. Analysis Types
- **RMSD**: Root Mean Square Deviation calculation
- **RMSF**: Root Mean Square Fluctuation analysis
- **Distances**: Inter-chain distance measurements
- **Radius of Gyration**: Compactness analysis
- **Secondary Structure**: Structural evolution analysis
- **All**: Comprehensive analysis suite

#### 3. Analysis Parameters
```python
params_group.add_argument("--reference-frame", type=int, default=0,
                         help="Reference frame for RMSD calculation")
params_group.add_argument("--atom-selection", choices=["all", "backbone", "ca", "heavy"], 
                         default="backbone", help="Atom selection for RMSD calculation")
params_group.add_argument("--skip-frames", type=int, default=1,
                         help="Skip every N frames for analysis")
```

#### 4. Output Options
- **Data Export**: CSV file generation
- **Plot Control**: Format, DPI, style options
- **Visualization**: Multiple matplotlib styles

### Validation System
Comprehensive validation ensures:
- **Directory Existence**: Trajectory directory validation
- **Required Files**: Topology and trajectory file checks
- **Parameter Ranges**: Reasonable analysis parameter values
- **Output Permissions**: Write access verification

## AnalysisRunner - Analysis Execution

### Initialization and Setup
```python
class AnalysisRunner:
    """
    Performs trajectory analysis for molecular dynamics simulations.
    
    Attributes:
        config: Configuration object containing analysis parameters
        trajectory_dir (Path): Path to trajectory directory
        output_dir (Path): Path to output directory
        topology: MDTraj topology object
        trajectory: Combined MDTraj trajectory object
    """
```

### Data Loading Workflow

#### 1. Topology Loading
```python
def _load_trajectory_data(self):
    """Load trajectory files and topology."""
    # Load OpenMM topology from pickle file
    with open(topology_file, 'rb') as f:
        topo_openmm = pickle.load(f)
    self.topology = md.Topology.from_openmm(topo_openmm)
    
    # Find and load all trajectory files
    dcd_files = sorted([f for f in os.listdir(self.trajectory_dir) 
                       if f.startswith('output_traj_') and f.endswith('.dcd')])
```

#### 2. Trajectory Combination
- **Multi-file Support**: Automatically combines multiple trajectory files
- **Frame Continuity**: Maintains proper temporal ordering
- **Memory Efficiency**: Optimized loading for large trajectories

#### 3. Frame Selection
```python
def _apply_frame_selection(self):
    """Apply frame selection based on start, end, and skip parameters."""
    frame_indices = list(range(start, min(end, self.trajectory.n_frames), skip))
    self.trajectory = self.trajectory[frame_indices]
```

## Analysis Methods

### 1. RMSD Analysis
```python
def _calculate_rmsd(self) -> Dict[str, Any]:
    """Calculate Root Mean Square Deviation (RMSD)."""
    # Get atom indices based on selection
    atom_indices = self._get_atom_indices(self.config.atom_selection)
    
    # Calculate RMSD against reference frame
    rmsd_values = md.rmsd(self.trajectory, self.trajectory, reference_frame, 
                         atom_indices=atom_indices)
    
    # Convert to Angstroms and return comprehensive results
    rmsd_values *= 10  # nm to Angstroms
```

**Features:**
- **Atom Selection**: All, backbone, CA, or heavy atoms
- **Reference Frame**: Configurable reference structure
- **Statistical Analysis**: Mean, std, min, max calculations
- **Unit Conversion**: Automatic nm to Angstrom conversion

**Output Data:**
- RMSD values over time
- Statistical summaries
- Atom selection metadata
- Reference frame information

### 2. RMSF Analysis
```python
def _calculate_rmsf(self) -> Dict[str, Any]:
    """Calculate Root Mean Square Fluctuation (RMSF)."""
    # Typically use CA atoms for proteins
    ca_indices = [atom.index for atom in self.trajectory.topology.atoms 
                 if atom.name == 'CA']
    
    # Calculate RMSF per residue
    rmsf_values = md.rmsf(self.trajectory, self.trajectory, reference_frame, 
                         atom_indices=ca_indices)
```

**Features:**
- **Per-Residue Analysis**: Individual residue flexibility
- **CA Atom Focus**: Standard protein analysis approach
- **Residue Mapping**: Links RMSF values to specific residues
- **Flexibility Profiling**: Identifies flexible and rigid regions

**Applications:**
- **Binding Site Analysis**: Identify flexible binding regions
- **Allosteric Studies**: Map conformational changes
- **Stability Assessment**: Evaluate structural stability

### 3. Inter-Chain Distance Analysis
```python
def _calculate_distances(self) -> Dict[str, Any]:
    """Calculate distances between centers of mass of chains."""
    # Get chains and calculate center of mass distances
    for i, chain1 in enumerate(chains):
        for j, chain2 in enumerate(chains[i+1:], i+1):
            # Calculate center of mass for each chain
            com1 = md.compute_center_of_mass(trajectory[frame].atom_slice(chain1_atoms))
            com2 = md.compute_center_of_mass(trajectory[frame].atom_slice(chain2_atoms))
            dist = np.linalg.norm(com1 - com2) * 10  # Convert to Angstroms
```

**Features:**
- **Multi-Chain Support**: Handles protein complexes
- **Center of Mass**: Robust distance measurements
- **Time Evolution**: Tracks distance changes over simulation
- **Complex Analysis**: Suitable for protein-protein interactions

**Applications:**
- **Complex Stability**: Monitor protein-protein interfaces
- **Conformational Changes**: Track domain movements
- **Binding Analysis**: Study ligand-protein interactions

### 4. Radius of Gyration
```python
def _calculate_radius_gyration(self) -> Dict[str, Any]:
    """Calculate radius of gyration."""
    # Calculate compactness measure
    rg_values = md.compute_rg(self.trajectory)
    rg_values *= 10  # Convert to Angstroms
```

**Features:**
- **Compactness Measure**: Overall protein compactness
- **Folding Analysis**: Monitor folding/unfolding events
- **Size Characterization**: Quantify protein size changes

**Applications:**
- **Folding Studies**: Track protein folding pathways
- **Denaturation**: Monitor unfolding processes
- **Conformational States**: Identify compact vs extended states

### 5. Secondary Structure Analysis
```python
def _analyze_secondary_structure(self) -> Dict[str, Any]:
    """Analyze secondary structure evolution."""
    # Use MDTraj's DSSP implementation
    ss_assignments = md.compute_dssp(self.trajectory)
```

**Features:**
- **DSSP Algorithm**: Standard secondary structure assignment
- **Time Evolution**: Track structural changes over time
- **Per-Residue Analysis**: Detailed structural mapping
- **Statistical Analysis**: Secondary structure content quantification

## Output Generation

### 1. Plot Generation
```python
def _generate_outputs(self, results):
    """Generate plots and save data files."""
    if not self.config.no_plots:
        self._create_plots(results)
    
    if self.config.save_data:
        self._save_data_files(results)
```

**Plot Types:**
- **RMSD Plots**: Time series with statistical annotations
- **RMSF Plots**: Per-residue flexibility profiles
- **Distance Plots**: Inter-chain distance evolution
- **Radius of Gyration**: Compactness over time
- **Secondary Structure**: Structural evolution heatmaps

### 2. Data Export
**CSV Files:**
- `rmsd_data.csv`: RMSD values and metadata
- `rmsf_data.csv`: Per-residue RMSF values
- `distances_data.csv`: Inter-chain distances
- `radius_gyration_data.csv`: Compactness measurements

**Plot Files:**
- High-resolution plots in PNG, PDF, or SVG format
- Customizable DPI and styling options
- Professional publication-ready figures

### 3. Statistical Summaries
Each analysis includes comprehensive statistics:
- **Mean and Standard Deviation**
- **Minimum and Maximum Values**
- **Percentile Analysis**
- **Trend Analysis**

## Integration Features

### 1. MDTraj Integration
- **Efficient Loading**: Optimized trajectory handling
- **Format Support**: Multiple trajectory formats
- **Analysis Tools**: Comprehensive analysis functions
- **Memory Management**: Efficient large trajectory processing

### 2. Matplotlib Integration
- **Flexible Plotting**: Multiple plot styles and formats
- **Publication Quality**: High-resolution output options
- **Customization**: Style and formatting control
- **Interactive Options**: Support for different backends

### 3. Pandas Integration
- **Data Management**: Structured data handling
- **Export Options**: Multiple file formats
- **Statistical Analysis**: Built-in statistical functions
- **Data Manipulation**: Flexible data processing

## Error Handling and Validation

### 1. Input Validation
- **File Existence**: Verify trajectory and topology files
- **Format Validation**: Check file formats and integrity
- **Parameter Validation**: Ensure reasonable analysis parameters

### 2. Runtime Error Handling
```python
try:
    results['rmsd'] = self._calculate_rmsd()
except Exception as e:
    print(f"❌ RMSD calculation failed: {str(e)}")
    return None
```

### 3. Data Quality Checks
- **Trajectory Continuity**: Verify frame consistency
- **Coordinate Validation**: Check for NaN or infinite values
- **Topology Consistency**: Ensure atom count consistency

## Performance Optimization

### 1. Memory Management
- **Efficient Loading**: Stream processing for large files
- **Memory Cleanup**: Proper object disposal
- **Chunked Processing**: Handle large trajectories in segments

### 2. Computational Efficiency
- **Vectorized Operations**: NumPy-based calculations
- **Parallel Processing**: Multi-core utilization where possible
- **Optimized Algorithms**: Efficient analysis implementations

### 3. I/O Optimization
- **Buffered Operations**: Efficient file reading/writing
- **Compressed Support**: Handle gzipped trajectories
- **Format Optimization**: Use efficient data formats

## Usage Examples

### Basic Analysis
```bash
# Perform RMSD and RMSF analysis
python -m EasyMD analyze out_0/ --rmsd --rmsf

# All analyses with data export
python -m EasyMD analyze out_0/ --all --save-data

# Custom parameters
python -m EasyMD analyze out_0/ --rmsd --atom-selection backbone --reference-frame 10
```

### Advanced Options
```bash
# High-quality plots with custom styling
python -m EasyMD analyze out_0/ --all --output-format pdf --dpi 600 --plot-style seaborn

# Frame selection and skipping
python -m EasyMD analyze out_0/ --rmsd --start-frame 1000 --end-frame 5000 --skip-frames 5
```

## Benefits

1. **Comprehensive Analysis**: Complete suite of trajectory analysis tools
2. **User-Friendly Interface**: Simple command-line interface with clear options
3. **Flexible Output**: Multiple formats and customization options
4. **Robust Error Handling**: Clear error messages and troubleshooting guidance
5. **Performance Optimized**: Efficient handling of large trajectories
6. **Publication Ready**: High-quality plots and data export
7. **Extensible Design**: Easy to add new analysis methods

The Analysis component provides researchers with powerful, easy-to-use tools for extracting meaningful insights from molecular dynamics simulations, supporting both routine analysis and advanced research applications.