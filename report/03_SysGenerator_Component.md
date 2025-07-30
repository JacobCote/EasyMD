# SysGenerator Component

## Overview
The SysGenerator component is responsible for preparing molecular systems for MD simulations. It handles protein structure fixing, ligand parameterization, solvation, and force field assignment. The component supports both new simulations and restart from previous states, with options for explicit solvation or implicit solvent (GBIS) models.

## Architecture

### Class Structure
```python
class SysGenerator:
    """
    System generator for molecular dynamics simulations.
    
    Attributes:
        config: Configuration object containing simulation parameters
        forcefield_kwargs: Force field configuration parameters
        last_state: State tracking for restart simulations
        modeller: OpenMM Modeller object with system topology and positions
        system: OpenMM System object with force field parameters
    """
```

### Key Features
1. **Automated Structure Fixing**: PDBFixer integration for missing residues/atoms
2. **Dual System Support**: Protein-only and protein-ligand complexes
3. **Flexible Solvation**: Explicit and implicit solvent models
4. **Restart Capability**: Resume from previous simulation states
5. **Comprehensive Validation**: Structure and coordinate validation
6. **Missing Residue Handling**: Advanced strategies for incomplete structures

## System Preparation Workflow

### 1. Initialization and Setup
```python
def __init__(self, config):
    """Initialize the SysGenerator and prepare the molecular system."""
    self.config = config
    self.forcefield_kwargs = {
        'constraints': app.HBonds, 
        'rigidWater': True, 
        'removeCMMotion': False, 
        'hydrogenMass': 4*unit.amu
    }
    
    if self.config.restart:
        self.modeller, self.system = self._restartSetup()
    else:
        self.modeller, self.system = self._setup()
```

### 2. New System Setup
The `_setup()` method handles new simulations:
- **Output Directory Management**: Auto-generates unique directories
- **Configuration Saving**: Stores setup parameters in YAML format
- **System Type Detection**: Routes to protein-only or complex preparation

### 3. Restart System Setup
The `_restartSetup()` method handles simulation continuation:
- **Configuration Loading**: Reads previous setup from `restart_setup.yml`
- **Ligand Detection**: Checks for `ligand.sdf` file
- **State Restoration**: Loads topology and positions from restart files

## Protein-Only System Preparation

### Structure Fixing Process
```python
def _prep_prot(self, pdb_in, ...):
    """Prepare protein-only system for molecular dynamics simulation."""
    
    # 1. Load and analyze structure
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    
    # 2. Fix structure issues
    fixer.replaceNonstandardResidues()
    self._handle_missing_residues(fixer)
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    
    # 3. Clean and prepare system
    # Remove unwanted molecules, add solvation if requested
```

### Key Processing Steps
1. **Structure Analysis**: Identify missing residues, atoms, and non-standard residues
2. **Issue Reporting**: Detailed console output of found problems
3. **Structure Repair**: Replace non-standard residues, add missing components
4. **Molecule Removal**: Remove specified artifacts (DMS, water, etc.)
5. **Solvation**: Add explicit solvent or configure implicit solvent
6. **System Generation**: Create OpenMM System with force field parameters

## Protein-Ligand Complex Preparation

### Complex Processing Workflow
```python
def _prep_complex(self, pdb_in, lig_name, ...):
    """Prepare protein-ligand complex for molecular dynamics simulation."""
    
    # 1. Fix protein structure (same as protein-only)
    # 2. Extract and validate ligand
    # 3. Process ligand separately
    # 4. Recombine with proper force fields
    # 5. Add solvation if requested
```

### Ligand Processing Steps
1. **Ligand Extraction**: Separate ligand from protein structure
2. **Structure Validation**: Check for coordinate issues and structural problems
3. **Charge Calculation**: Use OpenBabel for partial charge assignment
4. **Format Conversion**: Convert to SDF format for force field parameterization
5. **Topology Integration**: Add ligand topology to protein system

### Coordinate Validation
```python
def _validate_coordinates(self, positions, molecule_name="molecule"):
    """Validate that coordinates are finite and not NaN."""
    # Check for NaN values
    # Check for infinite values
    # Validate coordinate format
    # Report validation results
```

## Missing Residue Handling

### Advanced Missing Residue Strategies
The component provides sophisticated missing residue handling:

#### 1. Strategy Options
- **auto**: Automatic detection and intelligent handling
- **none**: Skip all missing residues
- **non-terminal**: Add only internal missing residues
- **terminal-only**: Add only N/C-terminal residues
- **all**: Add all detected missing residues

#### 2. Filtering Methods
```python
def _handle_missing_residues(self, fixer):
    """Handle missing residues based on user configuration."""
    missing_residues_strategy = getattr(self.config, 'missing_residues', 'auto')
    
    if missing_residues_strategy == "none":
        fixer.missingResidues = {}
    elif missing_residues_strategy == "terminal-only":
        self._filter_terminal_residues(fixer, max_terminal_residues, terminal_residue_types)
    # ... other strategies
```

#### 3. Conservative Filtering
- **Terminal Limits**: Restrict number of terminal residues added
- **Loop Detection**: Skip likely loop regions (>3 consecutive missing)
- **Type Validation**: Ensure appropriate terminal residue types

## Solvation Systems

### Explicit Solvation
```python
# Add explicit solvent with ions
modeller.addSolvent(
    system_generator.forcefield, 
    model=water_model, 
    padding=padding * unit.angstroms,
    positiveIon=positive_ion, 
    negativeIon=negative_ion,
    ionicStrength=ionic_strength * unit.molar, 
    neutralize=not no_neutralize
)
```

**Features:**
- Multiple water models (TIP3P, SPC/E, TIP4P-Ew, etc.)
- Automatic neutralization with ions
- Configurable ionic strength
- Customizable box padding

### Implicit Solvation (GBIS)
```python
# Configure implicit solvent system
system_generator = SystemGenerator(
    forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
    forcefield_kwargs=forcefield_kwargs,
    nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
)
```

**Benefits:**
- Faster simulation performance
- No periodic boundary conditions
- Reduced system size
- Suitable for conformational studies

## Force Field Integration

### SystemGenerator Configuration
The component uses OpenMM-ForceFields SystemGenerator for flexible force field assignment:

#### Protein-Only Systems
```python
system_generator = SystemGenerator(
    forcefields=[protein_force_field, water_force_field],
    forcefield_kwargs=forcefield_kwargs
)
```

#### Protein-Ligand Systems
```python
system_generator = SystemGenerator(
    forcefields=[protein_force_field, water_force_field],
    small_molecule_forcefield=ligand_force_field,
    molecules=[ligand_mol],
    forcefield_kwargs=forcefield_kwargs
)
```

### Supported Force Fields
- **Protein**: AMBER (14, 99sb-ildn, 03), CHARMM (36, 27)
- **Ligand**: OpenFF (2.2.0, 2.1.0, 2.0.0), GAFF (2.11, 1.81)
- **Water**: TIP3P, SPC/E, TIP4P-Ew, TIP5P

## Restart Functionality

### Restart File Management
Required files for restart:
- `restart_setup.yml`: Original configuration parameters
- `restart_model.pdb`: System topology and positions
- `last_state.xml`: OpenMM state information
- `ligand.sdf`: Ligand structure (if present)

### Restart System Recreation
```python
def _prep_restart_ligand(self, setup, outdir, forcefield_kwargs):
    """Prepare system for restart simulation with ligand present."""
    pdb = PDBFile(f'{outdir}/restart_model.pdb')
    ligand_mol = Molecule.from_file(f'{outdir}/ligand.sdf')
    
    # Recreate system with same parameters
    system_generator = SystemGenerator(...)
    system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
```

## Quality Assurance

### Structure Validation
1. **Coordinate Validation**: Check for NaN and infinite values
2. **Ligand Structure Validation**: Verify reasonable coordinate ranges
3. **Topology Consistency**: Ensure proper atom connectivity
4. **Force Field Compatibility**: Validate parameter assignment

### Error Handling
- **Graceful Degradation**: Continue with warnings for non-critical issues
- **Clear Error Messages**: Specific guidance for structure problems
- **Validation Checkpoints**: Multiple validation stages throughout preparation

### Output Files Generated
- `prot_receptor.pdb`: Cleaned protein structure
- `ligand.pdb`: Extracted ligand structure
- `ligand.sdf`: Ligand with charges and parameters
- `complex.pdb`: Combined protein-ligand system
- `solvated_complex.pdb`: Final solvated system
- `restart_model.pdb`: System ready for simulation
- `setup.yml`: Configuration backup

## Integration Points

### With ArgManager
- Receives validated configuration parameters
- Uses missing residue handling preferences
- Applies solvation and force field choices

### With SimRunner
- Provides prepared `modeller` and `system` objects
- Ensures compatibility with simulation parameters
- Maintains restart capability

### With Utils
- Uses utility functions for file operations
- Integrates with OpenBabel for ligand processing
- Leverages platform-specific optimizations

## Performance Considerations

### Optimization Strategies
1. **Efficient Structure Processing**: Minimize redundant operations
2. **Memory Management**: Clean up intermediate structures
3. **Validation Caching**: Avoid repeated coordinate checks
4. **Force Field Reuse**: Cache SystemGenerator instances

### Scalability
- **Large Systems**: Handles proteins with thousands of residues
- **Complex Ligands**: Supports multi-residue ligands
- **Multiple Chains**: Processes multi-chain protein complexes
- **Extensive Solvation**: Manages large solvent boxes

The SysGenerator component provides a robust, flexible foundation for MD system preparation, ensuring that structures are properly validated, parameterized, and ready for high-quality molecular dynamics simulations.