# Restart Component

## Overview
The Restart component enables seamless continuation of molecular dynamics simulations from previously saved states. It provides robust functionality for loading saved configurations, recreating system topologies, and resuming simulations with complete state preservation, supporting both protein-only and protein-ligand complex systems.

## Architecture

### Component Structure
```
restart/
├── __init__.py      # Module initialization
└── restarter.py     # Restart functionality implementation
```

### Key Classes
- **Restarter**: Handles simulation restart from saved states

## Restarter Class

### Class Structure
```python
class Restarter:
    """
    Handles restarting molecular dynamics simulations from previously saved states.
    
    This class provides functionality to restart MD simulations by loading saved states,
    topologies, and configurations from previous runs. It supports both protein-only
    and protein-ligand complex simulations with or without explicit solvent.
    
    Attributes:
        setup (dict): Configuration parameters from the original simulation setup
        outdir (str): Output directory containing restart files
        forcefield_kwargs (dict): Force field parameters for system generation
    """
```

### Initialization
```python
def __init__(self, setup: dict, outdir: str, forcefield_kwargs: dict):
    """
    Initialize the Restarter with simulation setup parameters.
    
    Args:
        setup (dict): Configuration parameters from the original simulation including
                     force field specifications, solvation settings, and simulation parameters
        outdir (str): Path to the output directory containing restart files
        forcefield_kwargs (dict): Force field parameters for system generation
    """
    self.setup = setup
    self.outdir = outdir
    self.forcefield_kwargs = {
        'constraints': app.HBonds, 
        'rigidWater': True, 
        'removeCMMotion': False, 
        'hydrogenMass': 4*unit.amu
    }
```

## Required Restart Files

### 1. Core Restart Files
The restart functionality requires specific files generated during the original simulation:

#### restart_setup.yml
```yaml
# Original simulation configuration
pdb: 'restart_model.pdb'
reporting_interval: 1000
step_size: 0.002
friction_coeff: 1.0
temperature: 300
solvate: true
state: 'last_state.xml'
protein_force_field: 'amber14-all.xml'
ligand_force_field: 'openff-2.2.0'
water_force_field: 'amber/tip3p_standard.xml'
last_state: 0
```

#### restart_model.pdb
- Complete system topology and positions
- Protein structure with all atoms
- Ligand coordinates (if present)
- Solvent molecules (if solvated)

#### last_state.xml
- Complete OpenMM simulation state
- Particle positions and velocities
- Box vectors and periodic boundary conditions
- Integrator state and random number generator state

#### ligand.sdf (if ligand present)
- Ligand structure with partial charges
- Force field parameters
- Connectivity information

### 2. File Validation
The restart process validates the presence and integrity of required files:
```python
# File existence checks
required_files = [
    'restart_setup.yml',
    'restart_model.pdb', 
    'last_state.xml'
]

# Additional files for ligand systems
if ligand_present:
    required_files.append('ligand.sdf')
```

## Restart Methods

### 1. Protein-Ligand Complex Restart
```python
def prep_restart_ligand(self) -> Tuple[Modeller, SystemGenerator]:
    """
    Prepare system for restart simulation with ligand present.
    
    This method loads the restart model PDB file and ligand SDF file, then creates
    a new system with the appropriate force fields for protein-ligand simulations.
    It handles both solvated and implicit solvent (GBIS) systems.
    
    Returns:
        Tuple[Modeller, SystemGenerator]: A tuple containing:
            - Modeller: OpenMM Modeller object with topology and positions
            - SystemGenerator: OpenMM System object with forces and parameters
    """
    # Load restart files
    pdb = PDBFile(f'{outdir}/restart_model.pdb')
    ligand_mol = Molecule.from_file(f'{outdir}/ligand.sdf')
    
    # Extract force field parameters from setup
    protein_force_field = setup['protein_force_field']
    water_force_field = setup['water_force_field']
    ligand_force_field = setup['ligand_force_field']
    
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # Create system based on solvation type
    if setup['solvate']:
        print('Generating system with solvent...')
        system_generator = SystemGenerator(
            forcefields=[protein_force_field, water_force_field],
            small_molecule_forcefield=ligand_force_field,
            molecules=[ligand_mol],
            forcefield_kwargs=forcefield_kwargs
        )
        system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
    else:
        print('Generating implicit solvent system...')
        system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
            small_molecule_forcefield=ligand_force_field,
            molecules=[ligand_mol],
            forcefield_kwargs=forcefield_kwargs,
            nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
        )
        system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
    
    return modeller, system
```

### 2. Protein-Only Restart
```python
def prep_restart(self) -> Tuple[Modeller, SystemGenerator]:
    """
    Prepare system for restart simulation without ligand (protein-only).
    
    Similar to ligand restart but without small molecule force field components.
    Handles both explicit and implicit solvation systems.
    """
    pdb = PDBFile(f'{outdir}/restart_model.pdb')
    protein_force_field = setup['protein_force_field']
    water_force_field = setup['water_force_field']
    modeller = Modeller(pdb.topology, pdb.positions)
    
    if setup['solvate']:
        system_generator = SystemGenerator(
            forcefields=[protein_force_field, water_force_field],
            forcefield_kwargs=forcefield_kwargs
        )
        system = system_generator.create_system(modeller.topology)
    else:
        system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
            forcefield_kwargs=forcefield_kwargs,
            nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
        )
        system = system_generator.create_system(modeller.topology)
    
    return modeller, system
```

## State Restoration Process

### 1. Configuration Loading
```python
# Load original simulation parameters
with open(f'{restart_dir}/restart_setup.yml', 'r') as f:
    setup = yaml.safe_load(f)

# Extract key parameters
temperature = setup['temperature']
step_size = setup['step_size']
friction_coeff = setup['friction_coeff']
reporting_interval = setup['reporting_interval']
```

### 2. System Reconstruction
The restart process recreates the exact same system as the original simulation:
- **Force Field Consistency**: Uses identical force field parameters
- **Topology Preservation**: Maintains exact atom ordering and connectivity
- **Parameter Matching**: Applies same simulation parameters

### 3. State Loading
```python
# Load complete simulation state
simulation.loadState(f'{restart_dir}/last_state.xml')

# State includes:
# - Particle positions and velocities
# - Box vectors and periodic boundary conditions
# - Integrator state and step count
# - Random number generator state
```

## Integration with RestartRunner

### RestartRunner Workflow
The Restart component integrates with the RestartRunner to provide complete restart functionality:

```python
class Restarter:  # In runners/restartRunner.py
    """
    Molecular dynamics simulation runner for restarting simulations from saved states.
    """
    
    def run(self):
        """Execute the restart molecular dynamics simulation."""
        # Load original parameters
        setup = self.setup
        
        # Create integrator with original settings
        integrator = LangevinIntegrator(
            self.temperature, 
            self.friction_coeff, 
            self.step_size
        )
        
        # Add barostat for solvated systems
        if setup['solvate']:
            self.system.addForce(openmm.MonteCarloBarostat(
                1 * unit.atmospheres, self.temperature, 25))
        
        # Create simulation and load state
        simulation = Simulation(self.modeller.topology, self.system, integrator)
        simulation.loadState(f'{self.config.restart}/last_state.xml')
        
        # Continue simulation with incremented trajectory numbering
        last_state = setup['last_state'] + 1
        output_traj_dcd = f'output_traj_{last_state}.dcd'
        
        # Run simulation (time-based or step-based)
        if self.config.clock is not None:
            simulation.runForClockTime(self.config.clock * unit.minute)
        else:
            simulation.step(self.config.steps)
        
        # Save new state and update configuration
        simulation.saveState(f'{self.config.restart}/last_state.xml')
        setup['last_state'] = last_state
        
        # Update restart configuration
        yaml.dump(setup, open(f'{self.config.restart}/restart_setup.yml', 'w'))
```

## Trajectory Management

### 1. Trajectory Numbering
The restart system automatically manages trajectory file numbering:
```python
# Original simulation: output_traj_0.dcd
# First restart: output_traj_1.dcd  
# Second restart: output_traj_2.dcd
# etc.

last_state = setup['last_state'] + 1
output_traj_dcd = f'output_traj_{last_state}.dcd'
```

### 2. Continuous Trajectory
- **Seamless Continuation**: No gaps in trajectory data
- **Time Consistency**: Maintains proper temporal ordering
- **Analysis Compatibility**: Combined trajectories work with analysis tools

### 3. State Management
```python
# Update state counter for next restart
setup['last_state'] = last_state

# Save updated configuration
yaml.dump(setup, open(f'{restart_dir}/restart_setup.yml', 'w'))
```

## Error Handling and Validation

### 1. File Validation
```python
def validate_restart_files(restart_dir):
    """Validate presence and integrity of restart files."""
    required_files = [
        'restart_setup.yml',
        'restart_model.pdb',
        'last_state.xml'
    ]
    
    for file in required_files:
        file_path = os.path.join(restart_dir, file)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Required restart file missing: {file}")
```

### 2. Configuration Validation
```python
def validate_restart_config(setup):
    """Validate restart configuration parameters."""
    required_params = [
        'protein_force_field', 'temperature', 'step_size',
        'friction_coeff', 'reporting_interval', 'solvate'
    ]
    
    for param in required_params:
        if param not in setup:
            raise KeyError(f"Required parameter missing from restart config: {param}")
```

### 3. State Consistency
- **Topology Matching**: Ensures topology consistency between runs
- **Parameter Validation**: Verifies simulation parameters match
- **Force Field Consistency**: Confirms identical force field setup

## Performance Considerations

### 1. Efficient Loading
- **Optimized File I/O**: Fast loading of large state files
- **Memory Management**: Efficient handling of system data
- **Minimal Overhead**: Quick restart with minimal setup time

### 2. State Preservation
- **Complete State**: Preserves all simulation state information
- **Numerical Precision**: Maintains full precision of coordinates and velocities
- **Reproducibility**: Ensures identical continuation of simulation

### 3. Scalability
- **Large Systems**: Handles restart of large protein complexes
- **Long Simulations**: Supports multi-segment long simulations
- **Multiple Restarts**: Efficient handling of multiple restart cycles

## Usage Examples

### Basic Restart
```bash
# Restart simulation from previous run
python -m EasyMD --restart out_0/

# Continue for additional time
python -m EasyMD --restart out_0/ --clock 120

# Continue for additional steps  
python -m EasyMD --restart out_0/ --steps 5000000
```

### Programmatic Usage
```python
from EasyMD.restart import Restarter

# Initialize restarter
setup = yaml.load(open('out_0/restart_setup.yml'))
restarter = Restarter(setup, 'out_0', forcefield_kwargs)

# Prepare system for restart
if os.path.exists('out_0/ligand.sdf'):
    modeller, system = restarter.prep_restart_ligand()
else:
    modeller, system = restarter.prep_restart()
```

## Benefits

1. **Seamless Continuation**: Resume simulations exactly where they left off
2. **Long Simulation Support**: Enable microsecond-scale simulations through segmentation
3. **Resource Management**: Optimize use of time-limited computational resources
4. **Fault Tolerance**: Recover from system failures or interruptions
5. **Flexible Scheduling**: Adapt to different computational environments
6. **Complete State Preservation**: Maintain full simulation state and history
7. **Analysis Compatibility**: Generate continuous trajectories for analysis

## Integration Points

### With ArgManager
- **Restart Detection**: Identifies restart mode from command-line arguments
- **Parameter Validation**: Validates restart directory and required files
- **Configuration Loading**: Loads original simulation parameters

### With SysGenerator
- **System Reconstruction**: Recreates identical molecular systems
- **Force Field Consistency**: Maintains same force field parameters
- **Topology Preservation**: Ensures exact topology matching

### With SimRunner
- **Runner Selection**: Automatically selects RestartRunner for restart mode
- **State Management**: Handles state loading and trajectory continuation
- **Output Management**: Manages incremented output file naming

The Restart component provides essential functionality for long-term molecular dynamics simulations, enabling researchers to conduct extended studies while efficiently managing computational resources and maintaining complete simulation continuity.