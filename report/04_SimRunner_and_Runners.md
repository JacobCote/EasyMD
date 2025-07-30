# SimRunner and Runners Components

## Overview
The SimRunner component serves as the central coordinator for molecular dynamics simulations, while the Runners package contains specialized simulation protocols. Together, they provide a flexible, modular system for executing different types of MD simulations with appropriate parameters and conditions.

## SimRunner Architecture

### Class Structure
```python
class SimRunner:
    """
    Central coordinator for molecular dynamics simulations.
    
    Attributes:
        config: Configuration object containing simulation parameters
        modeller: OpenMM Modeller object with system topology and positions
        system: OpenMM System object with forces and parameters
        setup: Setup configuration for restart simulations
        sim: Selected runner instance based on simulation type
    """
```

### Simulation Type Detection
The SimRunner automatically selects the appropriate runner based on configuration:

```python
def _setupSim(self):
    """Select appropriate simulation runner based on configuration."""
    if self.config.restart != None:
        return Restarter(self.setup, self.config, modeller=self.modeller, system=self.system)
    if self.config.simulated_annealing:
        return AnnealingRunner(self.config, modeller=self.modeller, system=self.system)
    elif self.config.GBIS:
        return GBISRunner(self.config, modeller=self.modeller, system=self.system)
    else:
        return SolvatedRunner(self.config, modeller=self.modeller, system=self.system)
```

### Runner Selection Logic
1. **Restart**: Continue from previous simulation state
2. **Simulated Annealing**: Temperature-based conformational sampling
3. **GBIS**: Implicit solvent simulations
4. **Solvated** (default): Explicit solvent simulations

## Runner Components

### 1. SolvatedRunner - Explicit Solvent Simulations

#### Purpose
Handles MD simulations with explicit water molecules and periodic boundary conditions, providing the most accurate representation of biological systems.

#### Key Features
- **NPT Ensemble**: Constant temperature and pressure control
- **Monte Carlo Barostat**: Pressure regulation at 1 atmosphere
- **Periodic Boundary Conditions**: Proper treatment of long-range interactions
- **Comprehensive Equilibration**: Multi-stage system preparation

#### Simulation Workflow
```python
def run(self):
    """Execute complete solvated MD simulation workflow."""
    # 1. Setup integrator and barostat
    self.system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, self.temperature, 25))
    
    # 2. Energy minimization (10,000 iterations)
    simulation.minimizeEnergy(maxIterations=10000)
    
    # 3. Equilibration at target temperature
    simulation.context.setVelocitiesToTemperature(self.temperature)
    simulation.step(self.config.equilibration_steps)
    
    # 4. Production simulation
    if self.config.clock is not None:
        simulation.runForClockTime(self.config.clock * unit.minute)
    else:
        simulation.step(self.config.steps)
```

#### Error Handling
Advanced error detection for common simulation failures:
```python
try:
    simulation.minimizeEnergy(maxIterations=10000)
except Exception as e:
    if "NaN" in str(e) or "coordinate is NaN" in str(e):
        print("❌ Energy minimization failed due to NaN coordinates")
        print("Troubleshooting suggestions:")
        print("  1. Check input PDB structure for overlapping atoms")
        print("  2. Verify ligand coordinates are reasonable")
        print("  3. Try using --keep-water flag")
        raise RuntimeError(f"Energy minimization failed: {e}")
```

#### Output Files
- `output_traj_0.dcd`: Trajectory data
- `log.txt`: Energy, temperature, and step information
- `minimised.pdb`: Structure after energy minimization
- `last_state.xml`: Complete simulation state for restart
- `last_state_0.pdb`: Final structure
- `topology.pkl`: Pickled OpenMM topology
- `restart_setup.yml`: Configuration for restart

### 2. GBISRunner - Implicit Solvent Simulations

#### Purpose
Provides faster MD simulations using the Generalized Born Implicit Solvent model, suitable for conformational sampling and preliminary studies.

#### Key Features
- **No Periodic Boundaries**: Simplified system without explicit solvent box
- **No Barostat**: Pressure control not needed for implicit solvent
- **Faster Performance**: Reduced computational overhead
- **Reduced Minimization**: Fewer iterations needed due to simpler system

#### Advantages
- **Speed**: 5-10x faster than explicit solvent
- **Memory Efficiency**: Smaller system size
- **Conformational Sampling**: Good for studying protein flexibility
- **Preliminary Studies**: Rapid assessment of system behavior

#### Limitations
- **Reduced Accuracy**: Less precise than explicit solvent
- **Missing Specific Interactions**: No explicit water-mediated effects
- **Limited Applicability**: Not suitable for all biological systems

### 3. AnnealingRunner - Simulated Annealing Protocol

#### Purpose
Implements simulated annealing for enhanced conformational sampling and optimization, helping systems escape local energy minima.

#### Annealing Protocol
```python
def run(self):
    """Execute simulated annealing protocol."""
    # Standard setup and minimization
    simulation.minimizeEnergy(maxIterations=10000)
    
    # Annealing cycles: 2000 cycles of 100 steps each
    # Temperature increases by 0.1K per cycle
    # Total temperature increase: 200K over protocol
```

#### Applications
- **Protein Folding**: Help proteins find native conformations
- **Conformational Search**: Explore multiple stable states
- **Optimization Problems**: Find global energy minima
- **Structure Refinement**: Improve initial structures

#### Protocol Details
- **Cycles**: 2000 annealing cycles
- **Steps per Cycle**: 100 MD steps
- **Temperature Increment**: 0.1K per cycle
- **Total Temperature Range**: Base temperature + 200K

### 4. Restarter - Simulation Continuation

#### Purpose
Enables continuation of long simulations from previously saved states, supporting segmented execution and recovery from interruptions.

#### Restart Workflow
```python
def run(self):
    """Continue simulation from saved state."""
    # 1. Load original simulation parameters
    # 2. Create integrator with same settings
    # 3. Load last saved state from XML
    # 4. Continue simulation with incremented trajectory numbering
    # 5. Update restart configuration
```

#### Required Files
- `restart_setup.yml`: Original simulation parameters
- `last_state.xml`: Complete OpenMM state
- `restart_model.pdb`: System topology and positions
- `ligand.sdf`: Ligand structure (if present)

#### State Management
- **Trajectory Numbering**: Automatic increment of output files
- **State Continuity**: Exact continuation from previous point
- **Parameter Consistency**: Uses original simulation settings
- **Progress Tracking**: Maintains cumulative simulation time

## Common Features Across Runners

### 1. Langevin Integration
All runners use Langevin dynamics for temperature control:
```python
integrator = LangevinIntegrator(
    temperature,      # Target temperature
    friction_coeff,   # Friction coefficient (1/ps)
    step_size        # Integration time step (ps)
)
```

### 2. Reporting System
Comprehensive output reporting:
- **Console Output**: Real-time progress with energy and temperature
- **Trajectory Files**: DCD format for analysis
- **Log Files**: Detailed simulation data
- **State Files**: Complete system state for restart

### 3. Error Handling
Robust error detection and user guidance:
- **NaN Detection**: Identify coordinate failures
- **Minimization Failures**: Provide troubleshooting guidance
- **File I/O Errors**: Handle missing or corrupted files
- **Parameter Validation**: Check for reasonable values

### 4. Performance Monitoring
Built-in performance tracking:
- **Timing**: Wall-clock time measurement
- **Progress**: Step count and simulation time
- **Efficiency**: Steps per second calculation
- **Resource Usage**: Memory and CPU monitoring

## Integration Architecture

### With SysGenerator
- Receives prepared `modeller` and `system` objects
- Uses validated force field parameters
- Maintains system integrity throughout simulation

### With ArgManager
- Uses validated configuration parameters
- Applies user-specified simulation conditions
- Respects output directory and file naming

### With Analysis Tools
- Generates compatible trajectory formats
- Provides structured output for post-processing
- Maintains metadata for analysis workflows

## Simulation Modes

### Time-Based Execution
```bash
python -m EasyMD --protein protein.pdb --clock 60 --solvate
```
- Runs for specified wall-clock time
- Useful for time-limited resources
- Automatic progress saving

### Step-Based Execution
```bash
python -m EasyMD --protein protein.pdb --steps 5000000 --solvate
```
- Runs for specified number of MD steps
- Precise control over simulation length
- Predictable computational requirements

### Restart Execution
```bash
python -m EasyMD --restart out_0/
```
- Continues from previous simulation
- Maintains all original parameters
- Seamless trajectory continuation

## Performance Considerations

### Optimization Strategies
1. **Efficient Integration**: Optimized time step selection
2. **Memory Management**: Proper cleanup of temporary objects
3. **I/O Optimization**: Buffered file operations
4. **Platform Selection**: Automatic GPU/CPU detection

### Scalability Features
- **Large Systems**: Handles proteins with >100k atoms
- **Long Simulations**: Supports microsecond-scale simulations
- **Multiple Trajectories**: Automatic file management
- **Resource Monitoring**: Memory and performance tracking

## Quality Assurance

### Validation Checks
- **Energy Conservation**: Monitor total energy drift
- **Temperature Control**: Verify thermostat performance
- **Pressure Stability**: Check barostat effectiveness
- **Structural Integrity**: Detect unfolding or artifacts

### Output Validation
- **File Integrity**: Verify complete trajectory files
- **State Consistency**: Validate restart file completeness
- **Metadata Accuracy**: Check simulation parameters
- **Format Compliance**: Ensure standard file formats

The SimRunner and Runners components provide a comprehensive, flexible framework for molecular dynamics simulations, supporting multiple protocols while maintaining ease of use and robust error handling. The modular design allows for easy extension and customization while ensuring reliable, high-quality simulations.