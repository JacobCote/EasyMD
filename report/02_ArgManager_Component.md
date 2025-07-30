# ArgManager Component

## Overview
The ArgManager component is responsible for command-line argument parsing, configuration file handling, and comprehensive input validation for EasyMD molecular dynamics simulations. It provides a robust interface between user input and the simulation engine.

## Architecture

### Class Structure
```python
class ArgManager:
    """
    Manages command-line arguments and configuration files for EasyMD molecular dynamics simulations.
    
    Attributes:
        parser (argparse.ArgumentParser): The argument parser instance
        args (argparse.Namespace): Parsed arguments after validation
    """
```

### Key Features
1. **Dual Input Support**: Command-line arguments and YAML configuration files
2. **Comprehensive Validation**: Extensive input validation with clear error messages
3. **Priority System**: CLI arguments override configuration file values
4. **User-Friendly Error Reporting**: Detailed error messages with quick fix suggestions

## Configuration Management

### YAML Configuration Support
The ArgManager supports YAML configuration files that can specify any command-line parameter:

```yaml
# Example configuration file
protein: "protein.pdb"
ligand: "ATP"
steps: 5000000
temperature: 310
solvate: true
padding: 12.0
ionic_strength: 0.15
```

### Configuration Processing Flow
1. **Early Config Detection**: Preliminary parse to detect `--config` argument
2. **YAML Loading**: Load and parse configuration file
3. **Argument Injection**: Convert YAML to CLI-style arguments
4. **Priority Handling**: Inject config args before user CLI args (allows CLI override)

```python
def _parse_config_file_if_provided(self):
    """Parse and inject configuration file arguments into sys.argv if --config is provided."""
    temp_args, _ = self.parser.parse_known_args()
    if temp_args.config:
        with open(temp_args.config, "r") as f:
            config_args = yaml.safe_load(f)
        
        # Convert YAML keys to CLI-style args and inject
        sys.argv = sys.argv[:1] + cli_args + sys.argv[1:]
```

## Argument Categories

### 1. Input/Output Parameters
- **Protein file**: Required PDB input (`--protein`)
- **Ligand specification**: Residue name in PDB (`--ligand`)
- **Output directory**: Auto-generated or custom (`--outdir`)

### 2. Simulation Parameters
- **Duration**: Steps (`--steps`) or time-based (`--clock`)
- **Integration**: Step size (`--step-size`), friction coefficient (`--friction-coeff`)
- **Temperature**: Simulation temperature (`--temperature`)
- **Reporting**: Output interval (`--interval`)
- **Equilibration**: Pre-production steps (`--equilibration-steps`)

### 3. Solvation Options (Mutually Exclusive)
- **Explicit Solvation**: `--solvate` with water models and ions
- **Implicit Solvation**: `--GBIS` for faster simulations

### 4. Solvation Parameters
- **Box setup**: Padding distance (`--padding`)
- **Water models**: TIP3P, SPC/E, TIP4P-Ew, etc. (`--water-model`)
- **Ions**: Positive/negative ion types (`--positive-ion`, `--negative-ion`)
- **Ionic strength**: Target concentration (`--ionic-strength`)

### 5. Force Field Selection
- **Protein**: AMBER, CHARMM force fields (`--protein-force-field`)
- **Ligand**: OpenFF, GAFF force fields (`--ligand-force-field`)
- **Water**: Corresponding water parameters (`--water-force-field`)

### 6. Structure Preparation
- **Molecule removal**: Artifacts to remove (`--remove`)
- **Water handling**: Keep crystal waters (`--keep-water`)
- **Protonation**: pH for state assignment (`--ph`)

### 7. Missing Residue Handling
- **Strategy**: Auto, none, terminal-only, etc. (`--missing-residues`)
- **Limits**: Maximum terminal residues (`--max-terminal-residues`)
- **Types**: Allowed terminal residue types (`--terminal-residue-types`)
- **Filtering**: Skip loops, conservative approach

### 8. Advanced Options
- **Restart**: Continue from previous simulation (`--restart`)
- **Simulated Annealing**: Alternative protocol (`--simulated-annealing`)

## Validation System

### Comprehensive Input Validation
The ArgManager performs extensive validation across multiple categories:

#### 1. Input File Validation
```python
def _validate_input_files(self, args, errors):
    """Validate input file arguments."""
    if not args.restart and not args.simulated_annealing:
        if not args.protein:
            errors.append("Protein PDB file is required. Use --protein <file.pdb>")
        elif not os.path.isfile(args.protein):
            errors.append(f"Protein file '{args.protein}' does not exist")
```

#### 2. Simulation Parameter Validation
- **Duration conflicts**: Prevents both `--steps` and `--clock`
- **Range checking**: Validates positive values and reasonable ranges
- **Physical constraints**: Temperature, step size, friction coefficient limits

#### 3. Solvation Validation
- **Method selection**: Ensures exactly one solvation method
- **Parameter ranges**: Padding, ionic strength, pH validation
- **Ion compatibility**: Validates supported ion types

#### 4. Restart Validation
- **Directory existence**: Checks restart directory
- **Required files**: Validates presence of restart files
  - `restart_setup.yml`
  - `restart_model.pdb`
  - `last_state.xml`

### Error Handling and User Experience

#### Error Classification
The validation system distinguishes between:
- **Errors**: Critical issues that prevent simulation
- **Warnings**: Potential issues that allow simulation to continue

#### User-Friendly Error Display
```python
def _display_validation_errors(self, actual_errors, warnings=None):
    """Display validation errors in a user-friendly format."""
    print("\n" + "="*60)
    print("❌ INPUT VALIDATION FAILED")
    print("="*60)
    
    for i, error in enumerate(actual_errors, 1):
        print(f"❌ {i}. {error}")
    
    self._provide_quick_fixes(actual_errors)
```

#### Quick Fix Suggestions
The system provides specific suggestions for common errors:
- Missing protein file → Add `--protein your_protein.pdb`
- Duration conflicts → Use either `--steps` OR `--clock`
- Solvation conflicts → Choose `--solvate` OR `--GBIS`

### Validation Success Display
For successful validation, the system provides:
- Confirmation of validation success
- Summary of warnings (if any)
- Simulation parameter overview

## Integration Points

### With Main Entry Point
```python
# In __main__.py
argManager = ArgManager(parser)
config = argManager.get_args()
```

### With Other Components
The validated configuration is passed to:
- **SysGenerator**: For system preparation
- **SimRunner**: For simulation execution
- **Utils**: For platform detection and setup

## Error Prevention Strategy

### 1. Early Validation
All validation occurs before any simulation setup, preventing wasted time on invalid configurations.

### 2. Clear Error Messages
Each error message includes:
- What went wrong
- Expected values or formats
- Specific fix suggestions

### 3. Comprehensive Coverage
Validation covers:
- File existence and accessibility
- Parameter ranges and physical constraints
- Logical consistency between options
- Force field compatibility

### 4. Warning System
Non-critical issues generate warnings that:
- Alert users to potential problems
- Allow simulation to continue
- Provide optimization suggestions

## Usage Examples

### Basic Usage
```python
import argparse
from EasyMD.argManager import ArgManager

parser = argparse.ArgumentParser()
manager = ArgManager(parser)
config = manager.get_args()
```

### Configuration File Usage
```bash
# Create config.yml
python -m EasyMD --config simulation_config.yml

# Override config with CLI
python -m EasyMD --config config.yml --temperature 310 --steps 10000000
```

## Benefits

1. **Robust Input Handling**: Prevents invalid simulations before they start
2. **User-Friendly Interface**: Clear error messages and suggestions
3. **Flexible Configuration**: Support for both CLI and file-based configuration
4. **Comprehensive Validation**: Covers all aspects of simulation setup
5. **Maintainable Code**: Modular validation functions for easy extension

The ArgManager component ensures that users can confidently set up simulations with proper validation and clear feedback, significantly reducing setup errors and improving the overall user experience.