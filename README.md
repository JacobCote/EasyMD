# EasyMD

Welcome to EasyMD, a comprehensive Python package for molecular dynamics simulations using OpenMM. EasyMD simplifies the process of setting up and running molecular dynamics simulations, allowing you to explore the behavior of biomolecules in different environments with ease.

## Key Features

- **Easy-to-use command-line interface** with comprehensive help and validation
- **Flexible solvation options**: explicit solvent with periodic boundary conditions or implicit solvent (GBIS)
- **Protein-ligand complex simulations** with automatic ligand parameterization
- **Robust structure preparation** including missing residue handling and protonation state assignment
- **Multiple force field support** (Amber, CHARMM, OpenFF)
- **Restart capabilities** for long simulations
- **Comprehensive validation system** with helpful error messages and warnings
- **GPU acceleration** via CUDA when available

For detailed information on EasyMD simulations, see [More Info](#more-info-on-easymd-simulations)
## Installation
To install EasyMD, simply run the following commands:
```bash
git clone https://github.com/JacobCote/EasyMD.git
cd EasyMD
conda env create -f requirements.yml
conda activate MdEnv
# if needed, install openmm with a specific cuda version
conda install -c conda-forge openmm cudatoolkit=11.8
```

## Usage

EasyMD provides a comprehensive command-line interface with organized argument groups and detailed help. Run EasyMD as a Python module:

```bash
python -m EasyMD [options]
```

### Getting Help

For comprehensive help with all available options:
```bash
python -m EasyMD --help
```

For version information:
```bash
python -m EasyMD --version
```

### Quick Start Examples

**Basic protein simulation (explicit solvent):**
```bash
python -m EasyMD --protein protein.pdb --steps 5000000 --solvate
```

**Protein-ligand complex:**
```bash
python -m EasyMD --protein complex.pdb --ligand ATP --steps 10000000 --solvate --temperature 310
```

**Fast implicit solvent simulation:**
```bash
python -m EasyMD --protein protein.pdb --steps 1000000 --GBIS
```

**Time-based simulation:**
```bash
python -m EasyMD --protein protein.pdb --clock 60 --solvate
```

**Using configuration file:**
```bash
python -m EasyMD --config simulation.yml
```

**Restart simulation:**
```bash
python -m EasyMD --restart out_0/
```

### Argument Groups

EasyMD organizes arguments into logical groups for better usability:

**Input/Output:**
- `-p, --protein`: Path to protein PDB file (required for new simulations)
- `-l, --ligand`: Ligand residue name as it appears in PDB (e.g., LIG, MOL, ATP)
- `-o, --outdir`: Output directory (auto-generated if not specified)

**Simulation Parameters:**
- `-s, --steps`: Number of simulation steps (mutually exclusive with --clock)
- `-z, --step-size`: Integration step size in picoseconds (default: 0.002 ps)
- `-f, --friction-coeff`: Langevin friction coefficient in 1/ps (default: 1.0)
- `-i, --interval`: Reporting interval for trajectory and log output (default: 1000)
- `-t, --temperature`: Simulation temperature in Kelvin (default: 300 K)
- `-e, --equilibration-steps`: Number of equilibration steps before production (default: 200)

**Solvation (choose one):**
- `--solvate`: Use explicit solvent with periodic boundary conditions
- `--GBIS`: Use Generalized Born implicit solvent (faster, less accurate)

**Solvation Parameters:**
- `--padding`: Solvent box padding around protein in Angstroms (default: 10 Å)
- `--water-model`: Water model for explicit solvation (default: tip3p)
- `--positive-ion`: Positive ion type for neutralization (default: Na+)
- `--negative-ion`: Negative ion type for neutralization (default: Cl-)
- `--ionic-strength`: Target ionic strength in Molar (default: 0.1 M)
- `--no-neutralize`: Skip automatic system neutralization

**Force Fields:**
- `--protein-force-field`: Protein force field (default: amber14-all.xml)
- `--ligand-force-field`: Small molecule force field (default: openff-2.2.0)
- `--water-force-field`: Water force field (default: amber/tip3p_standard.xml)

**Structure Preparation:**
- `--remove`: Molecule names to remove from structure (default: ['DMS'])
- `--keep-water`: Preserve crystal water molecules from PDB (default: remove all water)
- `--ph`: pH for protonation state assignment (default: 7.0)

**Advanced Options:**
- `-r, --restart`: Restart simulation from specified directory containing state files
- `--clock`: Simulation time duration in minutes - alternative to --steps
- `--simulated-annealing`: Use simulated annealing protocol instead of standard MD

**Missing Residue Handling:**
- `--missing-residues`: Strategy for handling missing residues (auto, none, non-terminal, terminal-only, all)
- `--max-terminal-residues`: Maximum number of terminal residues to add per chain (default: 5)
- `--conservative-missing`: Use conservative approach for missing residues

### Configuration Files

EasyMD supports YAML configuration files for complex setups:

```yaml
# simulation.yml
protein: "protein.pdb"
ligand: "ATP"
steps: 10000000
temperature: 310
solvate: true
ionic_strength: 0.15
keep_water: true
missing_residues: "terminal-only"
```

Run with: `python -m EasyMD --config simulation.yml`

### Output Files

EasyMD generates the following output files:
- **Trajectory files**: `output_traj_*.dcd` - simulation trajectories
- **Log file**: `log.txt` - energies, temperature, and other properties
- **Final structure**: `last_state_*.pdb` - final system configuration
- **System topology**: `topology.pkl` - pickled OpenMM topology
- **Restart files**: `last_state.xml`, `restart_setup.yml` - for continuing simulations

### Validation and Error Handling

EasyMD includes a comprehensive validation system that:
- **Validates input parameters** before starting simulations
- **Provides clear error messages** with suggested fixes
- **Shows warnings** for potentially problematic settings
- **Displays configuration summary** before starting

### Force Fields and Compatibility

- **Protein force fields**: Amber (amber14-all.xml, amber99sb-ildn.xml), CHARMM (charmm36.xml)
- **Ligand force fields**: OpenFF (openff-2.2.0, recommended), GAFF
- **Water models**: TIP3P (default), SPC/E, TIP4P-Ew, TIP5P
- **DNA systems**: AlphaFold3 structures work best with Amber14 force field

For available force fields, see the [OpenMM Force Fields documentation](https://ommprotocol.readthedocs.io/en/latest/forcefields.html).

### Cluster Usage

Example scripts for cluster usage are provided in the `examples/` folder. These demonstrate how to set up simulations on HPC systems with proper resource allocation and job management.

## Testing

EasyMD includes a comprehensive test suite to ensure reliability:

```bash
# Run all tests
python -m pytest src/EasyMD/tests/

# Run specific test categories
python -m pytest src/EasyMD/tests/test_arg_manager.py  # Argument parsing tests
python -m pytest src/EasyMD/tests/test_validation.py   # Validation system tests
python -m pytest src/EasyMD/tests/test_integration.py  # Integration tests
```

Test data files are organized in `src/EasyMD/tests/data/` for clean separation from the main codebase.

## Troubleshooting

### Common Issues

**Missing protein file error:**
```bash
❌ Protein PDB file is required. Use --protein <file.pdb>
```
**Solution:** Provide a valid PDB file path: `--protein your_protein.pdb`

**Solvation method not specified:**
```bash
❌ Must choose exactly one solvation method: either --solvate OR --GBIS
```
**Solution:** Choose either explicit (`--solvate`) or implicit (`--GBIS`) solvation

**Conflicting duration methods:**
```bash
❌ Cannot specify both --steps and --clock. Choose one simulation duration method
```
**Solution:** Use either `--steps 10000` OR `--clock 60` (not both)

### Getting Detailed Help

For parameter-specific help and examples:
```bash
python -m EasyMD --help | grep -A5 "missing-residues"  # Help for specific option
```

### Validation Messages

EasyMD provides detailed validation with:
- ✅ **Success messages** with configuration summary
- ⚠️ **Warnings** for potentially problematic settings (simulation continues)
- ❌ **Errors** with specific fixes required (simulation stops)

### Performance Tips

- Use `--GBIS` for faster simulations (less accurate)
- Use `--solvate` for production simulations (more accurate, slower)
- Adjust `--interval` to control output frequency
- Use `--clock` for time-limited cluster jobs

## More info on EasyMD simulations
EasyMD will automatically build a system and run the simulation based on the options provided. Protonation states of amino acids will be determined from the ph of the system. Disulfide bonds will be predicted based on proximity and orientation of Cys residus.

Solvated simulations are made with a PME periodic system using an NPT ensemble. 

Implicit solvent simulations use a generalized Born based implicit solvent.

Ligand charges are calculated using openbabel's python api. Forces for ligands are calculated using the openff-2.2.0 forcefield.

EasyMD uses the CUDA toolkit for GPU acceleration if available. 
## Roadmap

- [ ] Add an analysis tool documentation
- [ ] Add simulated annealing documentation
- [ ] Add possibility to install via pip
- [ ] add support for multiple ligands 
- [ ] Molecular dynamics based ligand docking

See the [open issues](https://github.com/JacobCote/EasyMD/issues) for a full list of proposed features (and known issues).

## Contributing
If you would like to contribute to EasyMD, please fork the repository and submit a pull request. We welcome any contributions, including bug fixes, new features, and documentation improvements.

## License
EasyMD is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Contact
If you have any questions or suggestions, feel free to open an issue or contact me at jacobcote@ulaval.ca.



