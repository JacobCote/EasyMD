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

**Analyze trajectory:**
```bash
python -m EasyMD analyze out_0/ --rmsd --rmsf --save-data
```

**Get PDB structure information:**
```bash
python -m EasyMD info protein.pdb
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

## Trajectory Analysis

EasyMD includes a comprehensive analysis module for post-simulation trajectory analysis. The analysis tools provide insights into structural dynamics, stability, and conformational changes.

### Analysis Usage

Run trajectory analysis using the `analyze` subcommand:

```bash
python -m EasyMD analyze [trajectory_directory] [options]
```

### Analysis Help

For comprehensive analysis help:
```bash
python -m EasyMD analyze --help
```

### Analysis Examples

**Basic RMSD and RMSF analysis:**
```bash
python -m EasyMD analyze out_0/ --rmsd --rmsf
```

**Perform all available analyses:**
```bash
python -m EasyMD analyze out_0/ --all --save-data
```

**Custom analysis with specific parameters:**
```bash
python -m EasyMD analyze out_0/ --rmsd --atom-selection backbone \
                 --reference-frame 10 --skip-frames 5
```

**High-quality publication plots:**
```bash
python -m EasyMD analyze out_0/ --rmsf --output-format pdf --dpi 600 \
                 --plot-style seaborn
```

**Analyze specific frame range:**
```bash
python -m EasyMD analyze out_0/ --rmsd --start-frame 1000 --end-frame 5000
```

### Available Analyses

**Structural Metrics:**
- `--rmsd`: Root Mean Square Deviation - measures structural similarity to reference
- `--rmsf`: Root Mean Square Fluctuation - identifies flexible regions
- `--radius-gyration`: Radius of gyration - measures protein compactness

**Inter-molecular Analysis:**
- `--distances`: Center-of-mass distances between protein chains
- `--secondary-structure`: Secondary structure evolution over time

**Comprehensive Analysis:**
- `--all`: Perform all available analyses

### Analysis Parameters

**Atom Selection:**
- `--atom-selection`: Choose atoms for analysis (`all`, `backbone`, `ca`, `heavy`)
- `--reference-frame`: Reference frame for RMSD calculation (default: 0)

**Frame Selection:**
- `--start-frame`: Starting frame for analysis (default: 0)
- `--end-frame`: Ending frame for analysis (default: all frames)
- `--skip-frames`: Skip every N frames (default: 1, no skipping)

**Output Options:**
- `--output-dir`: Custom output directory (default: trajectory_dir/analysis)
- `--output-format`: Plot format (`png`, `pdf`, `svg`)
- `--dpi`: Plot resolution (default: 300)
- `--plot-style`: Matplotlib style (`default`, `seaborn`, `ggplot`)
- `--save-data`: Export analysis data as CSV files
- `--no-plots`: Skip plot generation (data only)

### Analysis Output Files

**Generated Plots:**
- `rmsd.png` - RMSD vs time with statistics
- `rmsf.png` - RMSF per residue
- `distances.png` - Inter-chain distance evolution
- `radius_gyration.png` - Radius of gyration vs time
- `secondary_structure.png` - Secondary structure heatmap and evolution

**Data Files (with --save-data):**
- `rmsd_data.csv` - RMSD values and frame numbers
- `rmsf_data.csv` - RMSF values per residue
- `distances_data.csv` - Inter-chain distances
- `radius_gyration_data.csv` - Radius of gyration values
- `secondary_structure_data.csv` - Secondary structure counts

### Analysis Validation

The analysis module includes comprehensive validation:
- **Input validation**: Checks for required trajectory files and topology
- **Parameter validation**: Ensures analysis parameters are reasonable
- **Error handling**: Clear error messages with suggested fixes
- **Progress reporting**: Real-time analysis progress and statistics

## PDB Structure Information

EasyMD includes a comprehensive PDB structure analysis tool that provides detailed information about protein structures before running simulations. This helps users understand their structures and identify potential issues.

### Info Usage

Run PDB structure analysis using the `info` subcommand:

```bash
python -m EasyMD info [pdb_file] [options]
```

### Info Help

For comprehensive info help:
```bash
python -m EasyMD info --help
```

### Info Examples

**Basic structure analysis:**
```bash
python -m EasyMD info protein.pdb
```

**Custom output file:**
```bash
python -m EasyMD info protein.pdb --output protein_analysis.info
```

**Summary format only:**
```bash
python -m EasyMD info protein.pdb --format summary
```

**JSON output for programmatic use:**
```bash
python -m EasyMD info protein.pdb --format json
```

**Terminal output only (no file):**
```bash
python -m EasyMD info protein.pdb --no-file
```

**Custom disulfide bond detection:**
```bash
python -m EasyMD info protein.pdb --disulfide-distance 3.0
```

### Analysis Features

**Chain Analysis:**
- Number of chains and their composition
- Residue sequences and ranges
- Protein vs non-protein content

**Missing Residues:**
- Automatic detection of gaps in residue numbering
- Gap size and location analysis
- Context information (residues before/after gaps)

**Ligands and Small Molecules:**
- Identification of bound ligands
- Crystallization agents and buffer components
- Molecular composition and element analysis

**Water Molecules:**
- Water molecule count and distribution
- Chain-specific water analysis
- Spatial distribution statistics

**Disulfide Bonds:**
- Cysteine residue identification
- Distance-based disulfide bond prediction
- Inter-chain and intra-chain bond analysis

**Metal Ions:**
- Detection of metal ions (Mg²⁺, Zn²⁺, Ca²⁺, etc.)
- Spatial distribution and occupancy
- Element-specific statistics

**Modified Residues:**
- Non-standard amino acid detection
- Post-translational modifications
- MODRES record analysis

### Analysis Parameters

**Detection Thresholds:**
- `--disulfide-distance`: Maximum distance for disulfide bonds (default: 2.5 Å)
- `--water-threshold`: Minimum waters to show details (default: 10)
- `--missing-threshold`: Minimum gap size to report (default: 3)

**Output Options:**
- `--output`: Custom output file name
- `--format`: Output format (`detailed`, `summary`, `json`)
- `--no-file`: Terminal output only
- `--no-color`: Disable colored terminal output
- `--quiet`: Minimal output mode

### Output Files

**Generated Reports:**
- `<pdb_name>.info` - Detailed structure analysis report (default)
- Custom report file specified with `--output`
- Terminal display with colored formatting

**Report Contents:**
- Header information (PDB ID, classification, date)
- Chain composition and sequences
- Missing residue gaps and locations
- Ligand identification and properties
- Water molecule distribution
- Disulfide bond predictions
- Metal ion locations
- Modified residue catalog

### Info Validation

The info module includes comprehensive validation:
- **File validation**: Checks PDB file existence and format
- **Parameter validation**: Ensures analysis parameters are reasonable
- **Error handling**: Clear error messages for file issues
- **Progress reporting**: Real-time analysis progress

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



