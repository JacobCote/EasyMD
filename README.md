# EasyMD

Welcome to EasyMD, a Python package for simulating complex molecular systems using OpenMM. With EasyMD, you can easily perform molecular dynamics simulations and explore the behavior of biomolecules in different environments.

With EasyMD, you can easily create explicit or implicit solvent simulations. EasyMD can also simulate protein-ligand complexes and protein-DNA complexes from a pdb file.

For more info on EasyMD see [here](#more-info-on-easymd-simulations)
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
To use EasyMD, you can run the `simulate.py` script with the following options:
```bash
usage: simulate.py [-h] [-p PROTEIN] [-l LIGAND] [-o OUTPUT] [-s STEPS]
        [-z STEP_SIZE] [-f FRICTION_COEFF] [-i INTERVAL]
        [-t TEMPERATURE] [--solvate] [--GBIS]
        [--padding PADDING]
        [--water-model]
        [--positive-ion POSITIVE_ION]
        [--negative-ion NEGATIVE_ION]
        [--ionic-strength IONIC_STRENGTH] [--no-neutralize]
        [-e EQUILIBRATION_STEPS]
        [--protein-force-field PROTEIN_FORCE_FIELD]
        [--ligand-force-field LIGAND_FORCE_FIELD]
        [--water-force-field WATER_FORCE_FIELD]
        [--remove [REMOVE ...]] [--ph PH] [--restart]
        [--restart_dir RESTART_DIR] [--clock CLOCK]
        [--simulated-annealing]
```

You can customize the simulation by providing the necessary options. Here are some of the available options:

- `-p, --protein`: Specify the protein PDB file.
- `-l, --ligand`: Specify the ligand name in the PDB file.
- `-o, --output`: Specify the base name for output files.
- `-s, --steps`: Specify the number of simulation steps.
- `-z, --step-size`: Specify the step size in ps.
- `-f, --friction-coeff`: Specify the friction coefficient in ps.
- `-i, --interval`: Specify the reporting interval.
- `-t, --temperature`: Specify the temperature in K.
- `--solvate`: Add a solvent box to the system.
- `--GBIS`: Use Born generalize implicit solvent instead of adding a solvent box.
- `--padding`: Specify the padding for the solvent box in A.
- `--water-model`: Specify the water model for solvation.
- `--positive-ion`: Specify the positive ion for solvation.
- `--negative-ion`: Specify the negative ion for solvation.
- `--ionic-strength`: Specify the ionic strength for solvation.
- `--no-neutralize`: Do not add ions to neutralize the system.
- `-e, --equilibration-steps`: Specify the number of equilibration steps.
- `--protein-force-field`: Specify the protein force field.
- `--ligand-force-field`: Specify the ligand force field.
- `--water-force-field`: Specify the water force field.
- `--remove`: Specify molecules to remove from the system.
- `--ph`: Specify the pH for the protonation state of the residues.
- `--restart`: Use the restart files for simulation.
- `--restart_dir`: Specify the path to the restart files.
- `--simulated-annealing`: starts a simulated annealing simulation.
- `--clock`: Run the simulation based on clock time instead of steps, which can be a good option when running on clusters. You need to take into account time for system creation, minimization and equilibration since the clock argument only applies to the simulation.

By default, EasyMD uses the following settings:
- Protein: No default
- Ligand: `None`
- Output: `out_X`
- Steps: `10000`
- Step size: `2.0`
- Friction coefficient: `1.0`
- Reporting interval: `1000`
- Temperature: `300 K`
- Solvate: `False`
- GBIS: `False`
- Padding: `10 A`
- Water model: `tip3p`
- Positive ion: `Na+`
- Negative ion: `Cl-`
- Ionic strength: `0.1`
- Neutralize: `True`
- Equilibration steps: `200`
- Protein force field: `amber14-all.xml`
- Ligand force field: `openff-2.2.0` (take note that GAFF forcefield is not available at this time. See https://github.com/openmm/openmmforcefields/releases/tag/0.13.0  )
- Water force field: `tip3p`
- Remove: `['DMS']`
- pH: `7.0`
- Restart: `False`
- Restart directory: `None`
- Clock: `False`
- Simulated annealing: `False`

Feel free to explore and experiment with different options to suit your needs.

To explore available forcefields, check [here](https://ommprotocol.readthedocs.io/en/latest/forcefields.html).

For now, systems with DNA coming from AlphaFold3 only work with the amber14 forcefield.

There are examples provided in the `examples` folder to setup a simulation on clusters (you need to move the scripts to the root directory of EasyMD).

#### Here is a basic example for an explicit solvent simulation, which will run for 60 minutes:
```bash
python3 simulate.py --protein test.pdb --ligand GNP -o test_explicit --solvate --clock 60
```
At the end of the 60 minutes, a restart setup file will be created. It is possible to restart the simulation for 60 minutes from the last state by specifying the output directory of the initial simulation:
```bash
python3 simulate.py --restart  -restart_dir test_explicit --clock 60
```

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



