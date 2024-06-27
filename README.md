# EasyMD

usage: simulateComplex.py [-h] [-p PROTEIN] [-l LIGAND] [-o OUTPUT] [-s STEPS]
                          [-z STEP_SIZE] [-f FRICTION_COEFF] [-i INTERVAL]
                          [-t TEMPERATURE] [--solvate] [--GBIS]
                          [--padding PADDING]
                          [--water-model {tip3p,spce,tip4pew,tip5p,swm4ndp}]
                          [--positive-ion POSITIVE_ION]
                          [--negative-ion NEGATIVE_ION]
                          [--ionic-strength IONIC_STRENGTH] [--no-neutralize]
                          [-e EQUILIBRATION_STEPS]
                          [--protein-force-field PROTEIN_FORCE_FIELD]
                          [--ligand-force-field LIGAND_FORCE_FIELD]
                          [--water-force-field WATER_FORCE_FIELD] [--CUDA]
                          [--remove [REMOVE ...]] [--ph PH] [--restart]
                          [--restart_dir RESTART_DIR] [--clock CLOCK]

Simulate

options:
  -h, --help            show this help message and exit
  -p PROTEIN, --protein PROTEIN
                        Protein PDB file (default: None)
  -l LIGAND, --ligand LIGAND
                        Ligand name in pdb file (often LIG, check your pdb
                        file to be sure of the name) (default: None)
  -o OUTPUT, --output OUTPUT
                        Base name for output files (default: output)
  -s STEPS, --steps STEPS
                        Number of steps (default: None)
  -z STEP_SIZE, --step-size STEP_SIZE
                        Step size (ps (default: 0.002)
  -f FRICTION_COEFF, --friction-coeff FRICTION_COEFF
                        Friction coefficient (ps) (default: 1)
  -i INTERVAL, --interval INTERVAL
                        Reporting interval (default: 1000)
  -t TEMPERATURE, --temperature TEMPERATURE
                        Temperature (K) (default: 300)
  --solvate             Add solvent box (default: False)
  --GBIS                Don't add solvent box, use Born generalize implicit
                        solvent (default: False)
  --padding PADDING     Padding for solvent box (A) (default: 10)
  --water-model {tip3p,spce,tip4pew,tip5p,swm4ndp}
                        Water model for solvation (default: tip3p)
  --positive-ion POSITIVE_ION
                        Positive ion for solvation (default: Na+)
  --negative-ion NEGATIVE_ION
                        Negative ion for solvation (default: Cl-)
  --ionic-strength IONIC_STRENGTH
                        Ionic strength for solvation (default: 0)
  --no-neutralize       Don't add ions to neutralize (default: False)
  -e EQUILIBRATION_STEPS, --equilibration-steps EQUILIBRATION_STEPS
                        Number of equilibration steps (default: 200)
  --protein-force-field PROTEIN_FORCE_FIELD
                        Protein force field (default: amber/ff14SB.xml)
  --ligand-force-field LIGAND_FORCE_FIELD
                        Ligand force field (default: openff-2.2.0)
  --water-force-field WATER_FORCE_FIELD
                        Water force field (default: amber/tip3p_standard.xml)
  --CUDA                Use CUDA platform (default: False)
  --remove [REMOVE ...]
                        Space separated molecules name to remove ex: --remove
                        DMS LIG CA MG ... (default: ['DMS'])
  --ph PH               Ph for the protonation state of the residus (default:
                        7.0)
  --restart             Use CUDA platform (default: False)
  --restart_dir RESTART_DIR
                        path to the restart files (default: None)
  --clock CLOCK         Run the simulation based on clock time in minutes
                        instead of steps. (default: None)
