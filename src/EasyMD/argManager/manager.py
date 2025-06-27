import argparse
import yaml
import sys

class ArgManager():
    def __init__(self, parser):
        self.parser = parser
        self._add_config_argument_first()
        self._parse_config_file_if_provided()
        self._initialize()
        self.args = self.parser.parse_args()
        self._args_sanity_check()

    def _add_config_argument_first(self):
        # This allows early parsing of --config before full parser is built
        self.parser.add_argument("-c", "--config", type=str, help="Path to YAML config file", required=False)

    def _parse_config_file_if_provided(self):
        # Temporary parse to check if --config was provided
        temp_args, _ = self.parser.parse_known_args()
        if temp_args.config:
            with open(temp_args.config, "r") as f:
                config_args = yaml.safe_load(f)

            # Convert YAML keys to CLI-style args
            cli_args = []
            for k, v in config_args.items():
                key = f"--{k.replace('_', '-')}"
                if isinstance(v, bool):
                    if v:
                        cli_args.append(key)
                elif isinstance(v, list):
                    cli_args.append(key)
                    cli_args.extend(str(item) for item in v)
                else:
                    cli_args.append(key)
                    cli_args.append(str(v))

            # Inject YAML args into sys.argv *before* user CLI args (to allow CLI override)
            sys.argv = sys.argv[:1] + cli_args + sys.argv[1:]

    def _initialize(self):
        self.parser.add_argument("-p", "--protein", required=False, help="Protein PDB file")
        self.parser.add_argument("-l", "--ligand", required=False, help="Ligand name in pdb file (often LIG)")
        self.parser.add_argument("-o", "--output", default=None, help="Output directory name")
        self.parser.add_argument("-s", "--steps", type=int, default=None, help="Number of simulation steps")
        self.parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps)")
        self.parser.add_argument("-f", "--friction-coeff", type=float, default=1, help="Friction coefficient (ps)")
        self.parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
        self.parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
        self.parser.add_argument("--solvate", action='store_true', help="Add solvent box")
        self.parser.add_argument("--GBIS", action='store_true', help="Use Generalized Born implicit solvent")
        self.parser.add_argument("--padding", type=float, default=10, help="Padding for solvent box (Å)")
        self.parser.add_argument("--water-model", default="tip3p",
                            choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"],
                            help="Water model for solvation")
        self.parser.add_argument("--positive-ion", default="Na+", help="Positive ion for solvation")
        self.parser.add_argument("--negative-ion", default="Cl-", help="Negative ion for solvation")
        self.parser.add_argument("--ionic-strength", type=float, default=0.1, help="Ionic strength")
        self.parser.add_argument("--no-neutralize", action='store_true', help="Do not add neutralizing ions")
        self.parser.add_argument("-e", "--equilibration-steps", type=int, default=200, help="Equilibration steps")
        self.parser.add_argument("--protein-force-field", default='amber14-all.xml', help="Protein force field")
        self.parser.add_argument("--ligand-force-field", default='openff-2.2.0', help="Ligand force field")
        self.parser.add_argument("--water-force-field", default='amber/tip3p_standard.xml', help="Water force field")
        self.parser.add_argument("--remove", nargs='*', default=['DMS'], help="Molecules to remove (e.g. DMS LIG)")
        self.parser.add_argument("--ph", type=float, default=7.0, help="Protonation pH")
        self.parser.add_argument("-r","--restart", type=str, default=False, help="Use restart mode")
        self.parser.add_argument("--restart-dir", type=str, default=None, help="Path to restart files")
        self.parser.add_argument("--clock", type=float, default=None, help="Run simulation based on wall time (min)")
        self.parser.add_argument("--simulated-annealing", action='store_true', default=False, help="Run simulated annealing")

    def get_args(self):
        return self.args

    def _args_sanity_check(self):
        args = self.args

        if not args.simulated_annealing and not args.restart :
            if args.clock is not None and args.steps is not None:
                print("❌ Please choose either --steps or --clock, not both.")
                exit(1)
            if args.clock is None and args.steps is None:
                print("❌ Please provide either --steps or --clock.")
                exit(1)
            if not args.restart:
                if (args.solvate and args.GBIS) or (not args.solvate and not args.GBIS):
                    print("❌ Please choose either --solvate or --GBIS.")
                    exit(1)

        """ if args.restart and (args.restart_dir is None or args.restart_dir == 'None'):
            print("❌ Please provide the restart directory with --restart-dir <RESTARTDIR>")
            exit(1)"""
    



