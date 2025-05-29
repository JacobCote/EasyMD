
class ArgManager():
    def __init__(self,parser):
        self.parser = parser
        self.initialize()
        self._argsSanityCheck()

    def initialize(self,):
        self.parser.add_argument("-p", "--protein", required=False, help="Protein PDB file")
        self.parser.add_argument("-l", "--ligand", required=False, help="Ligand name in pdb file (often LIG, check your pdb file to be sure of the name)")
        self.parser.add_argument("-o", "--output", default=None, help="Name of an output directory")
        self.parser.add_argument("-s", "--steps", type=int, default=None, help="Number of steps")
        self.parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps")
        self.parser.add_argument("-f", "--friction-coeff", type=float, default=1, help="Friction coefficient (ps)")
        self.parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
        self.parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
        self.parser.add_argument("--solvate", action='store_true', help="Add solvent box")
        self.parser.add_argument("--GBIS", action='store_true', help="Doesn't add solvent box, use Born generalize implicit solvent")
        self.parser.add_argument("--padding", type=float, default=10, help="Padding for solvent box (A)")
        self.parser.add_argument("--water-model", default="tip3p",
                            choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"],
                            help="Water model for solvation")
        self.parser.add_argument("--positive-ion", default="Na+", help="Positive ion for solvation")
        self.parser.add_argument("--negative-ion", default="Cl-", help="Negative ion for solvation")
        self.parser.add_argument("--ionic-strength", type=float, default=0.1, help="Ionic strength for solvation")
        self.parser.add_argument("--no-neutralize", action='store_true', help="Don't add ions to neutralize")
        self.parser.add_argument("-e", "--equilibration-steps", type=int, default=200, help="Number of equilibration steps")
        self.parser.add_argument("--protein-force-field", default='amber14-all.xml', help="Protein force field")
        self.parser.add_argument("--ligand-force-field", default='openff-2.2.0', help="Ligand force field")
        self.parser.add_argument("--water-force-field", default='amber/tip3p_standard.xml', help="Water force field")
        self.parser.add_argument('--remove', nargs='*', help='Space separated molecules name to remove ex: --remove DMS LIG CA MG ... ', required=False, default=['DMS'])
        self.parser.add_argument('--ph', type=float, help='Ph for the protonation state of the residus', required=False, default=7.0)
        self.parser.add_argument("--restart", action='store_true', help="Use the program in restart mode.",default=False)
        self.parser.add_argument("--restart_dir",type=str, help="path to the restart files", required=False, default='None')
        self.arser.add_argument('--clock', type=float, help='Run the simulation based on clock time in minutes instead of steps.', required=False, default=None)
        self.parser.add_argument('--simulated-annealing', action='store_true', help='Run a simulated annealing simulation', required=False, default=False)


    def getargs(self,):
        return self.parser.parse_args()
    def _argsSanityCheck(self,):
        args = self.getargs()
        if not args.simulated_annealing :

            if args.clock is not None and args.steps is not None:
                print('Please choose either --steps or --clock')
                exit(1)
            if args.steps is None and args.clock is None:
                print('Please provide either --steps or --clock')
                exit(1)

            if (args.solvate + args.GBIS % 2 == 0) and not args.restart :
                print('Please choose either --solvate or --GBIS')
                exit(1)
                
            if args.restart and  args.restart_dir is None:
                print('Please  directory with the restart files with --restart_dir <RESTARTDIR>')
                exit(1)
            