import argparse
import yaml
import sys
import os
from pathlib import Path

class ArgManager:
    """
    Manages command-line arguments and configuration files for EasyMD molecular dynamics simulations.
    
    This class handles argument parsing with support for YAML configuration files, allowing users to
    specify simulation parameters either via command-line arguments or configuration files. CLI arguments
    take precedence over configuration file values.
    
    Attributes:
        parser (argparse.ArgumentParser): The argument parser instance.
        args (argparse.Namespace): Parsed arguments after validation.
    
    Example:
        >>> import argparse
        >>> parser = argparse.ArgumentParser()
        >>> manager = ArgManager(parser)
        >>> config = manager.get_args()
        >>> print(config.protein)  # Access protein file path
    """
    def __init__(self, parser):
        """
        Initialize the ArgManager with an ArgumentParser instance.
        
        Args:
            parser (argparse.ArgumentParser): The argument parser to configure.
            
        Raises:
            SystemExit: If argument validation fails during initialization.
        """
        self.parser = parser
        self._add_config_argument_first()
        self._parse_config_file_if_provided()
        self._initialize()
        self.args = self.parser.parse_args()
        self._args_sanity_check()

    def _add_config_argument_first(self):
        """
        Add the --config argument to the parser before other arguments.
        
        This method adds the configuration file argument early to allow for
        preliminary parsing and configuration file loading before the full
        argument parser is constructed.
        """
        self.parser.add_argument("-c", "--config", type=str, help="Path to YAML config file", required=False)

    def _parse_config_file_if_provided(self):
        """
        Parse and inject configuration file arguments into sys.argv if --config is provided.
        
        This method performs a preliminary parse to check if a configuration file was specified.
        If found, it loads the YAML configuration and converts the key-value pairs into
        command-line style arguments, injecting them into sys.argv before user-provided
        CLI arguments to allow CLI arguments to override config file values.
        
        Raises:
            FileNotFoundError: If the specified config file doesn't exist.
            yaml.YAMLError: If the config file contains invalid YAML.
        """
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
        """
        Initialize all command-line arguments for molecular dynamics simulation.
        
        This method defines all available command-line arguments including:
        - Input/output parameters (protein, ligand, output directory)
        - Simulation parameters (steps, temperature, step size, etc.)
        - Solvation options (water models, ions, padding)
        - Force field specifications
        - Special modes (restart, simulated annealing)
        """
        self.parser.add_argument("-p", "--protein", required=False, help="Protein PDB file")
        self.parser.add_argument("-l", "--ligand", required=False, help="Ligand name in pdb file (often LIG)")
        self.parser.add_argument("-o", "--outdir", default=None, help="Output directory name")
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
        self.parser.add_argument("--keep-water", action='store_true', help="Keep water molecules from PDB (default: remove)")
        self.parser.add_argument("--ph", type=float, default=7.0, help="Protonation pH")
        self.parser.add_argument("-r","--restart", type=str, default=None, help="Use restart mode, must give the directory to restart from")
        self.parser.add_argument("--clock", type=int, default=None, help="Run simulation based on wall time (min)")
        self.parser.add_argument("--simulated-annealing", action='store_true', default=False, help="Run simulated annealing")
        
        # Missing residue handling options
        self.parser.add_argument("--missing-residues", default="auto", 
                            choices=["auto", "none", "non-terminal", "terminal-only", "all"],
                            help="Strategy for adding missing residues: auto (default), none (skip all), non-terminal (only internal), terminal-only (only ends), all (add everything)")
        self.parser.add_argument("--max-terminal-residues", type=int, default=5,
                            help="Maximum number of missing residues to add at terminal ends (default: 5)")
        self.parser.add_argument("--terminal-residue-types", nargs='*', default=['ACE', 'NME', 'NH2', 'COOH'],
                            help="Allowed terminal residue types to add (default: ACE NME NH2 COOH)")
        self.parser.add_argument("--skip-missing-loops", action='store_true', default=False,
                            help="Skip adding missing residues in loop regions (may cause gaps)")
        self.parser.add_argument("--conservative-missing", action='store_true', default=False,
                            help="Use conservative approach: only add missing residues with high confidence")

    def get_args(self):
        return self.args

    def _args_sanity_check(self):
        """
        Perform comprehensive validation of all input arguments.
        
        This method validates all input parameters and provides clear, actionable
        error messages to help users fix their command-line arguments before
        starting the simulation.
        
        Raises:
            SystemExit: If any validation check fails, with descriptive error message.
        """
        args = self.args
        errors = []
        
        # Validate input files
        self._validate_input_files(args, errors)
        
        # Validate simulation parameters
        self._validate_simulation_parameters(args, errors)
        
        # Validate solvation settings
        self._validate_solvation_settings(args, errors)
        
        # Validate restart settings
        self._validate_restart_settings(args, errors)
        
        # Validate output settings
        self._validate_output_settings(args, errors)
        
        # Validate physical parameters
        self._validate_physical_parameters(args, errors)
        
        # Validate force field parameters
        self._validate_force_fields(args, errors)
        
        # Validate ligand parameters
        self._validate_ligand_parameters(args, errors)
        
        # Validate missing residue parameters
        self._validate_missing_residue_parameters(args, errors)
        
        # Separate errors from warnings
        actual_errors = [e for e in errors if not e.startswith("Warning:")]
        warnings = [e for e in errors if e.startswith("Warning:")]
        
        # Display warnings but continue execution
        if warnings:
            self._display_warnings(warnings)
        
        # Only exit on actual errors
        if actual_errors:
            self._display_validation_errors(actual_errors, warnings)
            sys.exit(1)
        
        # Display validation success message (even if there were warnings)
        self._display_validation_success(args, warnings)
    
    def _validate_input_files(self, args, errors):
        """Validate input file arguments."""
        # Check protein file for non-restart modes
        if not args.restart and not args.simulated_annealing:
            if not args.protein:
                errors.append("Protein PDB file is required. Use --protein <file.pdb>")
            elif not os.path.isfile(args.protein):
                errors.append(f"Protein file '{args.protein}' does not exist or is not accessible")
            elif not args.protein.lower().endswith(('.pdb', '.pdb.gz')):
                errors.append(f"Protein file '{args.protein}' should be a PDB file (.pdb or .pdb.gz)")
        
        # Check config file if provided
        if args.config:
            if not os.path.isfile(args.config):
                errors.append(f"Config file '{args.config}' does not exist")
            elif not args.config.lower().endswith(('.yml', '.yaml')):
                errors.append(f"Config file '{args.config}' should be a YAML file (.yml or .yaml)")
    
    def _validate_simulation_parameters(self, args, errors):
        """Validate simulation timing and step parameters."""
        if not args.simulated_annealing and not args.restart:
            # Check simulation duration
            if args.clock is not None and args.steps is not None:
                errors.append("Cannot specify both --steps and --clock. Choose one simulation duration method")
            elif args.clock is None and args.steps is None:
                errors.append("Must specify simulation duration using either --steps <number> or --clock <minutes>")
            
            # Validate step count
            if args.steps is not None:
                if args.steps <= 0:
                    errors.append(f"Number of steps must be positive, got {args.steps}")
                elif args.steps < 1000:
                    errors.append(f"Warning: Very short simulation ({args.steps} steps). Consider at least 1000 steps")
            
            # Validate clock time
            if args.clock is not None:
                if args.clock <= 0:
                    errors.append(f"Clock time must be positive, got {args.clock} minutes")
                elif args.clock < 1:
                    errors.append(f"Warning: Very short simulation ({args.clock} min). Consider at least 1 minute")
        
        # Validate equilibration steps
        if args.equilibration_steps < 0:
            errors.append(f"Equilibration steps cannot be negative, got {args.equilibration_steps}")
        elif args.equilibration_steps > 10000:
            errors.append(f"Warning: Very long equilibration ({args.equilibration_steps} steps). Typical range is 100-1000")
    
    def _validate_solvation_settings(self, args, errors):
        """Validate solvation and force field settings."""
        if not args.restart:
            # Check solvation method selection
            if (args.solvate and args.GBIS) or (not args.solvate and not args.GBIS):
                errors.append("Must choose exactly one solvation method: either --solvate OR --GBIS")
            
            # Validate solvation-specific parameters
            if args.solvate:
                if args.padding <= 0:
                    errors.append(f"Solvent padding must be positive, got {args.padding} Å")
                elif args.padding < 5:
                    errors.append(f"Warning: Small solvent padding ({args.padding} Å). Consider at least 10 Å")
                
                if args.ionic_strength < 0:
                    errors.append(f"Ionic strength cannot be negative, got {args.ionic_strength} M")
                elif args.ionic_strength > 2.0:
                    errors.append(f"Warning: Very high ionic strength ({args.ionic_strength} M). Typical range is 0.05-0.5 M")
    
    def _validate_restart_settings(self, args, errors):
        """Validate restart-specific settings."""
        if args.restart:
            if not os.path.isdir(args.restart):
                errors.append(f"Restart directory '{args.restart}' does not exist")
            else:
                # Check for required restart files
                restart_setup_file = os.path.join(args.restart, 'restart_setup.yml')
                restart_model_file = os.path.join(args.restart, 'restart_model.pdb')
                last_state_file = os.path.join(args.restart, 'last_state.xml')
                
                if not os.path.isfile(restart_setup_file):
                    errors.append(f"Missing restart setup file: {restart_setup_file}")
                if not os.path.isfile(restart_model_file):
                    errors.append(f"Missing restart model file: {restart_model_file}")
                if not os.path.isfile(last_state_file):
                    errors.append(f"Missing restart state file: {last_state_file}")
    
    def _validate_output_settings(self, args, errors):
        """Validate output directory and reporting settings."""
        # Check output directory
        if args.outdir:
            outdir_path = Path(args.outdir)
            if outdir_path.exists() and not outdir_path.is_dir():
                errors.append(f"Output path '{args.outdir}' exists but is not a directory")
            
            # Check if we can create the directory
            try:
                outdir_path.mkdir(parents=True, exist_ok=True)
            except PermissionError:
                errors.append(f"Cannot create output directory '{args.outdir}': Permission denied")
            except Exception as e:
                errors.append(f"Cannot create output directory '{args.outdir}': {str(e)}")
        
        # Validate reporting interval
        if args.interval <= 0:
            errors.append(f"Reporting interval must be positive, got {args.interval}")
        elif args.steps and args.interval > args.steps:
            errors.append(f"Reporting interval ({args.interval}) cannot be larger than total steps ({args.steps})")
    
    def _validate_physical_parameters(self, args, errors):
        """Validate physical simulation parameters."""
        # Temperature validation
        if args.temperature <= 0:
            errors.append(f"Temperature must be positive, got {args.temperature} K")
        elif args.temperature < 200:
            errors.append(f"Warning: Very low temperature ({args.temperature} K). Typical range is 250-400 K")
        elif args.temperature > 500:
            errors.append(f"Warning: Very high temperature ({args.temperature} K). Consider if this is intended")
        
        # Step size validation
        if args.step_size <= 0:
            errors.append(f"Step size must be positive, got {args.step_size} ps")
        elif args.step_size > 0.004:
            errors.append(f"Warning: Large step size ({args.step_size} ps). Consider 0.001-0.002 ps for stability")
        elif args.step_size < 0.0005:
            errors.append(f"Warning: Very small step size ({args.step_size} ps). This will be very slow")
        
        # Friction coefficient validation
        if args.friction_coeff <= 0:
            errors.append(f"Friction coefficient must be positive, got {args.friction_coeff} /ps")
        elif args.friction_coeff > 10:
            errors.append(f"Warning: Very high friction coefficient ({args.friction_coeff} /ps). Typical range is 0.5-5 /ps")
        
        # pH validation
        if args.ph < 0 or args.ph > 14:
            errors.append(f"pH must be between 0 and 14, got {args.ph}")
        
        # Ion validation
        valid_positive_ions = ['Na+', 'K+', 'Li+', 'Mg2+', 'Ca2+']
        valid_negative_ions = ['Cl-', 'Br-', 'I-', 'F-']
        
        if args.positive_ion not in valid_positive_ions:
            errors.append(f"Unsupported positive ion '{args.positive_ion}'. Supported: {', '.join(valid_positive_ions)}")
        
        if args.negative_ion not in valid_negative_ions:
            errors.append(f"Unsupported negative ion '{args.negative_ion}'. Supported: {', '.join(valid_negative_ions)}")
    
    def _validate_force_fields(self, args, errors):
        """Validate force field specifications."""
        # Common force field files that should exist
        common_protein_ffs = [
            'amber14-all.xml', 'amber99sb-ildn.xml', 'amber03.xml', 
            'charmm36.xml', 'charmm27.xml'
        ]
        
        common_water_ffs = [
            'amber/tip3p_standard.xml', 'amber/tip4pew_standard.xml',
            'amber/spce_standard.xml', 'charmm36/water.xml'
        ]
        
        # Validate protein force field
        if args.protein_force_field not in common_protein_ffs:
            errors.append(f"Warning: Uncommon protein force field '{args.protein_force_field}'. "
                         f"Common options: {', '.join(common_protein_ffs[:3])}")
        
        # Validate water force field for solvated systems
        if args.solvate and args.water_force_field not in common_water_ffs:
            errors.append(f"Warning: Uncommon water force field '{args.water_force_field}'. "
                         f"Common options: {', '.join(common_water_ffs[:3])}")
        
        # Validate ligand force field
        valid_ligand_ffs = ['openff-2.2.0', 'openff-2.1.0', 'openff-2.0.0', 'gaff-2.11', 'gaff-1.81']
        if args.ligand and args.ligand_force_field not in valid_ligand_ffs:
            errors.append(f"Warning: Uncommon ligand force field '{args.ligand_force_field}'. "
                         f"Recommended: {', '.join(valid_ligand_ffs[:3])}")
    
    def _validate_ligand_parameters(self, args, errors):
        """Validate ligand-specific parameters."""
        if args.ligand:
            # Check ligand name format
            if len(args.ligand) > 4:
                errors.append(f"Warning: Ligand name '{args.ligand}' is longer than 4 characters. "
                             "PDB format typically uses 3-letter codes")
            
            # Check if ligand should be removed
            if args.ligand in args.remove:
                errors.append(f"Conflicting settings: Ligand '{args.ligand}' is specified but also in removal list. "
                             f"Remove '{args.ligand}' from --remove list")
        
        # Validate molecules to remove
        if not args.remove:
            errors.append("Warning: No molecules specified for removal. Consider removing common artifacts like DMS, DMSO")
        
        common_artifacts = ['DMS', 'DMSO', 'EDO', 'PEG', 'GOL', 'SO4', 'PO4']
        for molecule in args.remove:
            if len(molecule) > 4:
                errors.append(f"Warning: Molecule name '{molecule}' is longer than 4 characters")
    
    def _validate_missing_residue_parameters(self, args, errors):
        """Validate missing residue handling parameters."""
        # Validate max terminal residues
        if args.max_terminal_residues < 0:
            errors.append(f"Maximum terminal residues cannot be negative, got {args.max_terminal_residues}")
        elif args.max_terminal_residues > 20:
            errors.append(f"Warning: Very high max terminal residues ({args.max_terminal_residues}). "
                         "Consider if this is intended - typical range is 1-10")
        
        # Validate terminal residue types
        common_terminal_types = ['ACE', 'NME', 'NH2', 'COOH', 'FOR', 'NH3+', 'COO-']
        for res_type in args.terminal_residue_types:
            if len(res_type) > 4:
                errors.append(f"Warning: Terminal residue type '{res_type}' is longer than 4 characters")
            if res_type not in common_terminal_types:
                errors.append(f"Warning: Uncommon terminal residue type '{res_type}'. "
                             f"Common types: {', '.join(common_terminal_types[:4])}")
        
        # Validate missing residue strategy combinations
        if args.missing_residues == "none" and (args.skip_missing_loops or args.conservative_missing):
            errors.append("Warning: Missing residue options (--skip-missing-loops, --conservative-missing) "
                         "have no effect when --missing-residues is set to 'none'")
        
        if args.missing_residues == "terminal-only" and args.skip_missing_loops:
            errors.append("Warning: --skip-missing-loops has no effect when --missing-residues is 'terminal-only'")
        
        # Provide guidance on missing residue strategies
        if args.missing_residues == "auto":
            errors.append("Warning: Using 'auto' missing residue strategy. "
                         "Consider specifying explicit strategy for reproducible results")
    
    def _display_warnings(self, warnings):
        """Display warnings but allow simulation to continue."""
        print("\n" + "="*60)
        print("⚠️  VALIDATION WARNINGS")
        print("="*60)
        print("\nThe following issues were detected but won't prevent simulation:")
        print("(Consider reviewing these for optimal results)\n")
        
        for i, warning in enumerate(warnings, 1):
            # Remove "Warning: " prefix for cleaner display
            clean_warning = warning.replace("Warning: ", "")
            print(f"⚠️  {i}. {clean_warning}")
        
        print(f"\nFound {len(warnings)} warning(s). Simulation will continue...")
        print("="*60)
    
    def _display_validation_errors(self, actual_errors, warnings=None):
        """Display validation errors in a user-friendly format."""
        print("\n" + "="*60)
        print("❌ INPUT VALIDATION FAILED")
        print("="*60)
        print("\nPlease fix the following issues before running the simulation:\n")
        
        for i, error in enumerate(actual_errors, 1):
            print(f"❌ {i}. {error}")
        
        warning_count = len(warnings) if warnings else 0
        print(f"\nSummary: {len(actual_errors)} errors, {warning_count} warnings")
        
        print("\n" + "="*60)
        print("QUICK FIXES:")
        print("="*60)
        self._provide_quick_fixes(actual_errors)
        
        print("\n" + "="*60)
        print("For detailed help, run: python -m EasyMD --help")
        print("="*60)
    
    def _provide_quick_fixes(self, errors):
        """Provide quick fix suggestions for common errors."""
        suggestions = []
        
        for error in errors:
            if "Protein PDB file is required" in error:
                suggestions.append("• Add protein file: --protein your_protein.pdb")
            elif "does not exist" in error and "protein" in error.lower():
                suggestions.append("• Check protein file path and ensure file exists")
            elif "Choose one simulation duration" in error:
                suggestions.append("• Use either --steps 10000 OR --clock 60 (not both)")
            elif "Must specify simulation duration" in error:
                suggestions.append("• Add simulation time: --steps 10000 or --clock 60")
            elif "exactly one solvation method" in error:
                suggestions.append("• Choose solvation: --solvate (explicit) OR --GBIS (implicit)")
            elif "restart directory" in error and "does not exist" in error:
                suggestions.append("• Check restart directory path or create new simulation")
        
        if suggestions:
            for suggestion in set(suggestions):  # Remove duplicates
                print(suggestion)
        else:
            print("• Check the error messages above for specific guidance")
            print("• Verify all file paths exist and are accessible")
            print("• Ensure parameter values are within reasonable ranges")
    
    def _display_validation_success(self, args, warnings=None):
        """Display validation success message with simulation summary."""
        if warnings:
            print("\n" + "="*60)
            print("✅ INPUT VALIDATION COMPLETED WITH WARNINGS")
            print("="*60)
            print(f"\nValidation passed with {len(warnings)} warning(s) (see above)")
        else:
            print("\n" + "="*60)
            print("✅ INPUT VALIDATION SUCCESSFUL")
            print("="*60)
        
        print("\nSimulation Configuration Summary:")
        print("-" * 40)
        
        if args.restart:
            print(f"Mode: Restart from {args.restart}")
        elif args.simulated_annealing:
            print("Mode: Simulated Annealing")
        else:
            print("Mode: Standard MD Simulation")
        
        if not args.restart:
            print(f"Protein: {args.protein}")
            if args.ligand:
                print(f"Ligand: {args.ligand}")
            
            if args.solvate:
                print(f"Solvation: Explicit ({args.water_model}, padding={args.padding}Å)")
            else:
                print("Solvation: Implicit (GBIS)")
        
        if args.steps:
            print(f"Duration: {args.steps:,} steps")
        elif args.clock:
            print(f"Duration: {args.clock} minutes")
        
        print(f"Temperature: {args.temperature} K")
        print(f"Output: {args.outdir or 'auto-generated directory'}")
        
        print("\n" + "="*60)
        if warnings:
            print("Proceeding with simulation despite warnings...")
        else:
            print("Ready to start simulation!")
        print("="*60 + "\n")

    def restart_setup(self, setup_file):
        # Temporary parse to check if --config was provided
    
        with open(setup_file, "r") as f:
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
    



