"""
Analysis Argument Manager for EasyMD

This module handles command-line arguments and configuration for trajectory analysis.
"""

import argparse
import os
import sys
from pathlib import Path


class AnalysisManager:
    """
    Manages command-line arguments and configuration for trajectory analysis.
    
    This class handles argument parsing with comprehensive validation for analysis
    parameters, ensuring that all required files exist and analysis options are valid.
    
    Attributes:
        parser (argparse.ArgumentParser): The argument parser instance.
        args (argparse.Namespace): Parsed arguments after validation.
    """
    
    def __init__(self, parser):
        """
        Initialize the AnalysisManager with an ArgumentParser instance.
        
        Args:
            parser (argparse.ArgumentParser): The argument parser to configure.
            
        Raises:
            SystemExit: If argument validation fails during initialization.
        """
        self.parser = parser
        self._initialize()
        self.args = self.parser.parse_args()
        self._validate_args()
    
    def _initialize(self):
        """
        Initialize all command-line arguments for trajectory analysis.
        
        This method defines all available command-line arguments organized into logical groups:
        - Input/output parameters (trajectory directory, output options)
        - Analysis types (RMSD, RMSF, distances, etc.)
        - Analysis parameters (reference frame, atom selection, etc.)
        - Output options (plot formats, data export)
        """
        # Configuration and version arguments
        self.parser.add_argument("--version", action="version", version="EasyMD Analysis 1.0.0", 
                                help="Show program version and exit")
        
        # Input/Output parameters
        io_group = self.parser.add_argument_group('Input/Output', 'Trajectory files and output options')
        io_group.add_argument("trajectory_directory", type=str,
                             help="Path to directory containing trajectory files and topology")
        io_group.add_argument("-o", "--output-dir", type=str, default=None,
                             help="Output directory for analysis results (default: trajectory_dir/analysis)")
        io_group.add_argument("--output-format", choices=["png", "pdf", "svg"], default="png",
                             help="Output format for plots (default: png)")
        
        # Analysis types
        analysis_group = self.parser.add_argument_group('Analysis Types', 'Select which analyses to perform')
        analysis_group.add_argument("--rmsd", action='store_true',
                                   help="Calculate Root Mean Square Deviation (RMSD)")
        analysis_group.add_argument("--rmsf", action='store_true', 
                                   help="Calculate Root Mean Square Fluctuation (RMSF)")
        analysis_group.add_argument("--distances", action='store_true',
                                   help="Calculate distances between centers of mass of chains")
        analysis_group.add_argument("--radius-gyration", action='store_true',
                                   help="Calculate radius of gyration")
        analysis_group.add_argument("--secondary-structure", action='store_true',
                                   help="Analyze secondary structure evolution")
        analysis_group.add_argument("--all", action='store_true',
                                   help="Perform all available analyses")
        
        # Analysis parameters
        params_group = self.parser.add_argument_group('Analysis Parameters', 'Control analysis behavior')
        params_group.add_argument("--reference-frame", type=int, default=0,
                                 help="Reference frame for RMSD calculation (default: 0)")
        params_group.add_argument("--atom-selection", choices=["all", "backbone", "ca", "heavy"], default="backbone",
                                 help="Atom selection for RMSD calculation (default: backbone)")
        params_group.add_argument("--skip-frames", type=int, default=1,
                                 help="Skip every N frames for analysis (default: 1, no skipping)")
        params_group.add_argument("--start-frame", type=int, default=0,
                                 help="Starting frame for analysis (default: 0)")
        params_group.add_argument("--end-frame", type=int, default=None,
                                 help="Ending frame for analysis (default: all frames)")
        
        # Output options
        output_group = self.parser.add_argument_group('Output Options', 'Control output generation')
        output_group.add_argument("--save-data", action='store_true',
                                 help="Save analysis data as CSV files")
        output_group.add_argument("--no-plots", action='store_true',
                                 help="Skip plot generation (only save data)")
        output_group.add_argument("--dpi", type=int, default=300,
                                 help="DPI for plot output (default: 300)")
        output_group.add_argument("--plot-style", choices=["default", "seaborn", "ggplot"], default="default",
                                 help="Matplotlib style for plots (default: default)")
    
    def get_args(self):
        """Return the parsed and validated arguments."""
        return self.args
    
    def _validate_args(self):
        """
        Perform comprehensive validation of all input arguments.
        
        This method validates all input parameters and provides clear, actionable
        error messages to help users fix their command-line arguments before
        starting the analysis.
        
        Raises:
            SystemExit: If any validation check fails, with descriptive error message.
        """
        args = self.args
        errors = []
        warnings = []
        
        # Validate trajectory directory
        self._validate_trajectory_directory(args, errors)
        
        # Validate analysis selection
        self._validate_analysis_selection(args, errors, warnings)
        
        # Validate analysis parameters
        self._validate_analysis_parameters(args, errors, warnings)
        
        # Validate output options
        self._validate_output_options(args, errors, warnings)
        
        # Display warnings but continue execution
        if warnings:
            self._display_warnings(warnings)
        
        # Only exit on actual errors
        if errors:
            self._display_validation_errors(errors, warnings)
            sys.exit(1)
        
        # Display validation success message
        self._display_validation_success(args, warnings)
    
    def _validate_trajectory_directory(self, args, errors):
        """Validate trajectory directory and required files."""
        traj_dir = Path(args.trajectory_directory)
        
        if not traj_dir.exists():
            errors.append(f"Trajectory directory '{args.trajectory_directory}' does not exist")
            return
        
        if not traj_dir.is_dir():
            errors.append(f"'{args.trajectory_directory}' is not a directory")
            return
        
        # Check for required files
        topology_file = traj_dir / 'topology.pkl'
        if not topology_file.exists():
            errors.append(f"Missing topology file: {topology_file}")
        
        # Check for trajectory files
        dcd_files = list(traj_dir.glob('output_traj_*.dcd'))
        if not dcd_files:
            errors.append(f"No trajectory files (output_traj_*.dcd) found in {traj_dir}")
        else:
            # Store the number of trajectory files for later use
            args._n_trajectory_files = len(dcd_files)
    
    def _validate_analysis_selection(self, args, errors, warnings):
        """Validate analysis type selection."""
        analysis_types = [args.rmsd, args.rmsf, args.distances, args.radius_gyration, 
                         args.secondary_structure, args.all]
        
        if not any(analysis_types):
            errors.append("No analysis type selected. Use --rmsd, --rmsf, --distances, "
                         "--radius-gyration, --secondary-structure, or --all")
        
        if args.all and any(analysis_types[:-1]):  # If --all is used with other options
            warnings.append("--all flag overrides individual analysis selections")
    
    def _validate_analysis_parameters(self, args, errors, warnings):
        """Validate analysis parameters."""
        # Validate reference frame
        if args.reference_frame < 0:
            errors.append(f"Reference frame must be non-negative, got {args.reference_frame}")
        
        # Validate frame range
        if args.start_frame < 0:
            errors.append(f"Start frame must be non-negative, got {args.start_frame}")
        
        if args.end_frame is not None and args.end_frame <= args.start_frame:
            errors.append(f"End frame ({args.end_frame}) must be greater than start frame ({args.start_frame})")
        
        # Validate skip frames
        if args.skip_frames < 1:
            errors.append(f"Skip frames must be at least 1, got {args.skip_frames}")
        elif args.skip_frames > 10:
            warnings.append(f"Large skip value ({args.skip_frames}) may result in insufficient data points")
    
    def _validate_output_options(self, args, errors, warnings):
        """Validate output directory and options."""
        # Set default output directory if not specified
        if args.output_dir is None:
            args.output_dir = os.path.join(args.trajectory_directory, 'analysis')
        
        # Check if we can create the output directory
        output_path = Path(args.output_dir)
        try:
            output_path.mkdir(parents=True, exist_ok=True)
        except PermissionError:
            errors.append(f"Cannot create output directory '{args.output_dir}': Permission denied")
        except Exception as e:
            errors.append(f"Cannot create output directory '{args.output_dir}': {str(e)}")
        
        # Validate DPI
        if args.dpi < 50 or args.dpi > 1200:
            warnings.append(f"Unusual DPI value ({args.dpi}). Typical range is 150-600")
        
        # Check for conflicting options
        if args.no_plots and not args.save_data:
            warnings.append("--no-plots specified without --save-data. No output will be generated")
    
    def _display_warnings(self, warnings):
        """Display warnings but allow analysis to continue."""
        print("\n" + "="*60)
        print("⚠️  ANALYSIS WARNINGS")
        print("="*60)
        print("\nThe following issues were detected but won't prevent analysis:")
        print("(Consider reviewing these for optimal results)\n")
        
        for i, warning in enumerate(warnings, 1):
            print(f"⚠️  {i}. {warning}")
        
        print(f"\nFound {len(warnings)} warning(s). Analysis will continue...")
        print("="*60)
    
    def _display_validation_errors(self, errors, warnings=None):
        """Display validation errors in a user-friendly format."""
        print("\n" + "="*60)
        print("❌ ANALYSIS VALIDATION FAILED")
        print("="*60)
        print("\nPlease fix the following issues before running analysis:\n")
        
        for i, error in enumerate(errors, 1):
            print(f"❌ {i}. {error}")
        
        warning_count = len(warnings) if warnings else 0
        print(f"\nSummary: {len(errors)} errors, {warning_count} warnings")
        
        print("\n" + "="*60)
        print("QUICK FIXES:")
        print("="*60)
        self._provide_quick_fixes(errors)
        
        print("\n" + "="*60)
        print("For detailed help, run: python -m EasyMD analyze --help")
        print("="*60)
    
    def _provide_quick_fixes(self, errors):
        """Provide quick fix suggestions for common errors."""
        suggestions = []
        
        for error in errors:
            if "does not exist" in error and "directory" in error:
                suggestions.append("• Check the trajectory directory path")
            elif "Missing topology file" in error:
                suggestions.append("• Ensure the simulation completed successfully and topology.pkl exists")
            elif "No trajectory files" in error:
                suggestions.append("• Check that trajectory files (output_traj_*.dcd) exist in the directory")
            elif "No analysis type selected" in error:
                suggestions.append("• Add analysis options: --rmsd, --rmsf, --distances, or --all")
            elif "Reference frame must be" in error:
                suggestions.append("• Use a non-negative reference frame number")
        
        if suggestions:
            for suggestion in set(suggestions):  # Remove duplicates
                print(suggestion)
        else:
            print("• Check the error messages above for specific guidance")
            print("• Verify all file paths exist and are accessible")
            print("• Ensure parameter values are within reasonable ranges")
    
    def _display_validation_success(self, args, warnings=None):
        """Display validation success message with analysis summary."""
        if warnings:
            print("\n" + "="*60)
            print("✅ ANALYSIS VALIDATION COMPLETED WITH WARNINGS")
            print("="*60)
            print(f"\nValidation passed with {len(warnings)} warning(s) (see above)")
        else:
            print("\n" + "="*60)
            print("✅ ANALYSIS VALIDATION SUCCESSFUL")
            print("="*60)
        
        print("\nAnalysis Configuration Summary:")
        print("-" * 40)
        
        print(f"Trajectory Directory: {args.trajectory_directory}")
        print(f"Output Directory: {args.output_dir}")
        
        # Show selected analyses
        analyses = []
        if args.all:
            analyses.append("All analyses")
        else:
            if args.rmsd:
                analyses.append("RMSD")
            if args.rmsf:
                analyses.append("RMSF")
            if args.distances:
                analyses.append("Chain distances")
            if args.radius_gyration:
                analyses.append("Radius of gyration")
            if args.secondary_structure:
                analyses.append("Secondary structure")
        
        print(f"Selected Analyses: {', '.join(analyses)}")
        print(f"Atom Selection: {args.atom_selection}")
        print(f"Reference Frame: {args.reference_frame}")
        
        if hasattr(args, '_n_trajectory_files'):
            print(f"Trajectory Files Found: {args._n_trajectory_files}")
        
        print("\n" + "="*60)
        if warnings:
            print("Proceeding with analysis despite warnings...")
        else:
            print("Ready to start analysis!")
        print("="*60 + "\n")