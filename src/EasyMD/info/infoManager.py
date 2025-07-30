"""
Info Argument Manager for EasyMD

This module handles command-line arguments and configuration for PDB structure analysis.
"""

import argparse
import os
import sys
from pathlib import Path


class InfoManager:
    """
    Manages command-line arguments and configuration for PDB structure analysis.
    
    This class handles argument parsing with comprehensive validation for PDB
    analysis parameters, ensuring that all required files exist and options are valid.
    
    Attributes:
        parser (argparse.ArgumentParser): The argument parser instance.
        args (argparse.Namespace): Parsed arguments after validation.
    """
    
    def __init__(self, parser):
        """
        Initialize the InfoManager with an ArgumentParser instance.
        
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
        Initialize all command-line arguments for PDB structure analysis.
        
        This method defines all available command-line arguments organized into logical groups:
        - Input/output parameters (PDB file, output options)
        - Analysis options (what information to extract)
        - Output formatting options
        """
        # Configuration and version arguments
        self.parser.add_argument("--version", action="version", version="EasyMD Info 1.0.0", 
                                help="Show program version and exit")
        
        # Input/Output parameters
        io_group = self.parser.add_argument_group('Input/Output', 'PDB file and output options')
        io_group.add_argument("pdb_file", type=str,
                             help="Path to PDB file to analyze")
        io_group.add_argument("-o", "--output", type=str, default=None,
                             help="Output file for detailed report (default: <pdb_name>.info)")
        io_group.add_argument("--no-file", action='store_true',
                             help="Don't create output file, only display to terminal")
        
        # Analysis options
        analysis_group = self.parser.add_argument_group('Analysis Options', 'Control what information to extract')
        analysis_group.add_argument("--chains", action='store_true', default=True,
                                   help="Analyze chain information (default: enabled)")
        analysis_group.add_argument("--missing-residues", action='store_true', default=True,
                                   help="Identify missing residues (default: enabled)")
        analysis_group.add_argument("--ligands", action='store_true', default=True,
                                   help="Identify ligands and small molecules (default: enabled)")
        analysis_group.add_argument("--water", action='store_true', default=True,
                                   help="Analyze water molecules (default: enabled)")
        analysis_group.add_argument("--disulfide", action='store_true', default=True,
                                   help="Identify potential disulfide bonds (default: enabled)")
        analysis_group.add_argument("--metals", action='store_true', default=True,
                                   help="Identify metal ions (default: enabled)")
        analysis_group.add_argument("--modifications", action='store_true', default=True,
                                   help="Identify modified residues (default: enabled)")
        analysis_group.add_argument("--all", action='store_true',
                                   help="Perform all available analyses (default behavior)")
        
        # Analysis parameters
        params_group = self.parser.add_argument_group('Analysis Parameters', 'Control analysis behavior')
        params_group.add_argument("--disulfide-distance", type=float, default=2.5,
                                 help="Maximum distance for disulfide bond detection in Angstroms (default: 2.5)")
        params_group.add_argument("--water-threshold", type=int, default=10,
                                 help="Minimum number of water molecules to report details (default: 10)")
        params_group.add_argument("--missing-threshold", type=int, default=3,
                                 help="Minimum gap size to report as missing residues (default: 3)")
        
        # Output formatting
        format_group = self.parser.add_argument_group('Output Formatting', 'Control output appearance')
        format_group.add_argument("--format", choices=["detailed", "summary", "json"], default="detailed",
                                 help="Output format (default: detailed)")
        format_group.add_argument("--no-color", action='store_true',
                                 help="Disable colored output in terminal")
        format_group.add_argument("--quiet", action='store_true',
                                 help="Minimal output, only show warnings and errors")
    
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
        
        # Validate PDB file
        self._validate_pdb_file(args, errors)
        
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
    
    def _validate_pdb_file(self, args, errors):
        """Validate PDB file existence and format."""
        pdb_path = Path(args.pdb_file)
        
        if not pdb_path.exists():
            errors.append(f"PDB file '{args.pdb_file}' does not exist")
            return
        
        if not pdb_path.is_file():
            errors.append(f"'{args.pdb_file}' is not a file")
            return
        
        # Check file extension
        valid_extensions = ['.pdb', '.pdb.gz', '.ent', '.ent.gz']
        if not any(args.pdb_file.lower().endswith(ext) for ext in valid_extensions):
            errors.append(f"File '{args.pdb_file}' does not appear to be a PDB file. "
                         f"Expected extensions: {', '.join(valid_extensions)}")
        
        # Check file size
        file_size = pdb_path.stat().st_size
        if file_size == 0:
            errors.append(f"PDB file '{args.pdb_file}' is empty")
        elif file_size > 100 * 1024 * 1024:  # 100 MB
            errors.append(f"Warning: PDB file '{args.pdb_file}' is very large ({file_size / (1024*1024):.1f} MB). "
                         "Analysis may take a long time")
    
    def _validate_analysis_parameters(self, args, errors, warnings):
        """Validate analysis parameters."""
        # Validate disulfide distance
        if args.disulfide_distance <= 0:
            errors.append(f"Disulfide distance must be positive, got {args.disulfide_distance}")
        elif args.disulfide_distance > 5.0:
            warnings.append(f"Large disulfide distance ({args.disulfide_distance} Å) may detect false positives")
        elif args.disulfide_distance < 1.5:
            warnings.append(f"Small disulfide distance ({args.disulfide_distance} Å) may miss valid bonds")
        
        # Validate water threshold
        if args.water_threshold < 0:
            errors.append(f"Water threshold must be non-negative, got {args.water_threshold}")
        elif args.water_threshold > 1000:
            warnings.append(f"High water threshold ({args.water_threshold}) may hide important information")
        
        # Validate missing residue threshold
        if args.missing_threshold < 1:
            errors.append(f"Missing residue threshold must be at least 1, got {args.missing_threshold}")
        elif args.missing_threshold > 20:
            warnings.append(f"High missing residue threshold ({args.missing_threshold}) may hide small gaps")
    
    def _validate_output_options(self, args, errors, warnings):
        """Validate output options."""
        # Set default output file if not specified
        if args.output is None and not args.no_file:
            pdb_path = Path(args.pdb_file)
            args.output = str(pdb_path.with_suffix('.info'))
        
        # Check if we can create the output file
        if not args.no_file and args.output:
            output_path = Path(args.output)
            try:
                # Try to create the file to check permissions
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.touch(exist_ok=True)
            except PermissionError:
                errors.append(f"Cannot create output file '{args.output}': Permission denied")
            except Exception as e:
                errors.append(f"Cannot create output file '{args.output}': {str(e)}")
        
        # Check for conflicting options
        if args.quiet and args.format == "detailed":
            warnings.append("--quiet flag may not show much with detailed format")
    
    def _display_warnings(self, warnings):
        """Display warnings but allow analysis to continue."""
        print("\n" + "="*60)
        print("⚠️  INFO ANALYSIS WARNINGS")
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
        print("❌ INFO ANALYSIS VALIDATION FAILED")
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
        print("For detailed help, run: python -m EasyMD info --help")
        print("="*60)
    
    def _provide_quick_fixes(self, errors):
        """Provide quick fix suggestions for common errors."""
        suggestions = []
        
        for error in errors:
            if "does not exist" in error:
                suggestions.append("• Check the PDB file path and ensure the file exists")
            elif "not a PDB file" in error:
                suggestions.append("• Ensure the file has a valid PDB extension (.pdb, .pdb.gz, .ent)")
            elif "is empty" in error:
                suggestions.append("• Check that the PDB file contains valid structure data")
            elif "Permission denied" in error:
                suggestions.append("• Check file permissions and write access to output directory")
            elif "must be positive" in error:
                suggestions.append("• Use positive values for distance parameters")
        
        if suggestions:
            for suggestion in set(suggestions):  # Remove duplicates
                print(suggestion)
        else:
            print("• Check the error messages above for specific guidance")
            print("• Verify the PDB file path and format")
            print("• Ensure parameter values are within reasonable ranges")
    
    def _display_validation_success(self, args, warnings=None):
        """Display validation success message with analysis summary."""
        if warnings:
            print("\n" + "="*60)
            print("✅ INFO ANALYSIS VALIDATION COMPLETED WITH WARNINGS")
            print("="*60)
            print(f"\nValidation passed with {len(warnings)} warning(s) (see above)")
        else:
            print("\n" + "="*60)
            print("✅ INFO ANALYSIS VALIDATION SUCCESSFUL")
            print("="*60)
        
        print("\nInfo Analysis Configuration Summary:")
        print("-" * 40)
        
        print(f"PDB File: {args.pdb_file}")
        if not args.no_file:
            print(f"Output File: {args.output}")
        else:
            print("Output File: Terminal only")
        
        print(f"Output Format: {args.format}")
        
        # Show analysis options
        analyses = []
        if args.all or (args.chains and args.missing_residues and args.ligands and 
                       args.water and args.disulfide and args.metals and args.modifications):
            analyses.append("All analyses")
        else:
            if args.chains:
                analyses.append("Chain analysis")
            if args.missing_residues:
                analyses.append("Missing residues")
            if args.ligands:
                analyses.append("Ligands")
            if args.water:
                analyses.append("Water molecules")
            if args.disulfide:
                analyses.append("Disulfide bonds")
            if args.metals:
                analyses.append("Metal ions")
            if args.modifications:
                analyses.append("Modified residues")
        
        print(f"Selected Analyses: {', '.join(analyses) if analyses else 'None'}")
        
        print("\n" + "="*60)
        if warnings:
            print("Proceeding with analysis despite warnings...")
        else:
            print("Ready to start PDB analysis!")
        print("="*60 + "\n")