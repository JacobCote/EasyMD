import pytest
import argparse
import tempfile
import os
from unittest.mock import patch, MagicMock
from io import StringIO
from EasyMD.argManager.manager import ArgManager


class TestValidationDisplay:
    """Test suite for validation display methods (warnings, errors, success messages)"""
    
    @pytest.fixture
    def basic_parser(self):
        """Create a basic ArgumentParser for testing"""
        return argparse.ArgumentParser()
    
    @pytest.fixture
    def temp_pdb_file(self):
        """Create a temporary PDB file for testing"""
        pdb_content = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)
    
    @pytest.fixture
    def arg_manager_instance(self, basic_parser):
        """Create an ArgManager instance for testing display methods"""
        with patch('sys.argv', ['test', '--help']):  # Use help to avoid validation
            try:
                return ArgManager(basic_parser)
            except SystemExit:
                # Help causes SystemExit, but we can still create the instance
                manager = ArgManager.__new__(ArgManager)
                manager.parser = basic_parser
                return manager

    def test_display_warnings_format(self, arg_manager_instance):
        """Test that warnings are displayed in correct format"""
        warnings = [
            "Warning: Very short simulation (100 steps). Consider at least 10,000 steps",
            "Warning: High temperature (450 K) is outside typical range"
        ]
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_warnings(warnings)
            
            # Check that print was called multiple times
            assert mock_print.call_count > 0
            
            # Check that warning symbols and formatting are used
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert '⚠️' in printed_output
            assert 'VALIDATION WARNINGS' in printed_output
            assert 'simulation will continue' in printed_output.lower()
    
    def test_display_validation_errors_format(self, arg_manager_instance):
        """Test that errors are displayed in correct format"""
        errors = [
            "Cannot specify both --steps and --clock. Choose one simulation duration method",
            "Protein PDB file is required. Use --protein <file.pdb>"
        ]
        warnings = ["Warning: Short simulation"]
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_errors(errors, warnings)
            
            # Check that print was called
            assert mock_print.call_count > 0
            
            # Check that error symbols and formatting are used
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert '❌' in printed_output
            assert 'INPUT VALIDATION FAILED' in printed_output
            assert 'QUICK FIXES' in printed_output
    
    def test_display_validation_success_no_warnings(self, arg_manager_instance, temp_pdb_file):
        """Test success display without warnings"""
        # Create mock args
        mock_args = MagicMock()
        mock_args.restart = None
        mock_args.simulated_annealing = False
        mock_args.protein = temp_pdb_file
        mock_args.ligand = None
        mock_args.solvate = True
        mock_args.water_model = 'tip3p'
        mock_args.padding = 10
        mock_args.steps = 10000
        mock_args.clock = None
        mock_args.temperature = 300
        mock_args.outdir = None
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_success(mock_args)
            
            # Check that success message is displayed
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert '✅' in printed_output
            assert 'INPUT VALIDATION SUCCESSFUL' in printed_output
            assert 'Ready to start simulation!' in printed_output
    
    def test_display_validation_success_with_warnings(self, arg_manager_instance, temp_pdb_file):
        """Test success display with warnings"""
        # Create mock args
        mock_args = MagicMock()
        mock_args.restart = None
        mock_args.simulated_annealing = False
        mock_args.protein = temp_pdb_file
        mock_args.ligand = 'LIG'
        mock_args.solvate = False  # GBIS
        mock_args.steps = 1000
        mock_args.clock = None
        mock_args.temperature = 300
        mock_args.outdir = 'test_output'
        
        warnings = ["Warning: Short simulation"]
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_success(mock_args, warnings)
            
            # Check that success with warnings message is displayed
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert '✅' in printed_output
            assert 'COMPLETED WITH WARNINGS' in printed_output
            assert 'Proceeding with simulation despite warnings' in printed_output
    
    def test_provide_quick_fixes_suggestions(self, arg_manager_instance):
        """Test that quick fixes provide relevant suggestions"""
        errors = [
            "Protein PDB file is required. Use --protein <file.pdb>",
            "Cannot specify both --steps and --clock. Choose one simulation duration method",
            "Must choose exactly one solvation method: --solvate OR --GBIS"
        ]
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._provide_quick_fixes(errors)
            
            # Check that relevant suggestions are provided
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert 'protein file' in printed_output.lower()
            assert '--steps' in printed_output or '--clock' in printed_output
            assert 'solvate' in printed_output.lower() or 'gbis' in printed_output.lower()
    
    def test_warning_message_cleaning(self, arg_manager_instance):
        """Test that 'Warning: ' prefix is removed from display"""
        warnings = [
            "Warning: This is a test warning message",
            "Warning: Another warning message"
        ]
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_warnings(warnings)
            
            # Check that "Warning: " prefix is removed in display
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            # Should contain the message but not the "Warning: " prefix in the numbered list
            assert 'This is a test warning message' in printed_output
            assert 'Another warning message' in printed_output
    
    def test_error_count_display(self, arg_manager_instance):
        """Test that error and warning counts are displayed correctly"""
        errors = ["Error 1", "Error 2"]
        warnings = ["Warning: Warning 1"]
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_errors(errors, warnings)
            
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert '2 errors' in printed_output
            assert '1 warnings' in printed_output
    
    def test_configuration_summary_restart_mode(self, arg_manager_instance):
        """Test configuration summary for restart mode"""
        mock_args = MagicMock()
        mock_args.restart = '/path/to/restart'
        mock_args.simulated_annealing = False
        mock_args.steps = 1000  # Provide actual integer for formatting
        mock_args.temperature = 300  # Provide actual integer for formatting
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_success(mock_args)
            
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert 'Restart from /path/to/restart' in printed_output
    
    def test_configuration_summary_annealing_mode(self, arg_manager_instance):
        """Test configuration summary for simulated annealing mode"""
        mock_args = MagicMock()
        mock_args.restart = None
        mock_args.simulated_annealing = True
        mock_args.steps = 1000  # Provide actual integer for formatting
        mock_args.temperature = 300  # Provide actual integer for formatting
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_success(mock_args)
            
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert 'Simulated Annealing' in printed_output
    
    def test_configuration_summary_explicit_solvation(self, arg_manager_instance, temp_pdb_file):
        """Test configuration summary shows explicit solvation details"""
        mock_args = MagicMock()
        mock_args.restart = None
        mock_args.simulated_annealing = False
        mock_args.protein = temp_pdb_file
        mock_args.ligand = None
        mock_args.solvate = True
        mock_args.water_model = 'tip4pew'
        mock_args.padding = 15
        mock_args.steps = 50000
        mock_args.clock = None
        mock_args.temperature = 310
        mock_args.outdir = 'my_simulation'
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_success(mock_args)
            
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert 'Explicit (tip4pew, padding=15Å)' in printed_output
            assert '50,000 steps' in printed_output
            assert '310 K' in printed_output
            assert 'my_simulation' in printed_output
    
    def test_configuration_summary_implicit_solvation(self, arg_manager_instance, temp_pdb_file):
        """Test configuration summary shows implicit solvation"""
        mock_args = MagicMock()
        mock_args.restart = None
        mock_args.simulated_annealing = False
        mock_args.protein = temp_pdb_file
        mock_args.ligand = 'ATP'
        mock_args.solvate = False  # GBIS
        mock_args.clock = 120
        mock_args.steps = None
        mock_args.temperature = 298
        mock_args.outdir = None
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_success(mock_args)
            
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert 'Implicit (GBIS)' in printed_output
            assert 'ATP' in printed_output
            assert '120 minutes' in printed_output
            assert 'auto-generated directory' in printed_output
    
    @pytest.mark.parametrize("error_type,expected_suggestion", [
        ("Protein PDB file is required", "Add protein file"),
        ("Choose one simulation duration", "either --steps"),
        ("exactly one solvation method", "Choose solvation"),
        ("restart directory does not exist", "Check restart directory"),
    ])
    def test_specific_quick_fix_suggestions(self, arg_manager_instance, error_type, expected_suggestion):
        """Test that specific error types generate appropriate quick fix suggestions"""
        errors = [error_type]
        
        with patch('builtins.print') as mock_print:
            arg_manager_instance._provide_quick_fixes(errors)
            
            printed_output = ' '.join([str(call) for call in mock_print.call_args_list])
            assert expected_suggestion.lower() in printed_output.lower()
    
    def test_display_methods_handle_empty_lists(self, arg_manager_instance):
        """Test that display methods handle empty error/warning lists gracefully"""
        # Test empty warnings
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_warnings([])
            # Should handle empty list without errors
        
        # Test empty errors
        with patch('builtins.print') as mock_print:
            arg_manager_instance._display_validation_errors([], [])
            # Should handle empty lists without errors
        
        # Test empty quick fixes
        with patch('builtins.print') as mock_print:
            arg_manager_instance._provide_quick_fixes([])
            # Should handle empty list without errors