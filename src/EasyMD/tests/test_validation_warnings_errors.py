import pytest
import argparse
import tempfile
import os
import sys
import yaml
from unittest.mock import patch, mock_open
from EasyMD.argManager.manager import ArgManager

# Add the test directory to path for importing test utilities
sys.path.append(os.path.dirname(__file__))
from test_data_utils import get_test_pdb_path


class TestValidationWarningsErrors:
    """Test suite for validation system warnings vs errors behavior"""
    
    @pytest.fixture
    def basic_parser(self):
        """Create a basic ArgumentParser for testing"""
        return argparse.ArgumentParser()
    
    @pytest.fixture
    def temp_pdb_file(self):
        """Create a temporary PDB file for testing"""
        pdb_content = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)
    
    @pytest.fixture
    def temp_config_with_warnings(self, temp_pdb_file):
        """Create a config file that should generate warnings but not errors"""
        config_data = {
            'protein': temp_pdb_file,
            'steps': 50,  # Very short - should warn
            'temperature': 450,  # High temp - should warn
            'step_size': 0.005,  # Large step - should warn
            'solvate': True
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
            yaml.dump(config_data, f)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)

    # Tests for scenarios that should PASS with warnings (not block execution)
    
    def test_short_simulation_warning_allows_execution(self, basic_parser, temp_pdb_file):
        """Test that short simulation generates warning but allows execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '50', '--interval', '10', '--solvate']):
            # Should not raise SystemExit - warnings don't block execution
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.steps == 50
            assert args.protein == temp_pdb_file
            assert args.solvate is True
    
    def test_high_temperature_warning_allows_execution(self, basic_parser, temp_pdb_file):
        """Test that high temperature generates warning but allows execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--temperature', '450']):
            # Should not raise SystemExit
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.temperature == 450
            assert args.steps == 1000
    
    def test_large_step_size_warning_allows_execution(self, basic_parser, temp_pdb_file):
        """Test that large step size generates warning but allows execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--step-size', '0.005']):
            # Should not raise SystemExit
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.step_size == 0.005
    
    def test_uncommon_force_field_warning_allows_execution(self, basic_parser, temp_pdb_file):
        """Test that uncommon force field generates warning but allows execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--protein-force-field', 'custom.xml']):
            # Should not raise SystemExit
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.protein_force_field == 'custom.xml'
    
    def test_long_ligand_name_warning_allows_execution(self, basic_parser, temp_pdb_file):
        """Test that long ligand name generates warning but allows execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--ligand', 'VERYLONGNAME', '--steps', '1000', '--GBIS']):
            # Should not raise SystemExit
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.ligand == 'VERYLONGNAME'
    
    def test_multiple_warnings_allow_execution(self, basic_parser, temp_pdb_file):
        """Test that multiple warnings don't block execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '50', '--interval', '10', '--solvate', 
                                '--temperature', '450', '--step-size', '0.005']):
            # Should not raise SystemExit even with multiple warnings
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.steps == 50
            assert args.temperature == 450
            assert args.step_size == 0.005
    
    def test_config_file_with_warnings_allows_execution(self, basic_parser, temp_config_with_warnings):
        """Test that config file with warning parameters allows execution"""
        with patch('sys.argv', ['test', '--config', temp_config_with_warnings, '--interval', '10']):
            # Should not raise SystemExit
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.steps == 50  # Short simulation
            assert args.temperature == 450  # High temperature
            assert args.step_size == 0.005  # Large step size

    # Tests for scenarios that should FAIL with errors (block execution)
    
    def test_missing_protein_file_blocks_execution(self, basic_parser):
        """Test that missing protein file blocks execution with error"""
        with patch('sys.argv', ['test', '--steps', '1000', '--solvate']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_conflicting_duration_methods_block_execution(self, basic_parser, temp_pdb_file):
        """Test that conflicting duration methods block execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--clock', '60', '--solvate']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_no_solvation_method_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that missing solvation method blocks execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_both_solvation_methods_block_execution(self, basic_parser, temp_pdb_file):
        """Test that specifying both solvation methods blocks execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--GBIS']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            # Exit code can be 1 (validation error) or 2 (argparse error)
            assert exc_info.value.code in [1, 2]
    
    def test_negative_temperature_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that negative temperature blocks execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--temperature', '-100']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_negative_steps_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that negative steps block execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '-1000', '--solvate']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_nonexistent_protein_file_blocks_execution(self, basic_parser):
        """Test that non-existent protein file blocks execution"""
        with patch('sys.argv', ['test', '--protein', 'nonexistent.pdb', '--steps', '1000', '--solvate']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_nonexistent_restart_directory_blocks_execution(self, basic_parser):
        """Test that non-existent restart directory blocks execution"""
        with patch('sys.argv', ['test', '--restart', '/nonexistent/directory']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_invalid_ph_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that invalid pH blocks execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--ph', '15']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1

    # Tests for mixed scenarios (warnings + errors)
    
    def test_warnings_with_errors_still_block_execution(self, basic_parser, temp_pdb_file):
        """Test that presence of errors blocks execution even with warnings"""
        # Short simulation (warning) + conflicting duration methods (error)
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '50', '--clock', '60', '--solvate']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1

    # Tests for edge cases
    
    def test_zero_steps_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that zero steps block execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '0', '--solvate']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_zero_temperature_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that zero temperature blocks execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--temperature', '0']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_negative_padding_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that negative padding blocks execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--padding', '-5']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_negative_ionic_strength_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that negative ionic strength blocks execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--ionic-strength', '-0.1']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1

    # Tests for boundary conditions
    
    def test_minimum_valid_steps_allows_execution(self, basic_parser, temp_pdb_file):
        """Test that minimum valid steps (1) allows execution but may warn"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1', '--interval', '1', '--solvate']):
            # Should not raise SystemExit (might warn but shouldn't block)
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.steps == 1
    
    def test_minimum_valid_temperature_allows_execution(self, basic_parser, temp_pdb_file):
        """Test that minimum valid temperature allows execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--temperature', '1']):
            # Should not raise SystemExit (might warn but shouldn't block)
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.temperature == 1
    
    def test_boundary_ph_values_allow_execution(self, temp_pdb_file):
        """Test that boundary pH values (0 and 14) allow execution"""
        # Test pH = 0
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--ph', '0']):
            parser = argparse.ArgumentParser()  # Fresh parser
            manager = ArgManager(parser)
            args = manager.get_args()
            assert args.ph == 0.0
        
        # Test pH = 14
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', '--ph', '14']):
            parser = argparse.ArgumentParser()  # Fresh parser
            manager = ArgManager(parser)
            args = manager.get_args()
            assert args.ph == 14.0

    # Integration tests
    
    def test_valid_configuration_with_no_warnings(self, basic_parser, temp_pdb_file):
        """Test that a completely valid configuration passes without warnings or errors"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '10000', '--solvate', 
                                '--temperature', '300', '--step-size', '0.002']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.steps == 10000
            assert args.temperature == 300
            assert args.step_size == 0.002
            assert args.solvate is True
    
    @pytest.mark.parametrize("warning_param,warning_value", [
        ('steps', '100'),  # Short simulation
        ('temperature', '450'),  # High temperature
        ('step-size', '0.005'),  # Large step size
        ('protein-force-field', 'custom.xml'),  # Uncommon force field
    ])
    def test_individual_warning_parameters_allow_execution(self, basic_parser, temp_pdb_file, warning_param, warning_value):
        """Test that individual warning parameters allow execution"""
        base_args = ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate']
        if warning_param != 'steps':  # Don't duplicate steps
            base_args.extend([f'--{warning_param}', warning_value])
        else:
            base_args[base_args.index('1000')] = warning_value
            # Add smaller interval for short simulations to avoid validation error
            base_args.extend(['--interval', '10'])
        
        with patch('sys.argv', base_args):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            # Should not raise SystemExit
            assert args is not None
    
    @pytest.mark.parametrize("error_scenario", [
        (['--steps', '1000', '--solvate']),  # Missing protein
        (['--protein', get_test_pdb_path(), '--steps', '1000', '--clock', '60', '--solvate']),  # Conflicting duration
        (['--protein', get_test_pdb_path(), '--steps', '1000']),  # No solvation method
        (['--protein', get_test_pdb_path(), '--steps', '1000', '--solvate', '--GBIS']),  # Both solvation methods
        (['--protein', get_test_pdb_path(), '--steps', '-100', '--solvate']),  # Negative steps
    ])
    def test_error_scenarios_block_execution(self, basic_parser, error_scenario):
        """Test that various error scenarios block execution"""
        with patch('sys.argv', ['test'] + error_scenario):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            # Exit code can be 1 (validation error) or 2 (argparse error)
            assert exc_info.value.code in [1, 2]