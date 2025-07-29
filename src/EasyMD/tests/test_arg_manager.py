import pytest
import argparse
import tempfile
import os
import yaml
from unittest.mock import patch, mock_open
from EasyMD.argManager.manager import ArgManager


class TestArgManager:
    """Test suite for ArgManager class"""
    
    @pytest.fixture
    def basic_parser(self):
        """Create a basic ArgumentParser for testing"""
        return argparse.ArgumentParser()
    
    @pytest.fixture
    def sample_config_data(self):
        """Sample configuration data for testing"""
        return {
            'protein': 'test.pdb',
            'ligand': 'LIG',
            'steps': 1000,
            'temperature': 300,
            'solvate': True
        }
    
    @pytest.fixture
    def temp_config_file(self, sample_config_data):
        """Create a temporary config file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
            yaml.dump(sample_config_data, f)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)
    
    def test_init_without_config(self, basic_parser):
        """Test ArgManager initialization without config file"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000', '--solvate']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.protein == 'test.pdb'
            assert args.steps == 1000
            assert args.solvate is True
    
    def test_init_with_config_file(self, basic_parser, temp_config_file):
        """Test ArgManager initialization with config file"""
        with patch('sys.argv', ['test', '--config', temp_config_file]):
            manager = ArgManager(basic_parser)
    
    def test_keep_water_argument_default(self, basic_parser):
        """Test that keep_water argument has correct default value"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000', '--solvate']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.keep_water is False
    
    def test_keep_water_argument_flag(self, basic_parser):
        """Test that --keep-water flag sets keep_water to True"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000', '--solvate', '--keep-water']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.keep_water is True
            assert args.protein == 'test.pdb'
            assert args.steps == 1000
            assert args.solvate is True
    
    def test_config_file_override_by_cli(self, basic_parser, temp_config_file):
        """Test that CLI arguments override config file values"""
        with patch('sys.argv', ['test', '--config', temp_config_file, '--steps', '2000']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.steps == 2000  # CLI should override config
            assert args.protein == 'test.pdb'  # Config value should remain
    
    def test_default_values(self, basic_parser):
        """Test default argument values"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000', '--solvate']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.step_size == 0.002
            assert args.friction_coeff == 1
            assert args.interval == 1000
            assert args.temperature == 300
            assert args.padding == 10
            assert args.water_model == "tip3p"
    
    def test_sanity_check_both_steps_and_clock(self, basic_parser):
        """Test sanity check fails when both steps and clock are provided"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000', '--clock', '60', '--solvate']):
            with pytest.raises(SystemExit):
                ArgManager(basic_parser)
    
    def test_sanity_check_neither_steps_nor_clock(self, basic_parser):
        """Test sanity check fails when neither steps nor clock are provided"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--solvate']):
            with pytest.raises(SystemExit):
                ArgManager(basic_parser)
    
    def test_sanity_check_both_solvate_and_gbis(self, basic_parser):
        """Test sanity check fails when both solvate and GBIS are provided"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000', '--solvate', '--GBIS']):
            with pytest.raises(SystemExit):
                ArgManager(basic_parser)
    
    def test_sanity_check_neither_solvate_nor_gbis(self, basic_parser):
        """Test sanity check fails when neither solvate nor GBIS are provided"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000']):
            with pytest.raises(SystemExit):
                ArgManager(basic_parser)
    

    @pytest.mark.parametrize("water_model", ["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"])
    def test_water_model_choices(self, basic_parser, water_model):
        """Test valid water model choices"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000', '--solvate', '--water-model', water_model]):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.water_model == water_model
    
    def test_remove_molecules_list(self, basic_parser):
        """Test remove molecules argument accepts multiple values"""
        with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '1000', '--solvate', '--remove', 'DMS', 'LIG', 'WAT']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.remove == ['DMS', 'LIG', 'WAT']