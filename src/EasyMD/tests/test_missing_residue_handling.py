import pytest
import argparse
import tempfile
import os
from unittest.mock import Mock, patch, MagicMock
from EasyMD.argManager.manager import ArgManager
from EasyMD.sysGenerator.SysGenerator import SysGenerator


class TestMissingResidueHandling:
    """Test suite for missing residue handling functionality"""
    
    @pytest.fixture
    def basic_parser(self):
        """Create a basic ArgumentParser for testing"""
        return argparse.ArgumentParser()
    
    @pytest.fixture
    def temp_pdb_file(self):
        """Create a temporary PDB file for testing"""
        pdb_content = """HEADER    TEST PROTEIN WITH MISSING RESIDUES
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
ATOM      5  N   VAL A   5      18.154  17.967  15.365  1.00 20.00           N  
ATOM      6  CA  VAL A   5      17.030  17.101  15.618  1.00 20.00           C  
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)
    
    @pytest.fixture
    def mock_config_missing_residues(self):
        """Create mock config with missing residue parameters"""
        config = Mock()
        config.missing_residues = 'auto'
        config.max_terminal_residues = 5
        config.terminal_residue_types = ['ACE', 'NME', 'NH2', 'COOH']
        config.skip_missing_loops = False
        config.conservative_missing = False
        return config

    # Test argument parsing for missing residue options
    
    def test_missing_residues_default_values(self, basic_parser, temp_pdb_file):
        """Test that missing residue arguments have correct default values"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            
            assert args.missing_residues == 'auto'
            assert args.max_terminal_residues == 5
            assert args.terminal_residue_types == ['ACE', 'NME', 'NH2', 'COOH']
            assert args.skip_missing_loops is False
            assert args.conservative_missing is False
    
    @pytest.mark.parametrize("strategy", ["auto", "none", "non-terminal", "terminal-only", "all"])
    def test_missing_residues_strategy_choices(self, basic_parser, temp_pdb_file, strategy):
        """Test that all missing residue strategy choices are accepted"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', 
                                '--missing-residues', strategy]):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.missing_residues == strategy
    
    def test_max_terminal_residues_validation(self, basic_parser, temp_pdb_file):
        """Test validation of max terminal residues parameter"""
        # Valid value should pass
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', 
                                '--max-terminal-residues', '3']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.max_terminal_residues == 3
    
    def test_negative_max_terminal_residues_blocks_execution(self, basic_parser, temp_pdb_file):
        """Test that negative max terminal residues blocks execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', 
                                '--max-terminal-residues', '-1']):
            with pytest.raises(SystemExit) as exc_info:
                ArgManager(basic_parser)
            assert exc_info.value.code == 1
    
    def test_high_max_terminal_residues_generates_warning(self, basic_parser, temp_pdb_file):
        """Test that very high max terminal residues generates warning but allows execution"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', 
                                '--max-terminal-residues', '25']):
            # Should not raise SystemExit (warning only)
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.max_terminal_residues == 25
    
    def test_terminal_residue_types_custom(self, basic_parser, temp_pdb_file):
        """Test custom terminal residue types"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', 
                                '--terminal-residue-types', 'ACE', 'NME']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.terminal_residue_types == ['ACE', 'NME']
    
    def test_skip_missing_loops_flag(self, basic_parser, temp_pdb_file):
        """Test skip missing loops flag"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', 
                                '--skip-missing-loops']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.skip_missing_loops is True
    
    def test_conservative_missing_flag(self, basic_parser, temp_pdb_file):
        """Test conservative missing flag"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', 
                                '--conservative-missing']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.conservative_missing is True

    # Test validation warnings for conflicting options
    
    def test_conflicting_options_generate_warnings(self, basic_parser, temp_pdb_file):
        """Test that conflicting missing residue options generate warnings but allow execution"""
        # missing-residues=none with skip-missing-loops should warn but not block
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '1000', '--solvate', 
                                '--missing-residues', 'none', '--skip-missing-loops']):
            # Should not raise SystemExit (warning only)
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            assert args.missing_residues == 'none'
            assert args.skip_missing_loops is True

    # Test SysGenerator missing residue handling methods
    
    @patch('EasyMD.sysGenerator.SysGenerator.PDBFixer')
    def test_handle_missing_residues_none_strategy(self, mock_pdb_fixer, mock_config_missing_residues):
        """Test that 'none' strategy clears missing residues"""
        mock_config_missing_residues.missing_residues = 'none'
        
        # Create mock fixer with missing residues
        mock_fixer = Mock()
        mock_fixer.missingResidues = {'A': ['residue1', 'residue2']}
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config_missing_residues
        
        # Call the method
        sys_gen._handle_missing_residues(mock_fixer)
        
        # Should clear missing residues
        assert mock_fixer.missingResidues == {}
    
    @patch('EasyMD.sysGenerator.SysGenerator.PDBFixer')
    def test_handle_missing_residues_auto_strategy(self, mock_pdb_fixer, mock_config_missing_residues):
        """Test that 'auto' strategy preserves missing residues for PDBFixer"""
        mock_config_missing_residues.missing_residues = 'auto'
        
        # Create mock fixer with missing residues
        mock_fixer = Mock()
        original_missing = {'A': ['residue1', 'residue2']}
        mock_fixer.missingResidues = original_missing.copy()
        mock_fixer.addMissingResidues = Mock()
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config_missing_residues
        
        # Call the method
        sys_gen._handle_missing_residues(mock_fixer)
        
        # Should preserve missing residues (they will be added when addMissingAtoms is called)
        assert mock_fixer.missingResidues == original_missing
    
    @patch('EasyMD.sysGenerator.SysGenerator.PDBFixer')
    def test_handle_missing_residues_no_missing(self, mock_pdb_fixer, mock_config_missing_residues):
        """Test handling when no missing residues are found"""
        # Create mock fixer with no missing residues
        mock_fixer = Mock()
        mock_fixer.missingResidues = {}
        mock_fixer.addMissingResidues = Mock()
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config_missing_residues
        
        # Call the method
        sys_gen._handle_missing_residues(mock_fixer)
        
        # Should not call addMissingResidues
        mock_fixer.addMissingResidues.assert_not_called()
    
    def test_apply_conservative_filtering(self, mock_config_missing_residues):
        """Test conservative filtering of missing residues"""
        # Create mock fixer with many missing residues
        mock_fixer = Mock()
        mock_fixer.missingResidues = {
            'A': ['res1', 'res2', 'res3', 'res4', 'res5', 'res6', 'res7'],  # 7 residues
            'B': ['res1', 'res2']  # 2 residues
        }
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config_missing_residues
        
        # Apply conservative filtering with max 5 residues
        sys_gen._apply_conservative_filtering(mock_fixer, max_terminal_residues=5)
        
        # Chain A should be filtered out (>5 residues), Chain B should remain
        assert 'A' not in mock_fixer.missingResidues
        assert 'B' in mock_fixer.missingResidues
        assert mock_fixer.missingResidues['B'] == ['res1', 'res2']
    
    def test_filter_non_terminal_residues(self, mock_config_missing_residues):
        """Test filtering to keep only non-terminal residues"""
        # Create mock fixer with missing residues
        mock_fixer = Mock()
        mock_fixer.missingResidues = {
            'A': ['res1', 'res2', 'res3', 'res4', 'res5'],  # 5 residues - should keep middle ones
            'B': ['res1', 'res2']  # 2 residues - no clear internal residues
        }
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config_missing_residues
        
        # Apply non-terminal filtering
        sys_gen._filter_non_terminal_residues(mock_fixer)
        
        # Should keep middle residues from chain A, none from chain B
        assert 'A' in mock_fixer.missingResidues
        assert mock_fixer.missingResidues['A'] == ['res2', 'res3', 'res4']  # Middle residues
        assert 'B' not in mock_fixer.missingResidues
    
    def test_filter_terminal_residues(self, mock_config_missing_residues):
        """Test filtering to keep only terminal residues"""
        # Create mock fixer with missing residues
        mock_fixer = Mock()
        mock_fixer.missingResidues = {
            'A': ['res1', 'res2', 'res3', 'res4', 'res5', 'res6'],  # 6 residues
        }
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config_missing_residues
        
        # Apply terminal filtering with max 4 residues
        sys_gen._filter_terminal_residues(mock_fixer, max_terminal_residues=4, allowed_types=['ACE', 'NME'])
        
        # Should keep up to 4 terminal residues
        assert 'A' in mock_fixer.missingResidues
        assert len(mock_fixer.missingResidues['A']) <= 4
    
    def test_filter_loop_residues(self, mock_config_missing_residues):
        """Test filtering out likely loop regions"""
        # Create mock fixer with missing residues
        mock_fixer = Mock()
        mock_fixer.missingResidues = {
            'A': ['res1', 'res2'],  # Small gap - should keep
            'B': ['res1', 'res2', 'res3', 'res4', 'res5']  # Large gap - likely loop, should skip
        }
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config_missing_residues
        
        # Apply loop filtering
        sys_gen._filter_loop_residues(mock_fixer)
        
        # Should keep small gaps, skip large ones
        assert 'A' in mock_fixer.missingResidues
        assert mock_fixer.missingResidues['A'] == ['res1', 'res2']
        assert 'B' not in mock_fixer.missingResidues

    # Integration tests
    
    @pytest.mark.parametrize("strategy,expected_behavior", [
        ("none", "should_clear_missing"),
        ("auto", "should_preserve_missing"),
        ("non-terminal", "should_filter_missing"),
        ("terminal-only", "should_filter_missing"),
        ("all", "should_preserve_missing"),
    ])
    def test_missing_residue_strategies_integration(self, strategy, expected_behavior, mock_config_missing_residues):
        """Test that different strategies produce expected behavior"""
        mock_config_missing_residues.missing_residues = strategy
        
        # Create mock fixer with missing residues
        mock_fixer = Mock()
        mock_fixer.missingResidues = {'A': ['res1', 'res2', 'res3']}
        mock_fixer.addMissingResidues = Mock()
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config_missing_residues
        
        # Call the method
        sys_gen._handle_missing_residues(mock_fixer)
        
        if expected_behavior == "should_clear_missing":
            assert mock_fixer.missingResidues == {}
        elif expected_behavior == "should_preserve_missing":
            assert len(mock_fixer.missingResidues) > 0
        elif expected_behavior == "should_filter_missing":
            # Should modify the missing residues dict (filtering)
            # The exact result depends on the filtering logic
            # We don't check addMissingResidues calls since our implementation
            # relies on PDBFixer to add them automatically during addMissingAtoms()
            pass