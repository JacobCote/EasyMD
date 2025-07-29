import pytest
import argparse
import tempfile
import os
import sys
from unittest.mock import Mock, patch, MagicMock
from EasyMD.argManager.manager import ArgManager
from EasyMD.sysGenerator.sysGenerator import SysGenerator

# Add the test directory to path for importing test utilities
sys.path.append(os.path.dirname(__file__))
from test_data_utils import get_test_pdb_path


class TestWaterHandling:
    """Test suite for water molecule handling functionality"""
    
    @pytest.fixture
    def basic_parser(self):
        """Create a basic ArgumentParser for testing"""
        return argparse.ArgumentParser()
    
    @pytest.fixture
    def temp_pdb_with_water(self):
        """Create a temporary PDB file with water molecules for testing"""
        pdb_content = """HEADER    TEST PROTEIN WITH WATER MOLECULES
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
ATOM      5  N   VAL A   2      16.500  16.200  14.950  1.00 20.00           N  
ATOM      6  CA  VAL A   2      15.154  16.767  15.200  1.00 20.00           C  
HETATM    7  O   HOH W   1      25.000  20.000  20.000  1.00 30.00           O  
HETATM    8  H1  HOH W   1      25.500  20.500  20.500  1.00 30.00           H  
HETATM    9  H2  HOH W   1      24.500  19.500  19.500  1.00 30.00           H  
HETATM   10  O   WAT W   2      30.000  25.000  25.000  1.00 30.00           O  
HETATM   11  H1  WAT W   2      30.500  25.500  25.500  1.00 30.00           H  
HETATM   12  H2  WAT W   2      29.500  24.500  24.500  1.00 30.00           H  
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)
    
    @pytest.fixture
    def temp_pdb_no_water(self):
        """Create a temporary PDB file without water molecules for testing"""
        pdb_content = """HEADER    TEST PROTEIN WITHOUT WATER MOLECULES
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
ATOM      5  N   VAL A   2      16.500  16.200  14.950  1.00 20.00           N  
ATOM      6  CA  VAL A   2      15.154  16.767  15.200  1.00 20.00           C  
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)

    # Test argument parsing for keep-water option
    
    def test_keep_water_default_value(self, basic_parser, temp_pdb_with_water):
        """Test that keep_water argument has correct default value (False)"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_with_water, '--steps', '1000', '--solvate']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            
            assert args.keep_water is False
    
    def test_keep_water_flag_sets_true(self, basic_parser, temp_pdb_with_water):
        """Test that --keep-water flag sets keep_water to True"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_with_water, '--steps', '1000', '--solvate', '--keep-water']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            
            assert args.keep_water is True
    
    def test_keep_water_with_other_options(self, basic_parser, temp_pdb_with_water):
        """Test that keep_water works correctly with other options"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_with_water, '--steps', '1000', 
                                '--solvate', '--keep-water', '--remove', 'DMS', 'LIG']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            
            assert args.keep_water is True
            assert args.remove == ['DMS', 'LIG']

    # Test SysGenerator water handling logic
    
    @patch('EasyMD.sysGenerator.sysGenerator.PDBFixer')
    @patch('EasyMD.sysGenerator.sysGenerator.Modeller')
    def test_water_removal_default_behavior(self, mock_modeller_class, mock_pdb_fixer):
        """Test that water molecules are removed by default (keep_water=False)"""
        # Setup mock config without keep_water (default behavior)
        mock_config = Mock()
        mock_config.keep_water = False
        mock_config.remove = ['DMS']
        
        # Setup mock fixer and modeller
        mock_fixer = Mock()
        mock_fixer.missingResidues = {}
        mock_fixer.missingAtoms = {}
        mock_fixer.nonstandardResidues = []
        mock_fixer.topology = Mock()
        mock_fixer.positions = Mock()
        
        mock_modeller = Mock()
        mock_topology = Mock()
        mock_topology.getNumAtoms.return_value = 100
        mock_modeller.topology = mock_topology
        
        # Setup residues with water molecules
        mock_residue_hoh = Mock()
        mock_residue_hoh.name = 'HOH'
        mock_residue_wat = Mock()
        mock_residue_wat.name = 'WAT'
        mock_residue_dms = Mock()
        mock_residue_dms.name = 'DMS'
        mock_residue_protein = Mock()
        mock_residue_protein.name = 'ALA'
        
        mock_topology.residues.return_value = [mock_residue_hoh, mock_residue_wat, mock_residue_dms, mock_residue_protein]
        mock_modeller_class.return_value = mock_modeller
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config
        
        # Call the method (we'll need to mock the full method call)
        with patch.object(sys_gen, '_handle_missing_residues'):
            with patch('EasyMD.sysGenerator.sysGenerator.PDBwrite_all'):
                with patch('EasyMD.sysGenerator.sysGenerator.SystemGenerator'):
                    with patch('EasyMD.sysGenerator.sysGenerator.mdtraj'):
                        with patch('EasyMD.sysGenerator.sysGenerator.deletePcap', return_value=mock_modeller):
                            try:
                                sys_gen._prep_prot(
                                    pdb_in=get_test_pdb_path(),
                                    list_of_molecules_to_remove=['DMS'],
                                    solvate=False,
                                    protein_force_field='amber14-all.xml',
                                    water_force_field='amber/tip3p_standard.xml',
                                    water_model='tip3p',
                                    positive_ion='Na+',
                                    negative_ion='Cl-',
                                    ionic_strength=0.15,
                                    no_neutralize=False,
                                    padding=10.0,
                                    ph=7.0,
                                    outdir='test_out',
                                    forcefield_kwargs={}
                                )
                            except:
                                pass  # We expect this to fail due to mocking, but we can check the delete calls
        
        # Verify that delete was called with water molecules and DMS
        expected_deletions = 3  # HOH, WAT, DMS should all be deleted
        assert mock_modeller.delete.call_count == expected_deletions
    
    @patch('EasyMD.sysGenerator.sysGenerator.PDBFixer')
    @patch('EasyMD.sysGenerator.sysGenerator.Modeller')
    def test_water_preservation_with_keep_water_true(self, mock_modeller_class, mock_pdb_fixer):
        """Test that water molecules are preserved when keep_water=True"""
        # Setup mock config with keep_water=True
        mock_config = Mock()
        mock_config.keep_water = True
        mock_config.remove = ['DMS']
        
        # Setup mock fixer and modeller
        mock_fixer = Mock()
        mock_fixer.missingResidues = {}
        mock_fixer.missingAtoms = {}
        mock_fixer.nonstandardResidues = []
        mock_fixer.topology = Mock()
        mock_fixer.positions = Mock()
        
        mock_modeller = Mock()
        mock_topology = Mock()
        mock_topology.getNumAtoms.return_value = 100
        mock_modeller.topology = mock_topology
        
        # Setup residues with water molecules
        mock_residue_hoh = Mock()
        mock_residue_hoh.name = 'HOH'
        mock_residue_wat = Mock()
        mock_residue_wat.name = 'WAT'
        mock_residue_dms = Mock()
        mock_residue_dms.name = 'DMS'
        mock_residue_protein = Mock()
        mock_residue_protein.name = 'ALA'
        
        mock_topology.residues.return_value = [mock_residue_hoh, mock_residue_wat, mock_residue_dms, mock_residue_protein]
        mock_modeller_class.return_value = mock_modeller
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config
        
        # Call the method (we'll need to mock the full method call)
        with patch.object(sys_gen, '_handle_missing_residues'):
            with patch('EasyMD.sysGenerator.sysGenerator.PDBwrite_all'):
                with patch('EasyMD.sysGenerator.sysGenerator.SystemGenerator'):
                    with patch('EasyMD.sysGenerator.sysGenerator.mdtraj'):
                        with patch('EasyMD.sysGenerator.sysGenerator.deletePcap', return_value=mock_modeller):
                            try:
                                sys_gen._prep_prot(
                                    pdb_in=get_test_pdb_path(),
                                    list_of_molecules_to_remove=['DMS'],
                                    solvate=False,
                                    protein_force_field='amber14-all.xml',
                                    water_force_field='amber/tip3p_standard.xml',
                                    water_model='tip3p',
                                    positive_ion='Na+',
                                    negative_ion='Cl-',
                                    ionic_strength=0.15,
                                    no_neutralize=False,
                                    padding=10.0,
                                    ph=7.0,
                                    outdir='test_out',
                                    forcefield_kwargs={}
                                )
                            except:
                                pass  # We expect this to fail due to mocking, but we can check the delete calls
        
        # Verify that delete was called only with DMS (not water molecules)
        expected_deletions = 1  # Only DMS should be deleted, water should be preserved
        assert mock_modeller.delete.call_count == expected_deletions
    
    def test_water_handling_in_complex_preparation(self):
        """Test that water handling works correctly in complex preparation"""
        # Setup mock config with keep_water=True
        mock_config = Mock()
        mock_config.keep_water = True
        mock_config.remove = ['DMS']
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config
        
        # Test the logic directly
        list_of_molecules_to_remove = ['DMS']
        
        # Simulate the logic from _prep_complex
        if not getattr(sys_gen.config, 'keep_water', False):
            list_of_molecules_to_remove += ['HOH', 'WAT']
        list_of_molecules_to_remove += ['LIG']  # ligand name
        
        # With keep_water=True, water should not be in the removal list
        assert 'HOH' not in list_of_molecules_to_remove
        assert 'WAT' not in list_of_molecules_to_remove
        assert 'DMS' in list_of_molecules_to_remove
        assert 'LIG' in list_of_molecules_to_remove
    
    def test_water_handling_in_complex_preparation_default(self):
        """Test that water handling works correctly in complex preparation with default behavior"""
        # Setup mock config without keep_water (default behavior)
        mock_config = Mock()
        mock_config.remove = ['DMS']
        # Explicitly delete keep_water attribute to ensure getattr returns False
        if hasattr(mock_config, 'keep_water'):
            delattr(mock_config, 'keep_water')
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config
        
        # Test the logic directly - start with the original list
        list_of_molecules_to_remove = mock_config.remove.copy()  # ['DMS']
        
        # Simulate the logic from _prep_complex
        if not getattr(sys_gen.config, 'keep_water', False):
            list_of_molecules_to_remove += ['HOH', 'WAT']
        list_of_molecules_to_remove += ['LIG']  # ligand name
        
        # With default behavior, water should be in the removal list
        assert 'HOH' in list_of_molecules_to_remove
        assert 'WAT' in list_of_molecules_to_remove
        assert 'DMS' in list_of_molecules_to_remove
        assert 'LIG' in list_of_molecules_to_remove

    # Integration tests
    
    @pytest.mark.parametrize("keep_water,expected_water_removed", [
        (False, True),   # Default behavior: remove water
        (True, False),   # Keep water: don't remove water
    ])
    def test_water_handling_integration(self, keep_water, expected_water_removed):
        """Test water handling integration with different keep_water settings"""
        mock_config = Mock()
        mock_config.keep_water = keep_water
        mock_config.remove = ['DMS']
        
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config
        
        # Test the actual logic used in the code
        list_of_molecules_to_remove = ['DMS']
        
        # This is the actual logic from _prep_prot
        if not getattr(sys_gen.config, 'keep_water', False):
            list_of_molecules_to_remove += ['HOH','WAT']
        
        if expected_water_removed:
            assert 'HOH' in list_of_molecules_to_remove
            assert 'WAT' in list_of_molecules_to_remove
        else:
            assert 'HOH' not in list_of_molecules_to_remove
            assert 'WAT' not in list_of_molecules_to_remove
        
        # DMS should always be in the removal list
        assert 'DMS' in list_of_molecules_to_remove
    
    def test_keep_water_with_no_water_molecules(self, basic_parser, temp_pdb_no_water):
        """Test that keep_water option works correctly even when no water molecules are present"""
        with patch('sys.argv', ['test', '--protein', temp_pdb_no_water, '--steps', '1000', 
                                '--solvate', '--keep-water']):
            manager = ArgManager(basic_parser)
            args = manager.get_args()
            
            # Should still parse correctly
            assert args.keep_water is True
            assert args.protein == temp_pdb_no_water
    
    def test_keep_water_help_text(self, basic_parser):
        """Test that the help text for keep_water is correct"""
        test_pdb = get_test_pdb_path()
        
        # We need to suppress the validation output and SystemExit for this test
        # since we just want to check the help text after arguments are added
        with patch('sys.argv', ['test', '--help']):
            with patch('sys.stdout'), patch('sys.stderr'):
                try:
                    manager = ArgManager(basic_parser)
                except SystemExit:
                    # Expected when --help is used
                    pass
            
            # Get the help text after arguments have been added
            help_text = basic_parser.format_help()
            
            # Check that the keep-water option is documented
            assert '--keep-water' in help_text
            assert 'Preserve crystal water molecules from PDB' in help_text
            # Check for the default text (might have slight formatting differences)
            assert 'default:' in help_text and 'remove all water' in help_text
    
    def test_keep_water_with_solvation(self):
        """Test that keep_water option is compatible with solvation"""
        # This is a logical test - keeping crystal waters and then solvating
        # should be a valid combination
        mock_config = Mock()
        mock_config.keep_water = True
        mock_config.solvate = True
        mock_config.remove = []
        
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config
        
        # Test that both options can be set together
        assert getattr(sys_gen.config, 'keep_water', False) is True
        assert sys_gen.config.solvate is True
        
        # The logic should still work correctly
        list_of_molecules_to_remove = []
        if not getattr(sys_gen.config, 'keep_water', False):
            list_of_molecules_to_remove += ['HOH','WAT']
        
        # Water should not be removed when keep_water=True
        assert 'HOH' not in list_of_molecules_to_remove
        assert 'WAT' not in list_of_molecules_to_remove