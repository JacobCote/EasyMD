"""
Test suite for EasyMD info functionality.
"""

import pytest
import argparse
import tempfile
import os
import sys
from unittest.mock import Mock, patch, MagicMock, mock_open
from pathlib import Path

# Add the test directory to path for importing test utilities
sys.path.append(os.path.dirname(__file__))
from test_data_utils import get_test_pdb_path, get_4zgm_pdb_path

from EasyMD.info.infoManager import InfoManager
from EasyMD.info.infoRunner import InfoRunner


class TestInfoManager:
    """Test suite for InfoManager class"""
    
    @pytest.fixture
    def basic_parser(self):
        """Create a basic ArgumentParser for testing"""
        return argparse.ArgumentParser(add_help=False)
    
    @pytest.fixture
    def temp_pdb_file(self):
        """Create a temporary PDB file for testing"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            # Write minimal PDB content
            f.write("HEADER    TEST PROTEIN\n")
            f.write("ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N\n")
            f.write("ATOM      2  CA  ALA A   1      18.729  16.334  14.618  1.00 20.00           C\n")
            f.write("END\n")
            temp_path = f.name
        
        yield temp_path
        
        # Cleanup
        try:
            os.unlink(temp_path)
        except OSError:
            pass
    
    def test_info_manager_initialization(self, basic_parser, temp_pdb_file):
        """Test InfoManager initialization with valid PDB file"""
        with patch('sys.argv', ['test', temp_pdb_file]):
            manager = InfoManager(basic_parser)
            args = manager.get_args()
            assert args.pdb_file == temp_pdb_file
            assert args.url is None
    
    def test_info_manager_url_option(self, basic_parser):
        """Test InfoManager initialization with URL option"""
        with patch('sys.argv', ['test', '--url', '1ABC']):
            manager = InfoManager(basic_parser)
            args = manager.get_args()
            assert args.url == '1ABC'
            assert args.pdb_file == '1abc.pdb'  # Should be set by validation
    
    def test_info_manager_url_validation_invalid_length(self, basic_parser):
        """Test URL validation with invalid PDB code length"""
        with patch('sys.argv', ['test', '--url', '1AB']):
            with pytest.raises(SystemExit):
                InfoManager(basic_parser)
    
    def test_info_manager_url_validation_invalid_characters(self, basic_parser):
        """Test URL validation with invalid characters"""
        with patch('sys.argv', ['test', '--url', '1AB@']):
            with pytest.raises(SystemExit):
                InfoManager(basic_parser)
    
    def test_info_manager_both_file_and_url_error(self, basic_parser, temp_pdb_file):
        """Test that providing both file and URL raises error"""
        with patch('sys.argv', ['test', temp_pdb_file, '--url', '1ABC']):
            with pytest.raises(SystemExit):
                InfoManager(basic_parser)
    
    def test_info_manager_neither_file_nor_url_error(self, basic_parser):
        """Test that providing neither file nor URL raises error"""
        with patch('sys.argv', ['test']):
            with pytest.raises(SystemExit):
                InfoManager(basic_parser)
    
    def test_missing_pdb_file_error(self, basic_parser):
        """Test error when PDB file doesn't exist"""
        with patch('sys.argv', ['test', '/nonexistent/file.pdb']):
            with pytest.raises(SystemExit):
                InfoManager(basic_parser)
    
    def test_invalid_file_extension_warning(self, basic_parser):
        """Test warning for invalid file extension"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("test content")
            temp_path = f.name
        
        try:
            with patch('sys.argv', ['test', temp_path]):
                with pytest.raises(SystemExit):
                    InfoManager(basic_parser)
        finally:
            os.unlink(temp_path)
    
    def test_output_file_generation(self, basic_parser, temp_pdb_file):
        """Test automatic output file generation"""
        with patch('sys.argv', ['test', temp_pdb_file]):
            manager = InfoManager(basic_parser)
            args = manager.get_args()
            expected_output = temp_pdb_file.replace('.pdb', '.info')
            assert args.output == expected_output
    
    def test_custom_output_file(self, basic_parser, temp_pdb_file):
        """Test custom output file specification"""
        custom_output = 'custom_analysis.info'
        with patch('sys.argv', ['test', temp_pdb_file, '--output', custom_output]):
            manager = InfoManager(basic_parser)
            args = manager.get_args()
            assert args.output == custom_output
    
    def test_no_file_option(self, basic_parser, temp_pdb_file):
        """Test --no-file option"""
        with patch('sys.argv', ['test', temp_pdb_file, '--no-file']):
            manager = InfoManager(basic_parser)
            args = manager.get_args()
            assert args.no_file is True
    
    def test_format_choices(self, basic_parser, temp_pdb_file):
        """Test format parameter choices"""
        valid_formats = ['detailed', 'summary', 'json']
        
        for fmt in valid_formats:
            parser = argparse.ArgumentParser(add_help=False)
            with patch('sys.argv', ['test', temp_pdb_file, '--format', fmt]):
                manager = InfoManager(parser)
                args = manager.get_args()
                assert args.format == fmt
    
    def test_analysis_parameters_validation(self, temp_pdb_file):
        """Test validation of analysis parameters"""
        # Test negative disulfide distance
        parser1 = argparse.ArgumentParser(add_help=False)
        with patch('sys.argv', ['test', temp_pdb_file, '--disulfide-distance', '-1.0']):
            with pytest.raises(SystemExit):
                InfoManager(parser1)
        
        # Test negative water threshold
        parser2 = argparse.ArgumentParser(add_help=False)
        with patch('sys.argv', ['test', temp_pdb_file, '--water-threshold', '-5']):
            with pytest.raises(SystemExit):
                InfoManager(parser2)
    
    def test_analysis_options(self, basic_parser, temp_pdb_file):
        """Test analysis option flags"""
        with patch('sys.argv', ['test', temp_pdb_file, '--all']):
            manager = InfoManager(basic_parser)
            args = manager.get_args()
            assert args.all is True


class TestInfoRunner:
    """Test suite for InfoRunner class"""
    
    @pytest.fixture
    def mock_config(self, tmp_path):
        """Create a mock configuration object"""
        config = Mock()
        config.pdb_file = str(tmp_path / 'test.pdb')
        config.output = str(tmp_path / 'test.info')
        config.no_file = False
        config.format = 'detailed'
        config.no_color = True
        config.quiet = False
        config.chains = True
        config.missing_residues = True
        config.ligands = True
        config.water = True
        config.disulfide = True
        config.metals = True
        config.modifications = True
        config.all = False
        config.disulfide_distance = 2.5
        config.water_threshold = 10
        config.missing_threshold = 3
        return config
    
    @pytest.fixture
    def sample_pdb_content(self):
        """Sample PDB file content for testing"""
        return """HEADER    TEST PROTEIN                            01-JAN-20   TEST
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N
ATOM      2  CA  ALA A   1      18.729  16.334  14.618  1.00 20.00           C
ATOM      3  C   ALA A   1      18.618  14.849  14.897  1.00 20.00           C
ATOM      4  O   ALA A   1      19.633  14.168  15.086  1.00 20.00           O
ATOM      5  CB  ALA A   1      17.951  16.634  13.339  1.00 20.00           C
ATOM      6  N   CYS A   2      17.378  14.378  14.924  1.00 20.00           N
ATOM      7  CA  CYS A   2      17.178  12.967  15.194  1.00 20.00           C
ATOM      8  C   CYS A   2      16.234  12.298  14.194  1.00 20.00           C
ATOM      9  O   CYS A   2      15.234  12.789  13.694  1.00 20.00           O
ATOM     10  CB  CYS A   2      16.578  12.789  16.584  1.00 20.00           C
ATOM     11  SG  CYS A   2      17.634  13.234  17.984  1.00 20.00           S
HETATM   12  O   HOH A 101      15.123  10.456  12.345  1.00 30.00           O
HETATM   13 MG    MG A 201      20.456  18.789  16.123  1.00 25.00          MG
HETATM   14  C1  LIG A 301      14.567  11.234  13.789  1.00 35.00           C
HETATM   15  C2  LIG A 301      13.234  10.567  14.123  1.00 35.00           C
END"""
    
    def test_info_runner_initialization(self, mock_config):
        """Test InfoRunner initialization"""
        runner = InfoRunner(mock_config)
        assert runner.config == mock_config
        assert str(runner.pdb_file) == mock_config.pdb_file
        assert str(runner.output_file) == mock_config.output
    
    @patch('builtins.open')
    @patch('gzip.open')
    def test_parse_pdb_file(self, mock_gzip_open, mock_open, mock_config, sample_pdb_content):
        """Test PDB file parsing"""
        # Setup the mock to return our sample content line by line
        lines = sample_pdb_content.split('\n')
        
        # Create a mock file object that can be iterated
        mock_file = MagicMock()
        mock_file.__enter__.return_value = lines
        mock_file.__exit__.return_value = None
        mock_open.return_value = mock_file
        
        runner = InfoRunner(mock_config)
        runner._parse_pdb_file()
        
        # Check that structure data was populated
        assert len(runner.structure_data['atoms']) > 0
        assert len(runner.structure_data['hetero_atoms']) > 0
        assert len(runner.structure_data['chains']) > 0
        assert 'A' in runner.structure_data['chains']
    
    def test_analyze_chains(self, mock_config):
        """Test chain analysis"""
        runner = InfoRunner(mock_config)
        
        # Mock structure data
        runner.structure_data = {
            'chains': {
                'A': {
                    'atoms': [
                        {'record_type': 'ATOM', 'residue_id': 1, 'residue_name': 'ALA', 'insertion_code': ''},
                        {'record_type': 'ATOM', 'residue_id': 2, 'residue_name': 'CYS', 'insertion_code': ''}
                    ]
                }
            }
        }
        
        result = runner._analyze_chains()
        
        assert 'A' in result
        assert result['A']['protein_residues'] == 2
        assert result['A']['sequence'] == ['ALA', 'CYS']
    
    def test_analyze_missing_residues(self, mock_config):
        """Test missing residue analysis"""
        runner = InfoRunner(mock_config)
        
        # Mock analysis results with a gap
        runner.analysis_results = {
            'chains': {
                'A': {
                    'residue_details': [
                        {'id': 1, 'name': 'ALA'},
                        {'id': 5, 'name': 'CYS'}  # Gap from 2-4
                    ]
                }
            }
        }
        
        result = runner._analyze_missing_residues()
        
        assert 'A' in result
        assert len(result['A']['gaps']) == 1
        assert result['A']['gaps'][0]['count'] == 3
        assert result['A']['gaps'][0]['start'] == 2
        assert result['A']['gaps'][0]['end'] == 4
    
    def test_analyze_ligands(self, mock_config):
        """Test ligand analysis"""
        runner = InfoRunner(mock_config)
        
        # Mock hetero atoms
        runner.structure_data = {
            'hetero_atoms': [
                {'residue_name': 'LIG', 'chain_id': 'A', 'residue_id': 301, 'element': 'C'},
                {'residue_name': 'HOH', 'chain_id': 'A', 'residue_id': 101, 'element': 'O'},
                {'residue_name': 'MG', 'chain_id': 'A', 'residue_id': 201, 'element': 'MG'}
            ]
        }
        
        result = runner._analyze_ligands()
        
        assert result['ligand_count'] == 1
        assert 'A_LIG_301' in result['ligands']
        assert result['ligands']['A_LIG_301']['name'] == 'LIG'
    
    def test_analyze_disulfide_bonds(self, mock_config):
        """Test disulfide bond analysis"""
        runner = InfoRunner(mock_config)
        
        # Mock atoms with two close cysteine sulfurs
        runner.structure_data = {
            'atoms': [
                {'element': 'S', 'atom_name': 'SG', 'residue_name': 'CYS', 
                 'chain_id': 'A', 'residue_id': 1, 'x': 0.0, 'y': 0.0, 'z': 0.0},
                {'element': 'S', 'atom_name': 'SG', 'residue_name': 'CYS', 
                 'chain_id': 'A', 'residue_id': 10, 'x': 2.0, 'y': 0.0, 'z': 0.0}
            ]
        }
        
        result = runner._analyze_disulfide_bonds()
        
        assert result['cysteine_count'] == 2
        assert result['bond_count'] == 1
        assert result['bonds'][0]['distance'] == 2.0
    
    def test_analyze_metal_ions(self, mock_config):
        """Test metal ion analysis"""
        runner = InfoRunner(mock_config)
        
        # Mock hetero atoms with metals
        runner.structure_data = {
            'hetero_atoms': [
                {'element': 'MG', 'residue_name': 'MG', 'chain_id': 'A', 'residue_id': 201,
                 'x': 1.0, 'y': 2.0, 'z': 3.0, 'occupancy': 1.0, 'b_factor': 20.0},
                {'element': 'ZN', 'residue_name': 'ZN', 'chain_id': 'A', 'residue_id': 202,
                 'x': 4.0, 'y': 5.0, 'z': 6.0, 'occupancy': 1.0, 'b_factor': 25.0}
            ]
        }
        
        result = runner._analyze_metal_ions()
        
        assert result['total_count'] == 2
        assert 'MG' in result['element_counts']
        assert 'ZN' in result['element_counts']
        assert result['element_counts']['MG'] == 1
        assert result['element_counts']['ZN'] == 1
    
    def test_analyze_water_molecules(self, mock_config):
        """Test water molecule analysis"""
        runner = InfoRunner(mock_config)
        
        # Mock water molecules
        runner.structure_data = {
            'water_molecules': [
                {'chain_id': 'A', 'x': 1.0, 'y': 2.0, 'z': 3.0},
                {'chain_id': 'A', 'x': 4.0, 'y': 5.0, 'z': 6.0},
                {'chain_id': 'B', 'x': 7.0, 'y': 8.0, 'z': 9.0}
            ]
        }
        
        result = runner._analyze_water_molecules()
        
        assert result['total_count'] == 3
        assert result['chain_counts']['A'] == 2
        assert result['chain_counts']['B'] == 1
    
    @patch('builtins.open', new_callable=mock_open)
    def test_save_results(self, mock_file, mock_config):
        """Test saving results to file"""
        runner = InfoRunner(mock_config)
        runner.structure_data = {'header_info': {'pdb_id': 'TEST'}}
        runner.analysis_results = {'chains': {'A': {'total_atoms': 10}}}
        
        runner._save_results()
        
        # Check that file was opened for writing
        mock_file.assert_called_once()
        # Check that write was called (content written to file)
        handle = mock_file.return_value
        assert handle.write.called
    
    @patch('urllib.request.urlopen')
    def test_download_pdb_file_success(self, mock_urlopen, mock_config):
        """Test successful PDB file download"""
        # Setup mock config with URL
        mock_config.url = '1ABC'
        mock_config.no_color = True
        
        # Mock the HTTP response with context manager support
        mock_response = Mock()
        mock_response.getcode.return_value = 200
        mock_response.read.return_value = b"HEADER    TEST PROTEIN\nATOM      1  N   ALA A   1\nEND\n"
        mock_response.__enter__ = Mock(return_value=mock_response)
        mock_response.__exit__ = Mock(return_value=None)
        mock_urlopen.return_value = mock_response
        
        runner = InfoRunner(mock_config)
        
        with patch('builtins.open', mock_open()) as mock_file:
            runner._download_pdb_file()
            
            # Check that the URL was called correctly
            expected_url = "https://files.rcsb.org/download/1ABC.pdb"
            mock_urlopen.assert_called_once_with(expected_url)
            
            # Check that file was written
            mock_file.assert_called_once_with('1abc.pdb', 'w')
            handle = mock_file.return_value
            handle.write.assert_called_once()
            
            # Check that pdb_file path was updated
            assert str(runner.pdb_file) == '1abc.pdb'
    
    @patch('urllib.request.urlopen')
    def test_download_pdb_file_not_found(self, mock_urlopen, mock_config):
        """Test PDB file download with 404 error"""
        mock_config.url = 'XXXX'
        
        # Mock 404 error
        from urllib.error import HTTPError
        mock_urlopen.side_effect = HTTPError(
            url="https://files.rcsb.org/download/XXXX.pdb",
            code=404,
            msg="Not Found",
            hdrs=None,
            fp=None
        )
        
        runner = InfoRunner(mock_config)
        
        with pytest.raises(RuntimeError, match="PDB code 'XXXX' not found"):
            runner._download_pdb_file()
    
    @patch('urllib.request.urlopen')
    def test_download_pdb_file_network_error(self, mock_urlopen, mock_config):
        """Test PDB file download with network error"""
        mock_config.url = '1ABC'
        
        # Mock network error
        from urllib.error import URLError
        mock_urlopen.side_effect = URLError("Network unreachable")
        
        runner = InfoRunner(mock_config)
        
        with pytest.raises(RuntimeError, match="Network error"):
            runner._download_pdb_file()


class TestInfoIntegration:
    """Integration tests for info functionality"""
    
    def test_info_help_command(self):
        """Test that info help command works"""
        with patch('sys.argv', ['test', 'dummy.pdb', '--help']):
            with patch('sys.stdout'), patch('sys.stderr'):
                try:
                    from EasyMD.info.infoManager import InfoManager
                    parser = argparse.ArgumentParser(add_help=True)
                    InfoManager(parser)
                except SystemExit as e:
                    # Help command should exit with code 0
                    assert e.code == 0
    
    @patch('EasyMD.info.infoRunner.InfoRunner.run')
    @patch('EasyMD.info.infoManager.InfoManager.__init__')
    def test_run_info_function(self, mock_manager_init, mock_runner_run):
        """Test the run_info function"""
        from EasyMD.__main__ import run_info
        
        # Mock the manager initialization to return without validation
        mock_manager_init.return_value = None
        mock_manager = Mock()
        
        # Create a proper mock config with string paths
        mock_config = Mock()
        mock_config.pdb_file = '/tmp/test.pdb'
        mock_config.output = '/tmp/test.info'
        mock_config.no_color = True
        mock_manager.get_args.return_value = mock_config
        
        with patch('EasyMD.__main__.InfoManager', return_value=mock_manager):
            with patch('sys.argv', ['test']):  # Minimal argv
                run_info()
        
        # Check that the runner was called
        mock_runner_run.assert_called_once()
    
    def test_main_function_info_routing(self):
        """Test that main function correctly routes to info"""
        with patch('EasyMD.__main__.run_info') as mock_run_info:
            with patch('sys.argv', ['test', 'info', 'some_file.pdb']):
                from EasyMD.__main__ import main
                main()
                mock_run_info.assert_called_once()