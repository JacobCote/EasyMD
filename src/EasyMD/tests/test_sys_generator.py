import pytest
import tempfile
import os
from unittest.mock import Mock, patch, MagicMock, mock_open
from EasyMD.sysGenerator.sysGenerator import SysGenerator


class TestSysGenerator:
    """Test suite for SysGenerator class"""
    
    @pytest.fixture
    def mock_config_no_restart(self):
        """Create a mock configuration without restart"""
        config = Mock()
        config.restart = None
        config.outdir = None
        config.protein = 'test.pdb'
        config.ligand = None
        config.remove = ['DMS']
        config.solvate = True
        config.protein_force_field = 'amber14-all.xml'
        config.water_force_field = 'amber/tip3p_standard.xml'
        config.water_model = 'tip3p'
        config.positive_ion = 'Na+'
        config.negative_ion = 'Cl-'
        config.ionic_strength = 0.1
        config.no_neutralize = False
        config.padding = 10
        config.ph = 7.0
        config.steps = 1000
        config.interval = 1000
        config.temperature = 300
        config.equilibration_steps = 200
        return config
    
    @pytest.fixture
    def mock_config_restart(self):
        """Create a mock configuration with restart"""
        config = Mock()
        config.restart = 'restart_dir'
        return config
    
    @pytest.fixture
    def mock_config_with_ligand(self):
        """Create a mock configuration with ligand"""
        config = Mock()
        config.restart = None
        config.outdir = 'test_out'
        config.protein = 'test.pdb'
        config.ligand = 'LIG'
        config.ligand_force_field = 'openff-2.2.0'
        config.remove = ['DMS']
        config.solvate = True
        config.protein_force_field = 'amber14-all.xml'
        config.water_force_field = 'amber/tip3p_standard.xml'
        config.water_model = 'tip3p'
        config.positive_ion = 'Na+'
        config.negative_ion = 'Cl-'
        config.ionic_strength = 0.1
        config.no_neutralize = False
        config.padding = 10
        config.ph = 7.0
        return config
    
    @patch('EasyMD.sysGenerator.sysGenerator.os.path.isdir')
    @patch('EasyMD.sysGenerator.sysGenerator.os.mkdir')
    @patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._prep_prot')
    def test_init_no_restart_no_ligand(self, mock_prep_prot, mock_mkdir, mock_isdir, mock_config_no_restart):
        """Test initialization without restart and without ligand"""
        mock_isdir.return_value = False
        mock_modeller = Mock()
        mock_system = Mock()
        mock_prep_prot.return_value = (mock_modeller, mock_system)
        
        sys_gen = SysGenerator(mock_config_no_restart)
        
        assert sys_gen.config == mock_config_no_restart
        assert sys_gen.modeller == mock_modeller
        assert sys_gen.system == mock_system
        assert sys_gen.last_state == 0
        mock_prep_prot.assert_called_once()
    
    @patch('EasyMD.sysGenerator.sysGenerator.os.path.isdir')
    @patch('EasyMD.sysGenerator.sysGenerator.os.mkdir')
    @patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._prep_complex')
    def test_init_no_restart_with_ligand(self, mock_prep_complex, mock_mkdir, mock_isdir, mock_config_with_ligand):
        """Test initialization without restart but with ligand"""
        mock_isdir.return_value = False
        mock_modeller = Mock()
        mock_system = Mock()
        mock_prep_complex.return_value = (mock_modeller, mock_system)
        
        sys_gen = SysGenerator(mock_config_with_ligand)
        
        assert sys_gen.config == mock_config_with_ligand
        assert sys_gen.modeller == mock_modeller
        assert sys_gen.system == mock_system
        mock_prep_complex.assert_called_once()
    
    @patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._restartSetup')
    def test_init_with_restart(self, mock_restart_setup, mock_config_restart):
        """Test initialization with restart"""
        mock_modeller = Mock()
        mock_system = Mock()
        mock_restart_setup.return_value = (mock_modeller, mock_system)
        
        sys_gen = SysGenerator(mock_config_restart)
        
        assert sys_gen.config == mock_config_restart
        assert sys_gen.modeller == mock_modeller
        assert sys_gen.system == mock_system
        mock_restart_setup.assert_called_once()
    
    @patch('EasyMD.sysGenerator.sysGenerator.os.path.exists')
    @patch('EasyMD.sysGenerator.sysGenerator.os.makedirs')
    def test_setup_creates_output_directory(self, mock_makedirs, mock_exists, mock_config_no_restart):
        """Test that setup creates output directory when it doesn't exist"""
        mock_exists.side_effect = [True, True, False]  # out_0, out_1 exist, out_2 doesn't
        
        with patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._prep_prot') as mock_prep:
            mock_prep.return_value = (Mock(), Mock())
            SysGenerator(mock_config_no_restart)
        
        mock_makedirs.assert_called_with('out_2', exist_ok=True)
    
    @patch('EasyMD.sysGenerator.sysGenerator.os.makedirs')
    def test_setup_uses_provided_outdir(self, mock_makedirs, mock_config_no_restart):
        """Test that setup uses provided output directory"""
        mock_config_no_restart.outdir = 'custom_out/'
        
        with patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._prep_prot') as mock_prep:
            mock_prep.return_value = (Mock(), Mock())
            SysGenerator(mock_config_no_restart)
        
        mock_makedirs.assert_called_with('custom_out/', exist_ok=True)
    
    @patch('builtins.open', new_callable=mock_open, read_data='protein_force_field: amber14-all.xml\nsolvate: true')
    @patch('EasyMD.sysGenerator.sysGenerator.yaml.load')
    @patch('EasyMD.sysGenerator.sysGenerator.os.path.isfile')
    @patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._prep_restart')
    def test_restart_setup_without_ligand(self, mock_prep_restart, mock_isfile, mock_yaml_load, mock_file, mock_config_restart):
        """Test restart setup without ligand file"""
        mock_setup = {'protein_force_field': 'amber14-all.xml', 'solvate': True}
        mock_yaml_load.return_value = mock_setup
        mock_isfile.return_value = False  # No ligand.sdf
        mock_modeller = Mock()
        mock_system = Mock()
        mock_prep_restart.return_value = (mock_modeller, mock_system)
        
        sys_gen = SysGenerator(mock_config_restart)
        
        mock_prep_restart.assert_called_once()
        assert sys_gen.modeller == mock_modeller
        assert sys_gen.system == mock_system
    
    @patch('builtins.open', new_callable=mock_open, read_data='protein_force_field: amber14-all.xml\nsolvate: true')
    @patch('EasyMD.sysGenerator.sysGenerator.yaml.load')
    @patch('EasyMD.sysGenerator.sysGenerator.os.path.isfile')
    @patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._prep_restart_ligand')
    def test_restart_setup_with_ligand(self, mock_prep_restart_ligand, mock_isfile, mock_yaml_load, mock_file, mock_config_restart):
        """Test restart setup with ligand file"""
        mock_setup = {'protein_force_field': 'amber14-all.xml', 'solvate': True}
        mock_yaml_load.return_value = mock_setup
        mock_isfile.return_value = True  # ligand.sdf exists
        mock_modeller = Mock()
        mock_system = Mock()
        mock_prep_restart_ligand.return_value = (mock_modeller, mock_system)
        
        sys_gen = SysGenerator(mock_config_restart)
        
        mock_prep_restart_ligand.assert_called_once()
        assert sys_gen.modeller == mock_modeller
        assert sys_gen.system == mock_system
    
    @patch('EasyMD.sysGenerator.sysGenerator.PDBFile')
    @patch('EasyMD.sysGenerator.sysGenerator.Molecule.from_file')
    @patch('EasyMD.sysGenerator.sysGenerator.Modeller')
    @patch('EasyMD.sysGenerator.sysGenerator.SystemGenerator')
    def test_prep_restart_ligand_solvated(self, mock_sys_gen, mock_modeller_class, mock_molecule, mock_pdb):
        """Test _prep_restart_ligand with solvated system"""
        setup = {
            'protein_force_field': 'amber14-all.xml',
            'water_force_field': 'tip3p.xml',
            'ligand_force_field': 'openff-2.2.0',
            'solvate': True
        }
        
        mock_pdb_instance = Mock()
        mock_pdb.return_value = mock_pdb_instance
        mock_ligand = Mock()
        mock_molecule.return_value = mock_ligand
        mock_modeller_instance = Mock()
        mock_modeller_class.return_value = mock_modeller_instance
        mock_system_generator = Mock()
        mock_sys_gen.return_value = mock_system_generator
        mock_system = Mock()
        mock_system_generator.create_system.return_value = mock_system
        
        sys_gen = SysGenerator.__new__(SysGenerator)  # Create without calling __init__
        result = sys_gen._prep_restart_ligand(setup, 'test_dir', {})
        
        assert result == (mock_modeller_instance, mock_system)
        mock_sys_gen.assert_called_once()
        mock_system_generator.create_system.assert_called_once()
    
    @pytest.mark.parametrize("solvate", [True, False])
    def test_forcefield_kwargs_initialization(self, solvate, mock_config_no_restart):
        """Test that forcefield_kwargs are properly initialized"""
        mock_config_no_restart.solvate = solvate
        
        with patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._prep_prot') as mock_prep:
            mock_prep.return_value = (Mock(), Mock())
            sys_gen = SysGenerator(mock_config_no_restart)
        
        expected_kwargs = {
            'constraints': 'HBonds',  # This will be the string representation
            'rigidWater': True,
            'removeCMMotion': False,
            'hydrogenMass': '4.0 amu'  # This will be the string representation
        }
        
        # Check that forcefield_kwargs is a dictionary with expected keys
        assert isinstance(sys_gen.forcefield_kwargs, dict)
        assert 'constraints' in sys_gen.forcefield_kwargs
        assert 'rigidWater' in sys_gen.forcefield_kwargs
        assert 'removeCMMotion' in sys_gen.forcefield_kwargs
        assert 'hydrogenMass' in sys_gen.forcefield_kwargs