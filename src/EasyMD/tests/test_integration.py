import pytest
import tempfile
import os
import sys
import yaml
from unittest.mock import Mock, patch, MagicMock

# Add the test directory to path for importing test utilities
sys.path.append(os.path.dirname(__file__))
from test_data_utils import get_test_pdb_path


class TestIntegration:
    """Integration tests for EasyMD workflow"""
    
    @pytest.fixture
    def sample_pdb_content(self):
        """Sample PDB file content for testing"""
        return """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
END
"""
    
    @pytest.fixture
    def temp_pdb_file(self, sample_pdb_content):
        """Create a temporary PDB file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(sample_pdb_content)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)
    
    @pytest.fixture
    def temp_config_file(self):
        """Create a temporary config file"""
        config_data = {
            'protein': get_test_pdb_path(),
            'steps': 100,
            'temperature': 300,
            'solvate': True,
            'interval': 10
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
            yaml.dump(config_data, f)
            temp_file = f.name
        yield temp_file
        os.unlink(temp_file)
    
    @patch('EasyMD.sysGenerator.SysGenerator.PDBFixer')
    @patch('EasyMD.sysGenerator.SysGenerator.SystemGenerator')
    @patch('EasyMD.sysGenerator.SysGenerator.Modeller')
    @patch('os.mkdir')
    @patch('os.path.isdir')
    def test_full_workflow_protein_only(self, mock_isdir, mock_mkdir, mock_modeller_class, 
                                       mock_sys_gen, mock_pdb_fixer, temp_pdb_file):
        """Test complete workflow for protein-only simulation"""
        from EasyMD.argManager.manager import ArgManager
        from EasyMD.sysGenerator.SysGenerator import SysGenerator
        from EasyMD.simRunner.simRunner import SimRunner
        import argparse
        
        # Setup mocks
        mock_isdir.return_value = False
        mock_fixer_instance = Mock()
        mock_pdb_fixer.return_value = mock_fixer_instance
        mock_modeller_instance = Mock()
        mock_modeller_class.return_value = mock_modeller_instance
        mock_system_generator = Mock()
        mock_sys_gen.return_value = mock_system_generator
        mock_system = Mock()
        mock_system_generator.create_system.return_value = mock_system
        
        # Test argument parsing
        parser = argparse.ArgumentParser()
        with patch('sys.argv', ['test', '--protein', temp_pdb_file, '--steps', '100', '--solvate','--interval','10']):
            arg_manager = ArgManager(parser)
            config = arg_manager.get_args()
        
        # Test system generation
        with patch('EasyMD.sysGenerator.SysGenerator.SysGenerator._prep_prot') as mock_prep:
            mock_prep.return_value = (mock_modeller_instance, mock_system)
            sys_gen = SysGenerator(config)
        
        # Test simulation runner setup
        with patch('EasyMD.simRunner.simRunner.SolvatedRunner') as mock_runner:
            mock_runner_instance = Mock()
            mock_runner.return_value = mock_runner_instance
            sim_runner = SimRunner(config, sys_gen.modeller, sys_gen.system)
        
        # Verify the workflow
        assert config.protein == temp_pdb_file
        assert config.steps == 100
        assert config.solvate is True
        assert sys_gen.modeller == mock_modeller_instance
        assert sys_gen.system == mock_system
        assert sim_runner.sim == mock_runner_instance
    
    def test_config_file_integration(self, temp_config_file):
        """Test integration with config file loading"""
        from EasyMD.argManager.manager import ArgManager
        import argparse
        
        parser = argparse.ArgumentParser()
        with patch('sys.argv', ['test', '--config', temp_config_file]):
            arg_manager = ArgManager(parser)
            config = arg_manager.get_args()
        
        assert config.protein == get_test_pdb_path()
        assert config.steps == 100
        assert config.temperature == 300
        assert config.solvate is True
        assert config.interval == 10
    
    @pytest.mark.slow
    def test_error_handling_missing_protein(self):
        """Test error handling when protein file is missing"""
        from EasyMD.argManager.manager import ArgManager
        import argparse
        
        parser = argparse.ArgumentParser()
        with patch('sys.argv', ['test', '--protein', 'nonexistent.pdb', '--steps', '100', '--solvate']):
            with pytest.raises(SystemExit) as e:
                arg_manager = ArgManager(parser)
                config = arg_manager.get_args()
            assert e.type == SystemExit
            assert e.value.code == 1
            
            # This should handle the missing file gracefully in the system generator
        

            
           
    
    def test_restart_workflow_integration(self):
        """Test restart workflow integration"""
        from EasyMD.argManager.manager import ArgManager
        from EasyMD.sysGenerator.SysGenerator import SysGenerator
        import argparse
        
        # Create mock restart directory structure with all required files
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create restart_setup.yml
            restart_setup = {
                'protein_force_field': 'amber14-all.xml',
                'solvate': True,
                'water_force_field': 'tip3p.xml',
                'last_state': 0,
                'temperature': 300,
                'step_size': 0.002,
                'friction_coeff': 1.0,
                'reporting_interval': 1000
            }
            
            setup_file = os.path.join(temp_dir, 'restart_setup.yml')
            with open(setup_file, 'w') as f:
                yaml.dump(restart_setup, f)
            
            # Create restart_model.pdb
            restart_model = os.path.join(temp_dir, 'restart_model.pdb')
            restart_pdb_content = """HEADER    RESTART MODEL
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
END
"""
            with open(restart_model, 'w') as f:
                f.write(restart_pdb_content)
            
            # Create last_state.xml (mock OpenMM state file)
            last_state_file = os.path.join(temp_dir, 'last_state.xml')
            last_state_content = """<?xml version="1.0" ?>
<State version="1" openmmVersion="8.1.1">
  <Parameters>
    <Parameter name="t" value="0.0"/>
  </Parameters>
  <Positions>
    <Position x="2.0154" y="1.6967" z="1.4365"/>
    <Position x="1.9030" y="1.6101" z="1.4618"/>
    <Position x="1.7664" y="1.6849" z="1.4897"/>
    <Position x="1.7764" y="1.8067" z="1.5086"/>
  </Positions>
  <Velocities>
    <Velocity x="0.0" y="0.0" z="0.0"/>
    <Velocity x="0.0" y="0.0" z="0.0"/>
    <Velocity x="0.0" y="0.0" z="0.0"/>
    <Velocity x="0.0" y="0.0" z="0.0"/>
  </Velocities>
</State>
"""
            with open(last_state_file, 'w') as f:
                f.write(last_state_content)
            
            # Test argument parsing with complete restart directory
            parser = argparse.ArgumentParser()
            with patch('sys.argv', ['test', '--restart', temp_dir]):
                arg_manager = ArgManager(parser)
                config = arg_manager.get_args()
            
            assert config.restart == temp_dir
            
            # Test system generator with restart
            with patch('EasyMD.sysGenerator.SysGenerator.PDBFile'), \
                 patch('EasyMD.sysGenerator.SysGenerator.SystemGenerator'), \
                 patch('EasyMD.sysGenerator.SysGenerator.Modeller'):
                sys_gen = SysGenerator(config)
                assert sys_gen.config.restart == temp_dir