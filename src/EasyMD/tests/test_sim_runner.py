import pytest
from unittest.mock import Mock, patch, MagicMock
from EasyMD.simRunner.simRunner import SimRunner


class TestSimRunner:
    """Test suite for SimRunner class"""
    
    @pytest.fixture
    def mock_config(self):
        """Create a mock configuration object"""
        config = Mock()
        config.restart = None
        config.simulated_annealing = False
        config.GBIS = False
        return config
    
    @pytest.fixture
    def mock_modeller(self):
        """Create a mock modeller object"""
        return Mock()
    
    @pytest.fixture
    def mock_system(self):
        """Create a mock system object"""
        return Mock()
    
    @pytest.fixture
    def mock_setup(self):
        """Create a mock setup object"""
        return Mock()
    
    def test_init_basic(self, mock_config, mock_modeller, mock_system):
        """Test basic SimRunner initialization"""
        with patch('EasyMD.simRunner.simRunner.SolvatedRunner') as mock_solvated:
            mock_solvated.return_value = Mock()
            
            runner = SimRunner(mock_config, mock_modeller, mock_system)
            
            assert runner.config == mock_config
            assert runner.modeller == mock_modeller
            assert runner.system == mock_system
            assert runner.setup is None
            assert runner.sim is not None
    
    def test_setup_sim_solvated_runner(self, mock_config, mock_modeller, mock_system):
        """Test _setupSim creates SolvatedRunner for default case"""
        with patch('EasyMD.simRunner.simRunner.SolvatedRunner') as mock_solvated:
            mock_instance = Mock()
            mock_solvated.return_value = mock_instance
            
            runner = SimRunner(mock_config, mock_modeller, mock_system)
            
            mock_solvated.assert_called_once_with(
                mock_config, 
                modeller=mock_modeller, 
                system=mock_system
            )
            assert runner.sim == mock_instance
    
    def test_setup_sim_gbis_runner(self, mock_config, mock_modeller, mock_system):
        """Test _setupSim creates GBISRunner when GBIS is True"""
        mock_config.GBIS = True
        
        with patch('EasyMD.simRunner.simRunner.GBISRunner') as mock_gbis:
            mock_instance = Mock()
            mock_gbis.return_value = mock_instance
            
            runner = SimRunner(mock_config, mock_modeller, mock_system)
            
            mock_gbis.assert_called_once_with(
                mock_config, 
                modeller=mock_modeller, 
                system=mock_system
            )
            assert runner.sim == mock_instance
    
    def test_setup_sim_annealing_runner(self, mock_config, mock_modeller, mock_system):
        """Test _setupSim creates AnnealingRunner when simulated_annealing is True"""
        mock_config.simulated_annealing = True
        
        with patch('EasyMD.simRunner.simRunner.AnnealingRunner') as mock_annealing:
            mock_instance = Mock()
            mock_annealing.return_value = mock_instance
            
            runner = SimRunner(mock_config, mock_modeller, mock_system)
            
            mock_annealing.assert_called_once_with(
                mock_config, 
                modeller=mock_modeller, 
                system=mock_system
            )
            assert runner.sim == mock_instance
    
    def test_setup_sim_restart_runner(self, mock_config, mock_modeller, mock_system, mock_setup):
        """Test _setupSim creates Restarter when restart is not None"""
        mock_config.restart = "some_restart_dir"
        
        with patch('EasyMD.simRunner.simRunner.Restarter') as mock_restarter:
            mock_instance = Mock()
            mock_restarter.return_value = mock_instance
            
            runner = SimRunner(mock_config, mock_modeller, mock_system, mock_setup)
            
            mock_restarter.assert_called_once_with(
                mock_setup,
                mock_config, 
                modeller=mock_modeller, 
                system=mock_system
            )
            assert runner.sim == mock_instance
    
    def test_run_calls_sim_run(self, mock_config, mock_modeller, mock_system):
        """Test that run() method calls the sim's run() method"""
        with patch('EasyMD.simRunner.simRunner.SolvatedRunner') as mock_solvated:
            mock_sim = Mock()
            mock_solvated.return_value = mock_sim
            
            runner = SimRunner(mock_config, mock_modeller, mock_system)
            runner.run()
            
            mock_sim.run.assert_called_once()
    
    def test_priority_restart_over_annealing(self, mock_config, mock_modeller, mock_system, mock_setup):
        """Test that restart takes priority over simulated annealing"""
        mock_config.restart = "restart_dir"
        mock_config.simulated_annealing = True
        
        with patch('EasyMD.simRunner.simRunner.Restarter') as mock_restarter:
            mock_instance = Mock()
            mock_restarter.return_value = mock_instance
            
            runner = SimRunner(mock_config, mock_modeller, mock_system, mock_setup)
            
            mock_restarter.assert_called_once()
            assert runner.sim == mock_instance
    
    def test_priority_annealing_over_gbis(self, mock_config, mock_modeller, mock_system):
        """Test that simulated annealing takes priority over GBIS"""
        mock_config.simulated_annealing = True
        mock_config.GBIS = True
        
        with patch('EasyMD.simRunner.simRunner.AnnealingRunner') as mock_annealing:
            mock_instance = Mock()
            mock_annealing.return_value = mock_instance
            
            runner = SimRunner(mock_config, mock_modeller, mock_system)
            
            mock_annealing.assert_called_once()
            assert runner.sim == mock_instance
    
    def test_priority_gbis_over_solvated(self, mock_config, mock_modeller, mock_system):
        """Test that GBIS takes priority over solvated runner"""
        mock_config.GBIS = True
        
        with patch('EasyMD.simRunner.simRunner.GBISRunner') as mock_gbis:
            mock_instance = Mock()
            mock_gbis.return_value = mock_instance
            
            runner = SimRunner(mock_config, mock_modeller, mock_system)
            
            mock_gbis.assert_called_once()
            assert runner.sim == mock_instance