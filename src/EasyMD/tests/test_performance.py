import pytest
import time
from unittest.mock import Mock, patch
import tempfile
import os


class TestPerformance:
    """Performance and benchmarking tests"""
    
    @pytest.mark.slow
    def test_arg_manager_config_parsing_performance(self):
        """Test that config file parsing is reasonably fast"""
        import yaml
        from EasyMD.argManager.manager import ArgManager
        import argparse
        
        # Create a large config file
        large_config = {f'param_{i}': f'value_{i}' for i in range(1000)}
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
            yaml.dump(large_config, f)
            config_file = f.name
        
        try:
            parser = argparse.ArgumentParser()
            start_time = time.time()
            
            with patch('sys.argv', ['test', '--config', config_file, '--protein', 'test.pdb', '--steps', '100', '--solvate']):
                ArgManager(parser)
            
            end_time = time.time()
            parsing_time = end_time - start_time
            
            # Should parse large config in under 1 second
            assert parsing_time < 1.0, f"Config parsing took {parsing_time:.2f}s, expected < 1.0s"
        finally:
            os.unlink(config_file)
    
    @pytest.mark.slow
    def test_system_generator_memory_usage(self):
        """Test that SysGenerator doesn't leak memory during initialization"""
        import psutil
        import gc
        from EasyMD.sysGenerator.sysGenerator import SysGenerator
        
        mock_config = Mock()
        mock_config.restart = None
        mock_config.outdir = None
        mock_config.protein = 'test.pdb'
        mock_config.ligand = None
        
        # Get initial memory usage
        process = psutil.Process()
        initial_memory = process.memory_info().rss
        
        # Create multiple instances
        instances = []
        with patch('EasyMD.sysGenerator.sysGenerator.SysGenerator._setup') as mock_setup:
            mock_setup.return_value = (Mock(), Mock())
            
            for _ in range(10):
                instances.append(SysGenerator(mock_config))
        
        # Force garbage collection
        del instances
        gc.collect()
        
        final_memory = process.memory_info().rss
        memory_increase = final_memory - initial_memory
        
        # Memory increase should be reasonable (less than 100MB)
        assert memory_increase < 100 * 1024 * 1024, f"Memory increased by {memory_increase / 1024 / 1024:.1f}MB"
    
    def test_concurrent_arg_parsing(self):
        """Test that argument parsing works with concurrent access"""
        import threading
        import argparse
        from EasyMD.argManager.manager import ArgManager
        
        results = []
        errors = []
        
        def parse_args():
            try:
                parser = argparse.ArgumentParser()
                with patch('sys.argv', ['test', '--protein', 'test.pdb', '--steps', '100', '--solvate']):
                    manager = ArgManager(parser)
                    results.append(manager.get_args())
            except Exception as e:
                errors.append(e)
        
        # Run multiple threads
        threads = [threading.Thread(target=parse_args) for _ in range(5)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        
        assert len(errors) == 0, f"Errors occurred: {errors}"
        assert len(results) == 5, f"Expected 5 results, got {len(results)}"