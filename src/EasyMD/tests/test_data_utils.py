"""
Utility functions for accessing test data files.
"""
import os
from pathlib import Path

def get_test_data_path(filename):
    """
    Get the full path to a test data file.
    
    Args:
        filename (str): Name of the test data file
        
    Returns:
        str: Full path to the test data file
    """
    test_dir = Path(__file__).parent
    data_dir = test_dir / "data"
    return str(data_dir / filename)

def get_test_pdb_path():
    """Get path to the main test PDB file."""
    return get_test_data_path("test.pdb")

def get_4zgm_pdb_path():
    """Get path to the 4zgm test PDB file."""
    return get_test_data_path("4zgm.pdb")

def get_3s8l_pdb_path():
    """Get path to the 3s8l test PDB file."""
    return get_test_data_path("3s8l.pdb")

def get_missing_structure_pdb_path():
    """Get path to the test structure with missing residues."""
    return get_test_data_path("test_structure_with_missing.pdb")