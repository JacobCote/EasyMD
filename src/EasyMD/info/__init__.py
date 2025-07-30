"""
EasyMD Info Module

This module provides detailed structural information about PDB files,
including chain analysis, missing residues, ligands, and other structural features.
"""

from .infoRunner import InfoRunner
from .infoManager import InfoManager

__all__ = ['InfoRunner', 'InfoManager']