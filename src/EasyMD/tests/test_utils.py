import pytest
from unittest.mock import Mock, patch, mock_open
import tempfile
import os


class TestUtils:
    """Test suite for utility functions"""
    
    def test_format_index(self):
        """Test _formatIndex utility function"""
        from EasyMD.utils.utils import _formatIndex
        
        assert _formatIndex(5) == "00005"
        assert _formatIndex(123) == "00123"
        assert _formatIndex(99999) == "99999"
    
    @patch('builtins.open', new_callable=mock_open)
    def test_write_footer(self, mock_file):
        """Test writeFooter utility function"""
        from EasyMD.utils.utils import writeFooter
        
        writeFooter("test_file.txt", "Test content")
        
        mock_file.assert_called_once_with("test_file.txt", "a")
        mock_file().write.assert_called_once_with("Test content")
    
    @patch('EasyMD.utils.utils.PDBFile.writeFile')
    def test_pdb_write_all(self, mock_write_file):
        """Test PDBwrite_all utility function"""
        from EasyMD.utils.utils import PDBwrite_all
        
        mock_modeller = Mock()
        mock_modeller.topology = Mock()
        mock_modeller.positions = Mock()
        
        with patch('builtins.open', mock_open()) as mock_file:
            PDBwrite_all(mock_modeller, "test.pdb")
            
            mock_file.assert_called_once_with("test.pdb", 'w')
            mock_write_file.assert_called_once()
    
    def test_delete_pcap(self):
        """Test deletePcap utility function"""
        from EasyMD.utils.utils import deletePcap
        
        mock_modeller = Mock()
        mock_topology = Mock()
        mock_modeller.topology = mock_topology
        
        # Mock residues with PCAP
        mock_residue1 = Mock()
        mock_residue1.name = "PCAP"
        mock_residue2 = Mock()
        mock_residue2.name = "ALA"
        
        mock_topology.residues.return_value = [mock_residue1, mock_residue2]
        
        result = deletePcap(mock_modeller)
        
        # Should call delete with PCAP residues
        mock_modeller.delete.assert_called_once_with([mock_residue1])
        assert result == mock_modeller