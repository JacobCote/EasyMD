import pytest
from unittest.mock import Mock, patch, mock_open
import tempfile
import os


class TestUtils:
    """Test suite for utility functions"""
    
    def test_format_index(self):
        """Test _formatIndex utility function"""
        from EasyMD.utils.utils import _formatIndex
        
        assert _formatIndex(5, 5) == "    5"
        assert _formatIndex(123, 5) == "  123"
        assert _formatIndex(99999, 5) == "99999"
    
    def test_write_footer(self):
        """Test writeFooter utility function"""
        from EasyMD.utils.utils import writeFooter
        
        # Create mock topology with bonds
        mock_topology = Mock()
        mock_topology.bonds.return_value = []  # No bonds for simplicity
        
        # Create mock file
        mock_file = Mock()
        
        # Test that writeFooter can be called without errors
        writeFooter(mock_topology, mock_file)
        
        # Should call bonds() on topology
        mock_topology.bonds.assert_called_once()
    
    @patch('EasyMD.utils.utils.PDBFile.writeFile')
    @patch('EasyMD.utils.utils.writeFooter')
    def test_pdb_write_all(self, mock_write_footer, mock_write_file):
        """Test PDBwrite_all utility function"""
        from EasyMD.utils.utils import PDBwrite_all
        
        mock_modeller = Mock()
        mock_modeller.topology = Mock()
        mock_modeller.positions = Mock()
        
        with patch('builtins.open', mock_open()) as mock_file:
            PDBwrite_all(mock_modeller, "test.pdb")
            
            # Should open file twice (once for writing, once for appending footer)
            assert mock_file.call_count >= 1
            mock_write_file.assert_called_once()
            mock_write_footer.assert_called_once()
    
    def test_delete_pcap(self):
        """Test deletePcap utility function"""
        from EasyMD.utils.utils import deletePcap
        
        mock_modeller = Mock()
        mock_topology = Mock()
        mock_modeller.topology = mock_topology
        
        # Create mock chain
        mock_chain = Mock()
        
        # Mock DNA residue with atoms to delete
        mock_dna_residue = Mock()
        mock_dna_residue.name = "DC"  # DNA cytosine
        mock_dna_residue.chain = mock_chain
        
        # Mock atom with phosphate group
        mock_atom = Mock()
        mock_atom.name = "P"  # Phosphorus atom
        mock_dna_residue.atoms.return_value = [mock_atom]
        
        # Mock regular residue
        mock_regular_residue = Mock()
        mock_regular_residue.name = "ALA"
        mock_regular_residue.chain = Mock()  # Different chain
        
        # Set up topology to return residues twice (as the function iterates twice)
        mock_topology.residues.return_value = [mock_dna_residue, mock_regular_residue]
        
        result = deletePcap(mock_modeller)
        
        # Should call delete with phosphate atoms from DNA residues
        mock_modeller.delete.assert_called_once()
        assert result == mock_modeller