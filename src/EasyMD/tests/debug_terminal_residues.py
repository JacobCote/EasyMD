#!/usr/bin/env python3
"""
Debug script to test terminal residue limiting functionality.
This will help identify why --max-terminal-residues isn't working.
"""

import sys
import argparse
import tempfile
import os
from unittest.mock import Mock, patch
from src.EasyMD.argManager.manager import ArgManager
from src.EasyMD.sysGenerator.sysGenerator import SysGenerator

def create_test_pdb():
    """Create a test PDB file with missing residues."""
    pdb_content = """HEADER    TEST PROTEIN WITH MISSING RESIDUES
REMARK   Missing residues: 1-5, 150-155 (simulated)
ATOM      1  N   ALA A   6      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   6      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   6      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   6      17.764  18.067  15.086  1.00 20.00           O  
ATOM      5  N   VAL A  149     18.154  17.967  15.365  1.00 20.00           N  
ATOM      6  CA  VAL A  149     17.030  17.101  15.618  1.00 20.00           C  
END
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_content)
        return f.name

def test_argument_parsing():
    """Test that arguments are parsed correctly."""
    print("1. Testing argument parsing...")
    
    test_pdb = create_test_pdb()
    
    try:
        parser = argparse.ArgumentParser()
        with patch('sys.argv', ['test', '--protein', test_pdb, '--steps', '1000', '--solvate', 
                                '--missing-residues', 'terminal-only', '--max-terminal-residues', '3']):
            manager = ArgManager(parser)
            args = manager.get_args()
            
            print(f"  ✅ missing_residues: {args.missing_residues}")
            print(f"  ✅ max_terminal_residues: {args.max_terminal_residues}")
            
            assert args.missing_residues == 'terminal-only'
            assert args.max_terminal_residues == 3
            
            return args
    finally:
        os.unlink(test_pdb)

def test_sys_generator_config():
    """Test that SysGenerator receives the config correctly."""
    print("\n2. Testing SysGenerator config...")
    
    # Create mock config with the exact values
    mock_config = Mock()
    mock_config.missing_residues = 'terminal-only'
    mock_config.max_terminal_residues = 3
    mock_config.terminal_residue_types = ['ACE', 'NME', 'NH2', 'COOH']
    mock_config.skip_missing_loops = False
    mock_config.conservative_missing = False
    
    # Test the getattr calls in _handle_missing_residues
    sys_gen = SysGenerator.__new__(SysGenerator)
    sys_gen.config = mock_config
    
    missing_residues_strategy = getattr(sys_gen.config, 'missing_residues', 'auto')
    max_terminal_residues = getattr(sys_gen.config, 'max_terminal_residues', 5)
    
    print(f"  ✅ Strategy from config: {missing_residues_strategy}")
    print(f"  ✅ Max terminal from config: {max_terminal_residues}")
    
    assert missing_residues_strategy == 'terminal-only'
    assert max_terminal_residues == 3
    
    return sys_gen

def test_filtering_logic():
    """Test the filtering logic directly."""
    print("\n3. Testing filtering logic...")
    
    # Create SysGenerator with test config
    sys_gen = test_sys_generator_config()
    
    # Create mock fixer with many missing residues
    mock_fixer = Mock()
    mock_fixer.missingResidues = {
        'A': [f'residue_{i}' for i in range(1, 8)]  # 7 residues
    }
    
    print(f"  Before filtering: Chain A has {len(mock_fixer.missingResidues['A'])} residues")
    print(f"    Residues: {mock_fixer.missingResidues['A']}")
    
    # Apply the filtering directly
    sys_gen._filter_terminal_residues(mock_fixer, max_terminal_residues=3, allowed_types=['ACE', 'NME'])
    
    print(f"  After filtering: Chain A has {len(mock_fixer.missingResidues['A'])} residues")
    print(f"    Residues: {mock_fixer.missingResidues['A']}")
    
    # Verify the result
    assert len(mock_fixer.missingResidues['A']) == 3, f"Expected 3 residues, got {len(mock_fixer.missingResidues['A'])}"
    
    return mock_fixer

def test_full_handle_missing_residues():
    """Test the complete _handle_missing_residues method."""
    print("\n4. Testing complete _handle_missing_residues method...")
    
    # Create SysGenerator with test config
    sys_gen = test_sys_generator_config()
    
    # Create mock fixer with many missing residues
    mock_fixer = Mock()
    mock_fixer.missingResidues = {
        'A': [f'residue_{i}' for i in range(1, 8)],  # 7 residues
        'B': [f'residue_{i}' for i in range(1, 3)]   # 2 residues
    }
    
    print(f"  Before _handle_missing_residues:")
    for chain_id, residues in mock_fixer.missingResidues.items():
        print(f"    Chain {chain_id}: {len(residues)} residues")
    
    # Call the complete method
    sys_gen._handle_missing_residues(mock_fixer)
    
    print(f"  After _handle_missing_residues:")
    for chain_id, residues in mock_fixer.missingResidues.items():
        print(f"    Chain {chain_id}: {len(residues)} residues")
    
    # Verify the results
    assert len(mock_fixer.missingResidues['A']) == 3, f"Chain A: Expected 3 residues, got {len(mock_fixer.missingResidues['A'])}"
    assert len(mock_fixer.missingResidues['B']) == 2, f"Chain B: Expected 2 residues, got {len(mock_fixer.missingResidues['B'])}"
    
    return mock_fixer

def test_with_real_pdbfixer():
    """Test with a real PDBFixer instance to see what happens."""
    print("\n5. Testing with real PDBFixer...")
    
    test_pdb = create_test_pdb()
    
    try:
        from pdbfixer import PDBFixer
        
        # Create real PDBFixer instance
        fixer = PDBFixer(filename=test_pdb)
        fixer.findMissingResidues()
        
        print(f"  Real PDBFixer found missing residues:")
        for chain_id, residues in fixer.missingResidues.items():
            print(f"    Chain {chain_id}: {len(residues)} residues")
            for residue in residues[:5]:  # Show first 5
                print(f"      {residue}")
        
        if not fixer.missingResidues:
            print("  No missing residues found by PDBFixer - this might be why filtering isn't working!")
            return None
        
        # Create SysGenerator and test filtering
        sys_gen = test_sys_generator_config()
        
        # Store original for comparison
        original_count = sum(len(residues) for residues in fixer.missingResidues.values())
        
        # Apply filtering
        sys_gen._handle_missing_residues(fixer)
        
        # Check result
        filtered_count = sum(len(residues) for residues in fixer.missingResidues.values())
        
        print(f"  Original total missing residues: {original_count}")
        print(f"  Filtered total missing residues: {filtered_count}")
        
        return fixer
        
    except ImportError:
        print("  ⚠️  PDBFixer not available for testing")
        return None
    finally:
        os.unlink(test_pdb)

def main():
    """Run all debug tests."""
    print("🔍 Debugging Terminal Residue Limiting Issue")
    print("=" * 60)
    
    try:
        # Test 1: Argument parsing
        args = test_argument_parsing()
        
        # Test 2: SysGenerator config
        sys_gen = test_sys_generator_config()
        
        # Test 3: Filtering logic
        filtered_fixer = test_filtering_logic()
        
        # Test 4: Complete method
        complete_fixer = test_full_handle_missing_residues()
        
        # Test 5: Real PDBFixer
        real_fixer = test_with_real_pdbfixer()
        
        print("\n" + "=" * 60)
        print("🎯 DIAGNOSIS:")
        print("=" * 60)
        
        if real_fixer is None:
            print("❌ ISSUE FOUND: PDBFixer might not be finding missing residues in your structure")
            print("   This could be because:")
            print("   1. The PDB file doesn't have missing residues")
            print("   2. PDBFixer can't detect missing residues in this structure")
            print("   3. The structure is already complete")
        elif not real_fixer.missingResidues:
            print("❌ ISSUE FOUND: No missing residues detected by PDBFixer")
            print("   The filtering only works if PDBFixer finds missing residues first")
        else:
            print("✅ All tests passed - the filtering logic is working correctly")
            print("   If you're still seeing all residues added, the issue might be:")
            print("   1. Different arguments being passed than expected")
            print("   2. Multiple PDBFixer instances being created")
            print("   3. The filtering being overridden somewhere else")
        
        print("\n💡 RECOMMENDATIONS:")
        print("1. Check that your PDB file actually has missing residues")
        print("2. Look at the console output for 'Chain X: keeping Y terminal missing residues'")
        print("3. Verify you're using --missing-residues terminal-only")
        print("4. Try with a structure known to have missing residues")
        
    except Exception as e:
        print(f"\n❌ DEBUG FAILED: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()