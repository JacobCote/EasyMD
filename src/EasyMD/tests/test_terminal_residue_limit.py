#!/usr/bin/env python3
"""
Test script to verify that --max-terminal-residues parameter works correctly.
This script demonstrates the terminal residue limiting functionality.
"""

import tempfile
import os
from unittest.mock import Mock
from EasyMD.sysGenerator.SysGenerator import SysGenerator

def test_terminal_residue_filtering():
    """Test that terminal residue filtering respects the max limit."""
    
    print("Testing terminal residue filtering...")
    
    # Create a mock config with terminal-only strategy and limit of 3
    mock_config = Mock()
    mock_config.missing_residues = 'terminal-only'
    mock_config.max_terminal_residues = 3
    mock_config.terminal_residue_types = ['ACE', 'NME', 'NH2', 'COOH']
    mock_config.skip_missing_loops = False
    mock_config.conservative_missing = False
    
    # Create a mock fixer with many missing residues
    mock_fixer = Mock()
    mock_fixer.missingResidues = {
        'A': ['res1', 'res2', 'res3', 'res4', 'res5', 'res6', 'res7'],  # 7 residues
        'B': ['res1', 'res2']  # 2 residues
    }
    
    # Create SysGenerator instance (without calling __init__)
    sys_gen = SysGenerator.__new__(SysGenerator)
    sys_gen.config = mock_config
    
    print(f"Before filtering:")
    print(f"  Chain A: {len(mock_fixer.missingResidues['A'])} residues")
    print(f"  Chain B: {len(mock_fixer.missingResidues['B'])} residues")
    
    # Apply terminal filtering with max 3 residues
    sys_gen._filter_terminal_residues(mock_fixer, max_terminal_residues=3, allowed_types=['ACE', 'NME'])
    
    print(f"\nAfter filtering (max 3 terminal residues):")
    for chain_id, residues in mock_fixer.missingResidues.items():
        print(f"  Chain {chain_id}: {len(residues)} residues - {residues}")
    
    # Verify the results
    assert 'A' in mock_fixer.missingResidues
    assert 'B' in mock_fixer.missingResidues
    
    # Chain A should have exactly 3 residues (limited from 7)
    assert len(mock_fixer.missingResidues['A']) == 3, f"Expected 3 residues for chain A, got {len(mock_fixer.missingResidues['A'])}"
    
    # Chain B should have 2 residues (all of them, since 2 < 3)
    assert len(mock_fixer.missingResidues['B']) == 2, f"Expected 2 residues for chain B, got {len(mock_fixer.missingResidues['B'])}"
    
    print("\n✅ Terminal residue filtering test PASSED!")
    print("The --max-terminal-residues parameter is working correctly.")

def test_different_limits():
    """Test different max terminal residue limits."""
    
    print("\nTesting different terminal residue limits...")
    
    test_cases = [
        (1, "Very restrictive"),
        (2, "Restrictive"), 
        (5, "Moderate"),
        (10, "Permissive")
    ]
    
    for max_limit, description in test_cases:
        print(f"\n{description} limit (max {max_limit} residues):")
        
        # Create mock config
        mock_config = Mock()
        mock_config.missing_residues = 'terminal-only'
        mock_config.max_terminal_residues = max_limit
        mock_config.terminal_residue_types = ['ACE', 'NME', 'NH2', 'COOH']
        mock_config.skip_missing_loops = False
        mock_config.conservative_missing = False
        
        # Create mock fixer with 8 missing residues
        mock_fixer = Mock()
        mock_fixer.missingResidues = {
            'A': [f'res{i}' for i in range(1, 9)]  # 8 residues
        }
        
        # Create SysGenerator instance
        sys_gen = SysGenerator.__new__(SysGenerator)
        sys_gen.config = mock_config
        
        # Apply filtering
        sys_gen._filter_terminal_residues(mock_fixer, max_terminal_residues=max_limit, allowed_types=['ACE', 'NME'])
        
        result_count = len(mock_fixer.missingResidues['A'])
        expected_count = min(max_limit, 8)  # Should not exceed the limit or original count
        
        print(f"  Original: 8 residues → Filtered: {result_count} residues")
        assert result_count == expected_count, f"Expected {expected_count}, got {result_count}"
    
    print("\n✅ Different limits test PASSED!")

def test_conservative_filtering():
    """Test conservative filtering with max terminal residues."""
    
    print("\nTesting conservative filtering...")
    
    # Create mock config with conservative filtering
    mock_config = Mock()
    mock_config.missing_residues = 'auto'
    mock_config.max_terminal_residues = 3
    mock_config.terminal_residue_types = ['ACE', 'NME', 'NH2', 'COOH']
    mock_config.skip_missing_loops = False
    mock_config.conservative_missing = True
    
    # Create mock fixer with different chain lengths
    mock_fixer = Mock()
    mock_fixer.missingResidues = {
        'A': ['res1', 'res2'],  # 2 residues - should keep
        'B': ['res1', 'res2', 'res3', 'res4', 'res5'],  # 5 residues - should skip (exceeds limit of 3)
        'C': ['res1', 'res2', 'res3']  # 3 residues - should keep (exactly at limit)
    }
    
    # Create SysGenerator instance
    sys_gen = SysGenerator.__new__(SysGenerator)
    sys_gen.config = mock_config
    
    print("Before conservative filtering:")
    for chain_id, residues in mock_fixer.missingResidues.items():
        print(f"  Chain {chain_id}: {len(residues)} residues")
    
    # Apply conservative filtering
    sys_gen._apply_conservative_filtering(mock_fixer, max_terminal_residues=3)
    
    print("\nAfter conservative filtering (max 3 residues per chain):")
    for chain_id, residues in mock_fixer.missingResidues.items():
        print(f"  Chain {chain_id}: {len(residues)} residues")
    
    # Verify results
    assert 'A' in mock_fixer.missingResidues  # Should keep (2 <= 3)
    assert 'B' not in mock_fixer.missingResidues  # Should remove (5 > 3)
    assert 'C' in mock_fixer.missingResidues  # Should keep (3 <= 3)
    
    print("\n✅ Conservative filtering test PASSED!")

if __name__ == "__main__":
    print("Testing EasyMD Terminal Residue Limiting Functionality")
    print("=" * 60)
    
    try:
        test_terminal_residue_filtering()
        test_different_limits()
        test_conservative_filtering()
        
        print("\n" + "=" * 60)
        print("🎉 ALL TESTS PASSED!")
        print("The --max-terminal-residues parameter is working correctly.")
        print("\nYou can now use commands like:")
        print("  python -m EasyMD --protein structure.pdb --steps 10000 --solvate \\")
        print("    --missing-residues terminal-only --max-terminal-residues 3")
        print("\nThis will limit terminal missing residues to 3 per chain.")
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()