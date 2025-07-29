#!/usr/bin/env python3
"""
Test script to check if your PDB structure has missing residues that PDBFixer can detect.
Run this with your actual PDB file to see what's happening.
"""

import sys
from pdbfixer import PDBFixer

def analyze_pdb_file(pdb_file):
    """Test a PDB file to see what missing residues PDBFixer finds."""
    
    print(f"Testing PDB file: {pdb_file}")
    print("=" * 60)
    
    try:
        # Create PDBFixer instance
        fixer = PDBFixer(filename=pdb_file)
        
        # Find missing residues
        print("1. Finding missing residues...")
        fixer.findMissingResidues()
        
        if fixer.missingResidues:
            print(f"✅ Found missing residues:")
            total_missing = 0
            for chain_id, residues in fixer.missingResidues.items():
                print(f"  Chain {chain_id}: {len(residues)} missing residues")
                total_missing += len(residues)
                for i, residue in enumerate(residues):
                    if i < 5:  # Show first 5
                        print(f"    {residue}")
                    elif i == 5:
                        print(f"    ... and {len(residues) - 5} more")
                        break
            
            print(f"\nTotal missing residues: {total_missing}")
            
            # Test what happens with max_terminal_residues = 3
            print(f"\n2. Testing with --max-terminal-residues 3:")
            
            # Simulate the filtering
            from src.EasyMD.sysGenerator.sysGenerator import SysGenerator
            from unittest.mock import Mock
            
            # Create mock config
            mock_config = Mock()
            mock_config.missing_residues = 'terminal-only'
            mock_config.max_terminal_residues = 3
            mock_config.terminal_residue_types = ['ACE', 'NME', 'NH2', 'COOH']
            mock_config.skip_missing_loops = False
            mock_config.conservative_missing = False
            
            # Create SysGenerator instance
            sys_gen = SysGenerator.__new__(SysGenerator)
            sys_gen.config = mock_config
            
            # Store original counts
            original_counts = {chain_id: len(residues) for chain_id, residues in fixer.missingResidues.items()}
            
            # Apply filtering
            sys_gen._handle_missing_residues(fixer)
            
            # Show results
            print(f"Results after filtering:")
            for chain_id, residues in fixer.missingResidues.items():
                original_count = original_counts.get(chain_id, 0)
                print(f"  Chain {chain_id}: {original_count} → {len(residues)} residues")
                if len(residues) > 3:
                    print(f"    ❌ ERROR: Still has {len(residues)} residues (should be ≤ 3)")
                else:
                    print(f"    ✅ OK: {len(residues)} residues (within limit)")
            
        else:
            print("❌ No missing residues found by PDBFixer")
            print("\nThis means:")
            print("1. Your structure is already complete, OR")
            print("2. PDBFixer can't detect missing residues in this format, OR") 
            print("3. The PDB file doesn't have SEQRES records or proper sequence info")
            print("\nSince no missing residues are found, the --max-terminal-residues")
            print("parameter has no effect (there's nothing to limit).")
        
        # Also check missing atoms
        print(f"\n3. Checking missing atoms...")
        fixer.findMissingAtoms()
        if fixer.missingAtoms:
            print(f"✅ Found missing atoms in {len(fixer.missingAtoms)} residues")
        else:
            print("❌ No missing atoms found")
        
        # Check non-standard residues
        print(f"\n4. Checking non-standard residues...")
        fixer.findNonstandardResidues()
        if fixer.nonstandardResidues:
            print(f"✅ Found {len(fixer.nonstandardResidues)} non-standard residues")
            for residue in list(fixer.nonstandardResidues)[:5]:
                print(f"    {residue}")
        else:
            print("❌ No non-standard residues found")
            
    except Exception as e:
        print(f"❌ Error testing PDB file: {e}")
        import traceback
        traceback.print_exc()

def main():
    if len(sys.argv) != 2:
        print("Usage: python test_your_structure.py <pdb_file>")
        print("\nExample:")
        print("  python test_your_structure.py src/EasyMD/tests/data/4zgm.pdb")
        print("  python test_your_structure.py your_protein.pdb")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    analyze_pdb_file(pdb_file)

if __name__ == "__main__":
    main()