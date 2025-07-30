"""
Info Runner for EasyMD

This module performs comprehensive PDB structure analysis including chain analysis,
missing residues detection, ligand identification, and structural feature analysis.
"""

import os
import sys
import json
import gzip
from pathlib import Path
from collections import defaultdict, Counter
from typing import List, Dict, Tuple, Set, Any, Optional
import numpy as np
from pdbfixer import PDBFixer


class InfoRunner:
    """
    Performs comprehensive PDB structure analysis.
    
    This class analyzes PDB files to extract detailed structural information
    including chains, missing residues, ligands, water molecules, disulfide bonds,
    and other structural features.
    """
    
    def __init__(self, config):
        """Initialize the InfoRunner with configuration."""
        self.config = config
        self.pdb_file = Path(config.pdb_file)
        self.output_file = Path(config.output) if config.output else None
        self.structure_data = {}
        self.analysis_results = {}
        
        # Color codes for terminal output
        self.colors = {
            'header': '\033[95m', 'blue': '\033[94m', 'green': '\033[92m',
            'yellow': '\033[93m', 'red': '\033[91m', 'bold': '\033[1m',
            'underline': '\033[4m', 'end': '\033[0m'
        } if not config.no_color else {k: '' for k in ['header', 'blue', 'green', 'yellow', 'red', 'bold', 'underline', 'end']}
    
    def run(self):
        """Execute the complete PDB analysis workflow."""
        try:
            print(f"{self.colors['bold']}{'='*60}{self.colors['end']}")
            print(f"{self.colors['header']}🔬 STARTING PDB STRUCTURE ANALYSIS{self.colors['end']}")
            print(f"{self.colors['bold']}{'='*60}{self.colors['end']}")
            
            self._parse_pdb_file()
            self._analyze_structure()
            self._display_results()
            
            if not self.config.no_file and self.output_file:
                self._save_results()
            
            self._display_completion_summary()
            
        except Exception as e:
            print(f"\n{self.colors['red']}❌ Analysis failed with error: {str(e)}{self.colors['end']}")
            sys.exit(1) 
   
    def _parse_pdb_file(self):
        """Parse PDB file and extract basic structural information."""
        print(f"{self.colors['blue']}📁 Parsing PDB file: {self.pdb_file.name}{self.colors['end']}")
        
        self.structure_data = {
            'atoms': [], 'residues': {}, 'chains': {}, 'hetero_atoms': [],
            'water_molecules': [], 'header_info': {}, 'connections': [], 'remarks': []
        }
        
        # Open file (handle gzipped files)
        file_opener = gzip.open if self.pdb_file.suffix.lower() == '.gz' else open
        mode = 'rt' if self.pdb_file.suffix.lower() == '.gz' else 'r'
        
        try:
            with file_opener(self.pdb_file, mode) as f:
                for line_num, line in enumerate(f, 1):
                    try:
                        self._parse_pdb_line(line.strip(), line_num)
                    except Exception as e:
                        if not self.config.quiet:
                            print(f"{self.colors['yellow']}⚠️  Warning: Error parsing line {line_num}: {str(e)}{self.colors['end']}")
            
            print(f"   {self.colors['green']}✅ PDB file parsed successfully{self.colors['end']}")
            print(f"   📊 Found {len(self.structure_data['atoms'])} atoms, "
                  f"{len(self.structure_data['chains'])} chains, "
                  f"{len(self.structure_data['hetero_atoms'])} hetero atoms")
            
        except Exception as e:
            raise RuntimeError(f"Failed to parse PDB file: {str(e)}")
    
    def _parse_pdb_line(self, line: str, line_num: int):
        """Parse a single line from the PDB file."""
        if not line:
            return
        
        record_type = line[:6].strip()
        
        if record_type in ['ATOM', 'HETATM']:
            self._parse_atom_line(line, record_type)
        elif record_type == 'HEADER':
            self._parse_header_line(line)
        elif record_type == 'REMARK':
            self.structure_data['remarks'].append(line[6:].strip())
    
    def _parse_atom_line(self, line: str, record_type: str):
        """Parse ATOM or HETATM line."""
        try:
            atom_data = {
                'record_type': record_type,
                'atom_id': int(line[6:11].strip()),
                'atom_name': line[12:16].strip(),
                'residue_name': line[17:20].strip(),
                'chain_id': line[21:22].strip(),
                'residue_id': int(line[22:26].strip()),
                'insertion_code': line[26:27].strip(),
                'x': float(line[30:38].strip()),
                'y': float(line[38:46].strip()),
                'z': float(line[46:54].strip()),
                'occupancy': float(line[54:60].strip()) if line[54:60].strip() else 1.0,
                'b_factor': float(line[60:66].strip()) if line[60:66].strip() else 0.0,
                'element': line[76:78].strip() if len(line) > 76 else '',
            }
            
            # Add to appropriate lists
            if record_type == 'ATOM':
                self.structure_data['atoms'].append(atom_data)
            else:  # HETATM
                self.structure_data['hetero_atoms'].append(atom_data)
                if atom_data['residue_name'] in ['HOH', 'WAT', 'TIP', 'H2O']:
                    self.structure_data['water_molecules'].append(atom_data)
            
            # Track chains
            chain_id = atom_data['chain_id']
            if chain_id not in self.structure_data['chains']:
                self.structure_data['chains'][chain_id] = {'residues': {}, 'atoms': []}
            
            self.structure_data['chains'][chain_id]['atoms'].append(atom_data)
            
        except (ValueError, IndexError) as e:
            raise ValueError(f"Invalid atom line format: {str(e)}")
    
    def _parse_header_line(self, line: str):
        """Parse HEADER line."""
        try:
            self.structure_data['header_info'] = {
                'classification': line[10:50].strip(),
                'deposition_date': line[50:59].strip(),
                'pdb_id': line[62:66].strip()
            }
        except IndexError:
            pass    

    def _analyze_structure(self):
        """Perform comprehensive structural analysis."""
        print(f"\n{self.colors['blue']}🔬 Performing structural analysis...{self.colors['end']}")
        
        self.analysis_results = {}
        
        if self.config.chains or self.config.all:
            self.analysis_results['chains'] = self._analyze_chains()
        
        if self.config.missing_residues or self.config.all:
            self.analysis_results['missing_residues'] = self._analyze_missing_residues()
        
        if self.config.ligands or self.config.all:
            self.analysis_results['ligands'] = self._analyze_ligands()
        
        if self.config.water or self.config.all:
            self.analysis_results['water'] = self._analyze_water_molecules()
        
        if self.config.disulfide or self.config.all:
            self.analysis_results['disulfide'] = self._analyze_disulfide_bonds()
        
        if self.config.metals or self.config.all:
            self.analysis_results['metals'] = self._analyze_metal_ions()
        
        print(f"   {self.colors['green']}✅ Structural analysis completed{self.colors['end']}")
    
    def _analyze_chains(self) -> Dict[str, Any]:
        """Analyze chain information."""
        print(f"   📊 Analyzing chains...")
        
        chain_analysis = {}
        
        for chain_id, chain_data in self.structure_data['chains'].items():
            protein_residues = []
            for atom in chain_data['atoms']:
                if (atom['record_type'] == 'ATOM' and 
                    atom['residue_name'] not in ['HOH', 'WAT', 'TIP', 'H2O']):
                    res_key = f"{atom['residue_id']}{atom['insertion_code']}"
                    if res_key not in [r['key'] for r in protein_residues]:
                        protein_residues.append({
                            'key': res_key,
                            'id': atom['residue_id'],
                            'name': atom['residue_name'],
                            'insertion': atom['insertion_code']
                        })
            
            protein_residues.sort(key=lambda x: x['id'])
            
            chain_analysis[chain_id] = {
                'total_atoms': len(chain_data['atoms']),
                'protein_residues': len(protein_residues),
                'residue_range': (
                    protein_residues[0]['id'] if protein_residues else None,
                    protein_residues[-1]['id'] if protein_residues else None
                ),
                'sequence': [r['name'] for r in protein_residues],
                'residue_details': protein_residues
            }
        
        return chain_analysis
    
    def _analyze_missing_residues(self) -> Dict[str, Any]:
        """Analyze missing residues using PDBFixer for accurate detection."""
        print(f"   🔍 Analyzing missing residues...")
        
        try:
            # Use PDBFixer to detect missing residues accurately
            fixer = PDBFixer(filename=str(self.pdb_file))
            fixer.findMissingResidues()
            
            missing_analysis = {}
            
            if fixer.missingResidues:
                for chain_tuple, missing_residues in fixer.missingResidues.items():
                    # Convert PDBFixer tuple format to chain ID
                    # chain_tuple is like (0, 0) where first number is chain index
                    chain_index = chain_tuple[0]
                    
                    # Map chain index to actual chain ID from our parsed data
                    chain_ids = list(self.structure_data['chains'].keys())
                    if chain_index < len(chain_ids):
                        chain_id = chain_ids[chain_index]
                    else:
                        chain_id = f"Chain_{chain_index}"
                    
                    # Initialize chain entry if not exists
                    if chain_id not in missing_analysis:
                        missing_analysis[chain_id] = {
                            'gaps': [],
                            'total_missing': 0,
                            'pdbfixer_data': []
                        }
                    
                    # Add this missing segment to the chain
                    if missing_residues:
                        gap_entry = {
                            'start': 'Unknown',  # PDBFixer doesn't give us residue numbers directly
                            'end': 'Unknown',
                            'count': len(missing_residues),
                            'residues': missing_residues,  # List of residue names
                            'type': 'pdbfixer_detected',
                            'chain_tuple': chain_tuple  # Keep original for reference
                        }
                        
                        missing_analysis[chain_id]['gaps'].append(gap_entry)
                        missing_analysis[chain_id]['total_missing'] += len(missing_residues)
                        missing_analysis[chain_id]['pdbfixer_data'].extend(missing_residues)
            
            # Add empty entries for chains with no missing residues
            for chain_id in self.structure_data['chains'].keys():
                if chain_id not in missing_analysis:
                    missing_analysis[chain_id] = {
                        'gaps': [],
                        'total_missing': 0,
                        'pdbfixer_data': []
                    }
            
            return missing_analysis
            
        except Exception as e:
            print(f"   {self.colors['yellow']}⚠️  Warning: PDBFixer analysis failed: {str(e)}{self.colors['end']}")
            print(f"   {self.colors['yellow']}   Falling back to simple gap detection...{self.colors['end']}")
            
            # Fallback to original method if PDBFixer fails
            return self._analyze_missing_residues_fallback()
    
    def _analyze_missing_residues_fallback(self) -> Dict[str, Any]:
        """Fallback method for missing residue detection using simple gap analysis."""
        missing_analysis = {}
        
        for chain_id, chain_data in self.analysis_results.get('chains', {}).items():
            missing_residues = []
            residue_details = chain_data['residue_details']
            
            if len(residue_details) < 2:
                missing_analysis[chain_id] = {
                    'gaps': [],
                    'total_missing': 0,
                    'pdbfixer_data': []
                }
                continue
            
            for i in range(len(residue_details) - 1):
                current_id = residue_details[i]['id']
                next_id = residue_details[i + 1]['id']
                
                gap_size = next_id - current_id - 1
                if gap_size >= self.config.missing_threshold:
                    missing_residues.append({
                        'start': current_id + 1,
                        'end': next_id - 1,
                        'count': gap_size,
                        'after_residue': residue_details[i]['name'],
                        'before_residue': residue_details[i + 1]['name'],
                        'type': 'gap'
                    })
            
            missing_analysis[chain_id] = {
                'gaps': missing_residues,
                'total_missing': sum(gap['count'] for gap in missing_residues),
                'pdbfixer_data': []
            }
        
        return missing_analysis
    
    def _analyze_ligands(self) -> Dict[str, Any]:
        """Analyze ligands and small molecules."""
        print(f"   💊 Analyzing ligands and small molecules...")
        
        exclude_list = {
            'HOH', 'WAT', 'TIP', 'H2O',  # Water
            'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU',  # Common ions
            'SO4', 'PO4', 'NO3', 'CO3',  # Common salts
            'GOL', 'EDO', 'PEG', 'DMS', 'DMSO'  # Common crystallization agents
        }
        
        ligands = {}
        crystallization_agents = {}
        
        for atom in self.structure_data['hetero_atoms']:
            res_name = atom['residue_name']
            
            if res_name in ['HOH', 'WAT', 'TIP', 'H2O']:
                continue
            
            res_key = f"{atom['chain_id']}_{res_name}_{atom['residue_id']}"
            
            if res_name in exclude_list:
                if res_key not in crystallization_agents:
                    crystallization_agents[res_key] = {
                        'name': res_name, 'chain': atom['chain_id'],
                        'residue_id': atom['residue_id'], 'atoms': []
                    }
                crystallization_agents[res_key]['atoms'].append(atom)
            else:
                if res_key not in ligands:
                    ligands[res_key] = {
                        'name': res_name, 'chain': atom['chain_id'],
                        'residue_id': atom['residue_id'], 'atoms': []
                    }
                ligands[res_key]['atoms'].append(atom)
        
        # Calculate properties
        for ligand_data in ligands.values():
            atoms = ligand_data['atoms']
            ligand_data['atom_count'] = len(atoms)
            ligand_data['elements'] = Counter(atom['element'] for atom in atoms if atom['element'])
        
        return {
            'ligands': ligands,
            'crystallization_agents': crystallization_agents,
            'ligand_count': len(ligands),
            'agent_count': len(crystallization_agents)
        }    
   
    def _analyze_water_molecules(self) -> Dict[str, Any]:
        """Analyze water molecules."""
        print(f"   💧 Analyzing water molecules...")
        
        water_molecules = self.structure_data['water_molecules']
        water_by_chain = defaultdict(list)
        for water in water_molecules:
            water_by_chain[water['chain_id']].append(water)
        
        return {
            'total_count': len(water_molecules),
            'by_chain': dict(water_by_chain),
            'chain_counts': {chain: len(waters) for chain, waters in water_by_chain.items()},
            'show_details': len(water_molecules) >= self.config.water_threshold
        }
    
    def _analyze_disulfide_bonds(self) -> Dict[str, Any]:
        """Analyze potential disulfide bonds."""
        print(f"   🔗 Analyzing disulfide bonds...")
        
        sulfur_atoms = []
        for atom in self.structure_data['atoms']:
            if (atom['element'] in ['S', 'SG'] or atom['atom_name'] == 'SG') and \
               atom['residue_name'] == 'CYS':
                sulfur_atoms.append(atom)
        
        disulfide_bonds = []
        for i, atom1 in enumerate(sulfur_atoms):
            for atom2 in sulfur_atoms[i+1:]:
                dist = np.sqrt(
                    (atom1['x'] - atom2['x'])**2 +
                    (atom1['y'] - atom2['y'])**2 +
                    (atom1['z'] - atom2['z'])**2
                )
                
                if dist <= self.config.disulfide_distance:
                    disulfide_bonds.append({
                        'residue1': {
                            'chain': atom1['chain_id'],
                            'residue_id': atom1['residue_id'],
                            'residue_name': atom1['residue_name']
                        },
                        'residue2': {
                            'chain': atom2['chain_id'],
                            'residue_id': atom2['residue_id'],
                            'residue_name': atom2['residue_name']
                        },
                        'distance': round(dist, 2)
                    })
        
        return {
            'cysteine_count': len(sulfur_atoms),
            'bonds': disulfide_bonds,
            'bond_count': len(disulfide_bonds)
        }
    
    def _analyze_metal_ions(self) -> Dict[str, Any]:
        """Analyze metal ions."""
        print(f"   ⚛️  Analyzing metal ions...")
        
        metal_elements = {
            'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'NI', 'CO', 'CD', 'HG',
            'NA', 'K', 'LI', 'RB', 'CS', 'AL', 'CR', 'MO', 'W', 'V'
        }
        
        metals = []
        for atom in self.structure_data['hetero_atoms']:
            if atom['element'] in metal_elements or atom['residue_name'] in metal_elements:
                metals.append({
                    'element': atom['element'] or atom['residue_name'],
                    'chain': atom['chain_id'],
                    'residue_id': atom['residue_id'],
                    'coordinates': [atom['x'], atom['y'], atom['z']],
                    'occupancy': atom['occupancy'],
                    'b_factor': atom['b_factor']
                })
        
        metals_by_element = defaultdict(list)
        for metal in metals:
            metals_by_element[metal['element']].append(metal)
        
        return {
            'metals': metals,
            'by_element': dict(metals_by_element),
            'element_counts': {element: len(metals) for element, metals in metals_by_element.items()},
            'total_count': len(metals)
        }    
  
    def _display_results(self):
        """Display analysis results to terminal."""
        if self.config.format == "json":
            self._display_json_results()
        elif self.config.format == "summary":
            self._display_summary_results()
        else:  # detailed
            self._display_detailed_results()
    
    def _display_detailed_results(self):
        """Display detailed analysis results."""
        print(f"\n{self.colors['bold']}{'='*60}{self.colors['end']}")
        print(f"{self.colors['header']}📋 PDB STRUCTURE ANALYSIS RESULTS{self.colors['end']}")
        print(f"{self.colors['bold']}{'='*60}{self.colors['end']}")
        
        # Header information
        if self.structure_data['header_info']:
            header = self.structure_data['header_info']
            print(f"\n{self.colors['bold']}📄 Header Information:{self.colors['end']}")
            print(f"   PDB ID: {header.get('pdb_id', 'Unknown')}")
            print(f"   Classification: {header.get('classification', 'Unknown')}")
            print(f"   Deposition Date: {header.get('deposition_date', 'Unknown')}")
        
        # Chain analysis
        if 'chains' in self.analysis_results:
            chains = self.analysis_results['chains']
            print(f"\n{self.colors['bold']}🔗 Chain Analysis:{self.colors['end']}")
            print(f"   Total chains: {len(chains)}")
            
            for chain_id, chain_data in chains.items():
                print(f"\n   {self.colors['blue']}Chain {chain_id}:{self.colors['end']}")
                print(f"      Total atoms: {chain_data['total_atoms']}")
                print(f"      Protein residues: {chain_data['protein_residues']}")
                
                if chain_data['residue_range'][0] is not None:
                    print(f"      Residue range: {chain_data['residue_range'][0]} - {chain_data['residue_range'][1]}")
                    
                    sequence = chain_data['sequence']
                    if sequence:
                        seq_display = ' '.join(sequence[:20])
                        if len(sequence) > 20:
                            seq_display += f" ... (+{len(sequence) - 20} more)"
                        print(f"      Sequence: {seq_display}")
        
        # Missing residues
        if 'missing_residues' in self.analysis_results:
            missing = self.analysis_results['missing_residues']
            print(f"\n{self.colors['bold']}❓ Missing Residues:{self.colors['end']}")
            
            total_missing = sum(chain_data['total_missing'] for chain_data in missing.values())
            print(f"   Total missing residues: {total_missing}")
            
            for chain_id, chain_data in missing.items():
                if chain_data['gaps']:
                    print(f"\n   {self.colors['yellow']}Chain {chain_id}:{self.colors['end']}")
                    for gap in chain_data['gaps']:
                        if gap['type'] == 'gap':
                            # Fallback format
                            print(f"      Gap: {gap['start']}-{gap['end']} ({gap['count']} residues)")
                            print(f"           Between {gap['after_residue']} and {gap['before_residue']}")
                        else:
                            # PDBFixer format
                            residue_list = ' '.join(gap['residues'])
                            print(f"      Missing residues: ({gap['count']} residues)")
                            print(f"           Residues: {residue_list}")
                            print(f"           Chain tuple: {gap.get('chain_tuple', 'Unknown')}")
        
        # Ligands
        if 'ligands' in self.analysis_results:
            ligands = self.analysis_results['ligands']
            print(f"\n{self.colors['bold']}💊 Ligands and Small Molecules:{self.colors['end']}")
            print(f"   Ligands found: {ligands['ligand_count']}")
            print(f"   Crystallization agents: {ligands['agent_count']}")
            
            if ligands['ligands']:
                print(f"\n   {self.colors['green']}Ligands:{self.colors['end']}")
                for ligand_key, ligand_data in ligands['ligands'].items():
                    print(f"      {ligand_data['name']} (Chain {ligand_data['chain']}, Residue {ligand_data['residue_id']})")
                    print(f"         Atoms: {ligand_data['atom_count']}")
                    if ligand_data['elements']:
                        elements = ', '.join(f"{elem}:{count}" for elem, count in ligand_data['elements'].items())
                        print(f"         Elements: {elements}")
        
        # Water molecules
        if 'water' in self.analysis_results:
            water = self.analysis_results['water']
            print(f"\n{self.colors['bold']}💧 Water Molecules:{self.colors['end']}")
            print(f"   Total water molecules: {water['total_count']}")
            
            if water['chain_counts']:
                print(f"   Distribution by chain:")
                for chain, count in water['chain_counts'].items():
                    print(f"      Chain {chain}: {count} waters")
        
        # Disulfide bonds
        if 'disulfide' in self.analysis_results:
            disulfide = self.analysis_results['disulfide']
            print(f"\n{self.colors['bold']}🔗 Disulfide Bonds:{self.colors['end']}")
            print(f"   Cysteine residues: {disulfide['cysteine_count']}")
            print(f"   Potential disulfide bonds: {disulfide['bond_count']}")
            
            if disulfide['bonds']:
                print(f"\n   {self.colors['green']}Detected bonds:{self.colors['end']}")
                for bond in disulfide['bonds']:
                    res1 = bond['residue1']
                    res2 = bond['residue2']
                    print(f"      {res1['chain']}:{res1['residue_name']}{res1['residue_id']} - "
                          f"{res2['chain']}:{res2['residue_name']}{res2['residue_id']} "
                          f"({bond['distance']} Å)")
        
        # Metal ions
        if 'metals' in self.analysis_results:
            metals = self.analysis_results['metals']
            print(f"\n{self.colors['bold']}⚛️  Metal Ions:{self.colors['end']}")
            print(f"   Total metal ions: {metals['total_count']}")
            
            if metals['element_counts']:
                print(f"   Distribution by element:")
                for element, count in metals['element_counts'].items():
                    print(f"      {element}: {count} ions")
    
    def _display_summary_results(self):
        """Display summary results."""
        print(f"\n{self.colors['bold']}📋 PDB STRUCTURE SUMMARY{self.colors['end']}")
        print(f"File: {self.pdb_file.name}")
        
        if 'chains' in self.analysis_results:
            chains = self.analysis_results['chains']
            print(f"Chains: {len(chains)}")
            total_residues = sum(chain['protein_residues'] for chain in chains.values())
            print(f"Total protein residues: {total_residues}")
        
        if 'missing_residues' in self.analysis_results:
            missing = self.analysis_results['missing_residues']
            total_missing = sum(chain_data['total_missing'] for chain_data in missing.values())
            print(f"Missing residues: {total_missing}")
        
        if 'ligands' in self.analysis_results:
            ligands = self.analysis_results['ligands']
            print(f"Ligands: {ligands['ligand_count']}")
        
        if 'water' in self.analysis_results:
            water = self.analysis_results['water']
            print(f"Water molecules: {water['total_count']}")
        
        if 'disulfide' in self.analysis_results:
            disulfide = self.analysis_results['disulfide']
            print(f"Disulfide bonds: {disulfide['bond_count']}")
        
        if 'metals' in self.analysis_results:
            metals = self.analysis_results['metals']
            print(f"Metal ions: {metals['total_count']}")
    
    def _display_json_results(self):
        """Display results in JSON format."""
        output_data = {
            'pdb_file': str(self.pdb_file),
            'header_info': self.structure_data['header_info'],
            'analysis_results': self.analysis_results
        }
        print(json.dumps(output_data, indent=2))
    
    def _save_results(self):
        """Save analysis results to file."""
        print(f"\n{self.colors['blue']}💾 Saving results to {self.output_file.name}...{self.colors['end']}")
        
        try:
            with open(self.output_file, 'w') as f:
                # Write header
                f.write("="*80 + "\n")
                f.write("PDB STRUCTURE ANALYSIS REPORT\n")
                f.write("="*80 + "\n")
                f.write(f"Generated by EasyMD Info Tool\n")
                f.write(f"PDB File: {self.pdb_file}\n")
                f.write(f"Analysis Date: {self._get_current_timestamp()}\n")
                f.write("="*80 + "\n\n")
                
                # Write header information
                if self.structure_data['header_info']:
                    header = self.structure_data['header_info']
                    f.write("HEADER INFORMATION\n")
                    f.write("-" * 40 + "\n")
                    f.write(f"PDB ID: {header.get('pdb_id', 'Unknown')}\n")
                    f.write(f"Classification: {header.get('classification', 'Unknown')}\n")
                    f.write(f"Deposition Date: {header.get('deposition_date', 'Unknown')}\n\n")
                
                # Write analysis results
                if 'chains' in self.analysis_results:
                    self._write_chain_analysis(f)
                
                if 'missing_residues' in self.analysis_results:
                    self._write_missing_residues(f)
                
                if 'ligands' in self.analysis_results:
                    self._write_ligand_analysis(f)
                
                if 'water' in self.analysis_results:
                    self._write_water_analysis(f)
                
                if 'disulfide' in self.analysis_results:
                    self._write_disulfide_analysis(f)
                
                if 'metals' in self.analysis_results:
                    self._write_metal_analysis(f)
                
                # Write footer
                f.write("\n" + "="*80 + "\n")
                f.write("END OF REPORT\n")
                f.write("="*80 + "\n")
            
            print(f"   {self.colors['green']}✅ Results saved to {self.output_file}{self.colors['end']}")
            
        except Exception as e:
            print(f"   {self.colors['red']}❌ Failed to save results: {str(e)}{self.colors['end']}")
    
    def _write_chain_analysis(self, f):
        """Write chain analysis to file."""
        chains = self.analysis_results['chains']
        f.write("CHAIN ANALYSIS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total chains: {len(chains)}\n\n")
        
        for chain_id, chain_data in chains.items():
            f.write(f"Chain {chain_id}:\n")
            f.write(f"  Total atoms: {chain_data['total_atoms']}\n")
            f.write(f"  Protein residues: {chain_data['protein_residues']}\n")
            
            if chain_data['residue_range'][0] is not None:
                f.write(f"  Residue range: {chain_data['residue_range'][0]} - {chain_data['residue_range'][1]}\n")
                
                # Write full sequence
                sequence = chain_data['sequence']
                if sequence:
                    f.write(f"  Sequence: {' '.join(sequence)}\n")
            f.write("\n")
    
    def _write_missing_residues(self, f):
        """Write missing residue analysis to file."""
        missing = self.analysis_results['missing_residues']
        f.write("MISSING RESIDUES\n")
        f.write("-" * 40 + "\n")
        
        total_missing = sum(chain_data['total_missing'] for chain_data in missing.values())
        f.write(f"Total missing residues: {total_missing}\n\n")
        
        for chain_id, chain_data in missing.items():
            if chain_data['gaps']:
                f.write(f"Chain {chain_id}:\n")
                for gap in chain_data['gaps']:
                    if gap['type'] == 'gap':
                        # Fallback format
                        f.write(f"  Gap: {gap['start']}-{gap['end']} ({gap['count']} residues)\n")
                        f.write(f"       Between {gap['after_residue']} and {gap['before_residue']}\n")
                    else:
                        # PDBFixer format
                        residue_list = ' '.join(gap['residues'])
                        f.write(f"  Missing residues: ({gap['count']} residues)\n")
                        f.write(f"       Residues: {residue_list}\n")
                        f.write(f"       Chain tuple: {gap.get('chain_tuple', 'Unknown')}\n")
                f.write("\n")
    
    def _write_ligand_analysis(self, f):
        """Write ligand analysis to file."""
        ligands = self.analysis_results['ligands']
        f.write("LIGANDS AND SMALL MOLECULES\n")
        f.write("-" * 40 + "\n")
        f.write(f"Ligands found: {ligands['ligand_count']}\n")
        f.write(f"Crystallization agents: {ligands['agent_count']}\n\n")
        
        if ligands['ligands']:
            f.write("Ligands:\n")
            for ligand_key, ligand_data in ligands['ligands'].items():
                f.write(f"  {ligand_data['name']} (Chain {ligand_data['chain']}, Residue {ligand_data['residue_id']})\n")
                f.write(f"    Atoms: {ligand_data['atom_count']}\n")
                if ligand_data['elements']:
                    elements = ', '.join(f"{elem}:{count}" for elem, count in ligand_data['elements'].items())
                    f.write(f"    Elements: {elements}\n")
            f.write("\n")
        
        if ligands['crystallization_agents']:
            f.write("Crystallization Agents:\n")
            agent_counts = Counter(agent['name'] for agent in ligands['crystallization_agents'].values())
            for agent_name, count in agent_counts.items():
                f.write(f"  {agent_name}: {count} molecules\n")
            f.write("\n")
    
    def _write_water_analysis(self, f):
        """Write water analysis to file."""
        water = self.analysis_results['water']
        f.write("WATER MOLECULES\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total water molecules: {water['total_count']}\n")
        
        if water['chain_counts']:
            f.write("Distribution by chain:\n")
            for chain, count in water['chain_counts'].items():
                f.write(f"  Chain {chain}: {count} waters\n")
        f.write("\n")
    
    def _write_disulfide_analysis(self, f):
        """Write disulfide analysis to file."""
        disulfide = self.analysis_results['disulfide']
        f.write("DISULFIDE BONDS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Cysteine residues: {disulfide['cysteine_count']}\n")
        f.write(f"Potential disulfide bonds: {disulfide['bond_count']}\n")
        
        if disulfide['bonds']:
            f.write("\nDetected bonds:\n")
            for bond in disulfide['bonds']:
                res1 = bond['residue1']
                res2 = bond['residue2']
                f.write(f"  {res1['chain']}:{res1['residue_name']}{res1['residue_id']} - "
                       f"{res2['chain']}:{res2['residue_name']}{res2['residue_id']} "
                       f"({bond['distance']} Å)\n")
        f.write("\n")
    
    def _write_metal_analysis(self, f):
        """Write metal analysis to file."""
        metals = self.analysis_results['metals']
        f.write("METAL IONS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total metal ions: {metals['total_count']}\n")
        
        if metals['element_counts']:
            f.write("Distribution by element:\n")
            for element, count in metals['element_counts'].items():
                f.write(f"  {element}: {count} ions\n")
        f.write("\n")
    
    def _get_current_timestamp(self):
        """Get current timestamp for report."""
        from datetime import datetime
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    def _display_completion_summary(self):
        """Display analysis completion summary."""
        print(f"\n{self.colors['bold']}{'='*60}{self.colors['end']}")
        print(f"{self.colors['green']}🎉 PDB ANALYSIS COMPLETED SUCCESSFULLY{self.colors['end']}")
        print(f"{self.colors['bold']}{'='*60}{self.colors['end']}")
        
        print(f"\n📁 Analyzed file: {self.pdb_file}")
        if not self.config.no_file and self.output_file:
            print(f"📄 Report saved to: {self.output_file}")
        
        # Summary statistics
        stats = []
        if 'chains' in self.analysis_results:
            chains = self.analysis_results['chains']
            stats.append(f"{len(chains)} chains")
        
        if 'ligands' in self.analysis_results:
            ligands = self.analysis_results['ligands']
            if ligands['ligand_count'] > 0:
                stats.append(f"{ligands['ligand_count']} ligands")
        
        if 'missing_residues' in self.analysis_results:
            missing = self.analysis_results['missing_residues']
            total_missing = sum(chain_data['total_missing'] for chain_data in missing.values())
            if total_missing > 0:
                stats.append(f"{total_missing} missing residues")
        
        if 'disulfide' in self.analysis_results:
            disulfide = self.analysis_results['disulfide']
            if disulfide['bond_count'] > 0:
                stats.append(f"{disulfide['bond_count']} disulfide bonds")
        
        if stats:
            print(f"📊 Key findings: {', '.join(stats)}")
        
        print(f"\n{self.colors['bold']}{'='*60}{self.colors['end']}")
        print("Analysis complete! Check the output for detailed structural information.")
        print(f"{self.colors['bold']}{'='*60}{self.colors['end']}")