# Info Component

## Overview
The Info component provides comprehensive PDB structure analysis, offering detailed insights into protein structures before MD simulation setup. It features advanced missing residue detection using PDBFixer, ligand identification, water analysis, and structural validation, making it an essential tool for structure preparation and quality assessment.

## Architecture

### Component Structure
```
info/
├── __init__.py        # Module initialization
├── infoManager.py     # Argument parsing and validation
└── infoRunner.py      # Structure analysis and reporting
```

### Key Classes
- **InfoManager**: Handles command-line arguments and validation
- **InfoRunner**: Performs PDB analysis and generates reports

## InfoManager - Argument Management

### Command-Line Interface
```bash
# Basic structure analysis
python -m EasyMD info protein.pdb

# Custom output format
python -m EasyMD info protein.pdb --format summary --output analysis.info

# JSON output for programmatic use
python -m EasyMD info protein.pdb --format json --no-file
```

### Argument Categories

#### 1. Input/Output Parameters
```python
io_group.add_argument("pdb_file", type=str, help="Path to PDB file to analyze")
io_group.add_argument("-o", "--output", type=str, default=None,
                     help="Output file for detailed report (default: <pdb_name>.info)")
io_group.add_argument("--no-file", action='store_true',
                     help="Don't create output file, only display to terminal")
```

#### 2. Analysis Options
- **Chains**: Chain composition and sequence analysis
- **Missing Residues**: PDBFixer-based missing residue detection
- **Ligands**: Small molecule and ligand identification
- **Water**: Water molecule distribution analysis
- **Disulfide**: Disulfide bond prediction
- **Metals**: Metal ion detection
- **Modifications**: Modified residue identification

#### 3. Analysis Parameters
```python
params_group.add_argument("--disulfide-distance", type=float, default=2.5,
                         help="Maximum distance for disulfide bond detection in Angstroms")
params_group.add_argument("--water-threshold", type=int, default=10,
                         help="Minimum number of water molecules to report details")
params_group.add_argument("--missing-threshold", type=int, default=3,
                         help="Minimum gap size to report as missing residues")
```

#### 4. Output Formatting
- **Format Options**: Detailed, summary, or JSON output
- **Color Control**: Terminal color output management
- **Verbosity**: Quiet mode for minimal output

### Validation System
Comprehensive validation includes:
- **File Existence**: PDB file accessibility checks
- **Format Validation**: PDB format verification
- **Parameter Ranges**: Analysis parameter validation
- **Output Permissions**: Write access verification

## InfoRunner - Structure Analysis

### Initialization and Setup
```python
class InfoRunner:
    """
    Performs comprehensive PDB structure analysis.
    
    Attributes:
        config: Configuration object containing analysis parameters
        pdb_file (Path): Path to PDB file
        output_file (Path): Output file path
        structure_data (dict): Parsed PDB structure data
        analysis_results (dict): Analysis results
        colors (dict): Terminal color codes
    """
```

### PDB File Parsing

#### 1. File Format Support
- **Standard PDB**: Regular PDB files
- **Compressed PDB**: Gzipped PDB files (.pdb.gz)
- **Error Handling**: Graceful handling of format issues

#### 2. Record Type Processing
```python
def _parse_pdb_line(self, line: str, line_num: int):
    """Parse a single line from the PDB file."""
    record_type = line[:6].strip()
    
    if record_type in ['ATOM', 'HETATM']:
        self._parse_atom_line(line, record_type)
    elif record_type == 'HEADER':
        self._parse_header_line(line)
    elif record_type == 'REMARK':
        self.structure_data['remarks'].append(line[6:].strip())
```

#### 3. Data Structure Organization
```python
self.structure_data = {
    'atoms': [],           # Protein atoms
    'residues': {},        # Residue information
    'chains': {},          # Chain organization
    'hetero_atoms': [],    # Ligands and other molecules
    'water_molecules': [], # Water molecules
    'header_info': {},     # PDB header information
    'connections': [],     # Connectivity information
    'remarks': []          # PDB remarks
}
```

## Analysis Methods

### 1. Missing Residue Detection (PDBFixer Integration)

#### Advanced Missing Residue Analysis
```python
def _analyze_missing_residues_pdbfixer(self) -> Dict[str, Any]:
    """Use PDBFixer for accurate missing residue detection."""
    try:
        fixer = PDBFixer(filename=str(self.pdb_file))
        fixer.findMissingResidues()
        
        if not fixer.missingResidues:
            return {'total_missing': 0, 'chains': {}, 'method': 'PDBFixer'}
        
        # Process missing residues by chain
        missing_by_chain = {}
        total_missing = 0
        
        for chain_key, missing_residues in fixer.missingResidues.items():
            chain_id = self._convert_chain_key_to_id(chain_key)
            
            # Categorize missing residues
            n_terminal, c_terminal, internal = self._categorize_missing_residues(
                missing_residues, chain_id)
            
            missing_by_chain[chain_id] = {
                'total': len(missing_residues),
                'n_terminal': n_terminal,
                'c_terminal': c_terminal,
                'internal': internal,
                'residues': missing_residues
            }
            total_missing += len(missing_residues)
```

**Key Features:**
- **PDBFixer Integration**: Uses same detection method as sysGenerator
- **Chain-Specific Analysis**: Detailed per-chain missing residue information
- **Terminal Classification**: Distinguishes N-terminal, C-terminal, and internal gaps
- **Residue Details**: Provides specific residue names and positions

**Accuracy Improvement:**
- **Before**: Simple gap detection (often 0 missing residues)
- **After**: PDBFixer detection (e.g., 25 missing residues for 4zgm.pdb)

#### Fallback Gap Detection
```python
def _analyze_missing_residues_gaps(self) -> Dict[str, Any]:
    """Fallback method using simple gap detection."""
    # Identify gaps in residue numbering
    # Less accurate but always available
```

### 2. Chain Analysis
```python
def _analyze_chains(self) -> Dict[str, Any]:
    """Analyze chain composition and properties."""
    chain_info = {}
    
    for chain_id, chain_data in self.structure_data['chains'].items():
        # Calculate chain statistics
        residue_count = len(set(atom['residue_id'] for atom in chain_data['atoms']))
        atom_count = len(chain_data['atoms'])
        
        # Determine chain type (protein, DNA, RNA, etc.)
        chain_type = self._determine_chain_type(chain_data['atoms'])
        
        # Get sequence information
        sequence = self._extract_sequence(chain_data['atoms'])
        
        chain_info[chain_id] = {
            'type': chain_type,
            'residue_count': residue_count,
            'atom_count': atom_count,
            'sequence': sequence,
            'first_residue': min(atom['residue_id'] for atom in chain_data['atoms']),
            'last_residue': max(atom['residue_id'] for atom in chain_data['atoms'])
        }
```

**Analysis Features:**
- **Chain Type Detection**: Protein, DNA, RNA, or other
- **Sequence Extraction**: Amino acid or nucleotide sequences
- **Residue Counting**: Accurate residue and atom counts
- **Range Analysis**: First and last residue identification

### 3. Ligand and Small Molecule Analysis
```python
def _analyze_ligands(self) -> Dict[str, Any]:
    """Identify and analyze ligands and small molecules."""
    ligands = {}
    
    # Common solvent and buffer molecules to exclude
    common_solvents = {'HOH', 'WAT', 'TIP', 'H2O', 'SO4', 'PO4', 'CL', 'NA', 'MG', 'CA'}
    
    for hetero_atom in self.structure_data['hetero_atoms']:
        residue_name = hetero_atom['residue_name']
        
        if residue_name not in common_solvents:
            if residue_name not in ligands:
                ligands[residue_name] = {
                    'count': 0,
                    'chains': set(),
                    'atoms': [],
                    'molecular_weight': 0
                }
            
            ligands[residue_name]['count'] += 1
            ligands[residue_name]['chains'].add(hetero_atom['chain_id'])
            ligands[residue_name]['atoms'].append(hetero_atom)
```

**Features:**
- **Automatic Detection**: Identifies non-standard residues as potential ligands
- **Solvent Filtering**: Excludes common solvents and ions
- **Multi-Chain Support**: Tracks ligands across multiple chains
- **Molecular Properties**: Calculates basic molecular properties

### 4. Water Molecule Analysis
```python
def _analyze_water(self) -> Dict[str, Any]:
    """Analyze water molecule distribution."""
    water_info = {
        'total_count': len(self.structure_data['water_molecules']),
        'by_chain': {},
        'occupancy_stats': {},
        'b_factor_stats': {}
    }
    
    if water_info['total_count'] > 0:
        # Analyze water distribution by chain
        # Calculate occupancy and B-factor statistics
        # Identify high-occupancy waters (likely structural)
```

**Analysis Features:**
- **Distribution Analysis**: Water molecules per chain
- **Quality Assessment**: Occupancy and B-factor statistics
- **Structural Waters**: Identification of likely structural waters
- **Crystallographic Quality**: Assessment of water quality

### 5. Disulfide Bond Prediction
```python
def _analyze_disulfide_bonds(self) -> Dict[str, Any]:
    """Predict potential disulfide bonds."""
    cysteine_residues = self._find_cysteine_residues()
    potential_bonds = []
    
    for i, cys1 in enumerate(cysteine_residues):
        for cys2 in cysteine_residues[i+1:]:
            # Calculate distance between sulfur atoms
            distance = self._calculate_distance(cys1['sg_coords'], cys2['sg_coords'])
            
            if distance <= self.config.disulfide_distance:
                potential_bonds.append({
                    'residue1': cys1,
                    'residue2': cys2,
                    'distance': distance,
                    'confidence': self._assess_bond_confidence(distance)
                })
```

**Features:**
- **Distance-Based Prediction**: Uses configurable distance threshold
- **Confidence Assessment**: Evaluates bond likelihood
- **Cross-Chain Bonds**: Detects inter-chain disulfide bonds
- **Structural Validation**: Considers geometric constraints

### 6. Metal Ion Detection
```python
def _analyze_metal_ions(self) -> Dict[str, Any]:
    """Identify and analyze metal ions."""
    common_metals = {'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'NI', 'CO', 'K', 'NA'}
    metal_ions = {}
    
    for hetero_atom in self.structure_data['hetero_atoms']:
        if hetero_atom['residue_name'] in common_metals:
            # Analyze metal coordination
            # Identify potential binding sites
            # Calculate coordination geometry
```

**Analysis Features:**
- **Common Metal Detection**: Identifies biologically relevant metals
- **Coordination Analysis**: Examines metal binding environment
- **Binding Site Identification**: Locates potential metal binding sites
- **Structural Role Assessment**: Evaluates metal structural importance

## Output Generation

### 1. Terminal Display
```python
def _display_results(self):
    """Display analysis results to terminal with color formatting."""
    print(f"{self.colors['header']}📊 STRUCTURE ANALYSIS RESULTS{self.colors['end']}")
    
    # Display each analysis section with appropriate formatting
    self._display_chain_info()
    self._display_missing_residues()
    self._display_ligand_info()
    # ... other sections
```

**Features:**
- **Color-Coded Output**: Visual distinction of different information types
- **Hierarchical Display**: Organized presentation of results
- **Statistical Summaries**: Key metrics prominently displayed
- **Warning Highlights**: Important issues clearly marked

### 2. File Output
```python
def _save_results(self):
    """Save detailed analysis results to file."""
    if self.config.format == 'json':
        self._save_json_results()
    elif self.config.format == 'summary':
        self._save_summary_results()
    else:
        self._save_detailed_results()
```

**Output Formats:**
- **Detailed**: Comprehensive human-readable report
- **Summary**: Concise overview of key findings
- **JSON**: Machine-readable structured data

### 3. Report Structure
```
PDB Structure Analysis Report
=============================

File Information:
- PDB ID: 4ZGM
- Classification: HYDROLASE/HYDROLASE INHIBITOR
- Deposition Date: 2015-04-15

Chain Analysis:
- Chain A: Protein, 245 residues, 1,876 atoms
- Chain B: Protein, 245 residues, 1,876 atoms

Missing Residues (PDBFixer Detection):
- Total: 25 missing residues
- Chain A: 22 missing (5 N-terminal, 17 C-terminal)
- Chain B: 3 missing (3 N-terminal)

Ligands and Small Molecules:
- ATP: 2 instances (Chains A, B)
- MG: 4 instances (Chains A, B)

Water Molecules:
- Total: 156 water molecules
- Average occupancy: 0.65
- High-occupancy waters: 23

Disulfide Bonds:
- Predicted bonds: 2
- Intra-chain: 2, Inter-chain: 0

Quality Assessment:
- Structure completeness: 95.2%
- Average B-factor: 35.4 Å²
- Resolution: 2.1 Å (from header)
```

## Integration Features

### 1. PDBFixer Integration
```python
# Same missing residue detection as sysGenerator
fixer = PDBFixer(filename=str(self.pdb_file))
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.findNonstandardResidues()
```

**Benefits:**
- **Consistency**: Same detection method as system preparation
- **Accuracy**: Superior to simple gap detection
- **Reliability**: Handles complex missing residue patterns

### 2. Error Handling and Validation
```python
def _validate_coordinates(self, atoms):
    """Validate atom coordinates for common issues."""
    for atom in atoms:
        if any(coord > 9999 or coord < -999 for coord in [atom['x'], atom['y'], atom['z']]):
            self.warnings.append(f"Unusual coordinates for atom {atom['atom_id']}")
```

**Validation Features:**
- **Coordinate Validation**: Check for reasonable coordinate values
- **Format Validation**: Ensure proper PDB format compliance
- **Completeness Checks**: Identify incomplete or corrupted structures

### 3. Performance Optimization
- **Memory Efficient**: Processes large PDB files without excessive memory use
- **Fast Parsing**: Optimized PDB parsing algorithms
- **Lazy Loading**: Loads only necessary data for analysis

## Usage Examples

### Basic Analysis
```bash
# Analyze structure with default settings
python -m EasyMD info protein.pdb

# Generate summary report
python -m EasyMD info protein.pdb --format summary

# JSON output for scripts
python -m EasyMD info protein.pdb --format json --no-file
```

### Advanced Options
```bash
# Custom disulfide bond detection
python -m EasyMD info protein.pdb --disulfide-distance 3.0

# Detailed water analysis
python -m EasyMD info protein.pdb --water-threshold 5

# Quiet mode with custom output
python -m EasyMD info protein.pdb --quiet --output detailed_analysis.info
```

## Benefits

1. **Comprehensive Analysis**: Complete structural assessment before simulation
2. **PDBFixer Integration**: Accurate missing residue detection
3. **User-Friendly Output**: Clear, informative reports
4. **Multiple Formats**: Flexible output options for different needs
5. **Quality Assessment**: Identifies potential structural issues
6. **Preparation Guidance**: Helps users understand structure requirements
7. **Consistency**: Uses same methods as system preparation tools

## Real-World Impact

### Before Enhancement
- Simple gap detection: 0 missing residues detected
- Limited structural insight
- Potential simulation failures due to undetected issues

### After Enhancement
- PDBFixer detection: 25 missing residues accurately identified
- Comprehensive structural analysis
- Informed decision-making for simulation setup
- Reduced simulation failures through better preparation

The Info component serves as an essential quality control and analysis tool, providing researchers with detailed insights into their protein structures and helping ensure successful molecular dynamics simulations through comprehensive structural assessment.