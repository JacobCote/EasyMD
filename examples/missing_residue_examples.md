# Missing Residue Handling in EasyMD

EasyMD provides fine-grained control over how missing residues are handled during protein structure preparation. This is crucial for obtaining reliable simulation results, as different strategies work better for different types of structures and research goals.

## 🎯 **Missing Residue Strategies**

### 1. **Auto (Default)**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate --missing-residues auto
```
- **Behavior**: Uses PDBFixer's automatic detection
- **Best for**: General use when you trust the structure quality
- **Note**: May add many residues, including uncertain loop regions

### 2. **None (Skip All)**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate --missing-residues none
```
- **Behavior**: Skips all missing residues
- **Best for**: High-resolution structures with minimal missing data
- **Warning**: May leave gaps that could affect simulation stability

### 3. **Non-Terminal Only**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate --missing-residues non-terminal
```
- **Behavior**: Only adds missing residues in the middle of chains
- **Best for**: When you want to preserve chain integrity but avoid uncertain terminal regions
- **Use case**: Structures with missing internal loops but uncertain termini

### 4. **Terminal Only**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate --missing-residues terminal-only --max-terminal-residues 3
```
- **Behavior**: Only adds missing residues at chain ends
- **Best for**: Structures with well-defined cores but missing terminal regions
- **Customizable**: Control how many terminal residues to add

### 5. **All (Complete Structure)**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate --missing-residues all
```
- **Behavior**: Attempts to add all missing residues
- **Best for**: Low-resolution structures requiring extensive completion
- **Warning**: May introduce many uncertain regions

## 🔧 **Advanced Options**

### **Limit Terminal Residues**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate \
    --missing-residues terminal-only --max-terminal-residues 2
```
- Controls maximum number of residues added at each terminus
- Default: 5 residues
- Range: 0-20 (higher values generate warnings)

### **Conservative Approach**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate \
    --missing-residues auto --conservative-missing
```
- Only adds missing residues with high confidence
- Limits terminal additions to `--max-terminal-residues`
- Reduces risk of adding incorrect structures

### **Skip Loop Regions**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate \
    --missing-residues all --skip-missing-loops
```
- Avoids adding residues in likely loop regions
- Uses heuristics to identify potential loops (>3 consecutive missing residues)
- Preserves well-defined secondary structure elements

### **Specify Terminal Types**
```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate \
    --missing-residues terminal-only --terminal-residue-types ACE NME NH2
```
- Controls which types of terminal residues are allowed
- Default: ACE, NME, NH2, COOH
- Useful for specific capping strategies

## 📋 **Practical Examples**

### **High-Resolution Crystal Structure**
```bash
# Conservative approach for high-quality structures
python -m EasyMD --protein high_res.pdb --steps 50000 --solvate \
    --missing-residues non-terminal --conservative-missing
```

### **NMR Structure with Flexible Termini**
```bash
# Add only a few terminal residues
python -m EasyMD --protein nmr_structure.pdb --steps 30000 --GBIS \
    --missing-residues terminal-only --max-terminal-residues 2
```

### **Low-Resolution Structure**
```bash
# Complete structure but avoid uncertain loops
python -m EasyMD --protein low_res.pdb --steps 20000 --solvate \
    --missing-residues all --skip-missing-loops --conservative-missing
```

### **Homology Model**
```bash
# Skip all missing residues to avoid model artifacts
python -m EasyMD --protein homology_model.pdb --steps 25000 --solvate \
    --missing-residues none
```

### **Membrane Protein**
```bash
# Conservative terminal addition for membrane proteins
python -m EasyMD --protein membrane_protein.pdb --steps 40000 --solvate \
    --missing-residues terminal-only --max-terminal-residues 1 \
    --terminal-residue-types ACE NME
```

## ⚠️ **Validation and Warnings**

EasyMD provides helpful warnings for missing residue settings:

### **Strategy Conflicts**
```bash
# This generates a warning:
python -m EasyMD --protein structure.pdb --steps 10000 --solvate \
    --missing-residues none --skip-missing-loops
# Warning: Loop options have no effect when missing-residues is 'none'
```

### **Excessive Terminal Residues**
```bash
# This generates a warning:
python -m EasyMD --protein structure.pdb --steps 10000 --solvate \
    --max-terminal-residues 25
# Warning: Very high max terminal residues (25). Consider if this is intended
```

### **Uncommon Terminal Types**
```bash
# This generates a warning:
python -m EasyMD --protein structure.pdb --steps 10000 --solvate \
    --terminal-residue-types CUSTOM_CAP
# Warning: Uncommon terminal residue type 'CUSTOM_CAP'
```

## 🎯 **Decision Guide**

### **Choose Based on Structure Quality:**

| Structure Type | Recommended Strategy | Reasoning |
|----------------|---------------------|-----------|
| High-res X-ray (< 2.0 Å) | `non-terminal` or `none` | Minimize artifacts |
| Medium-res X-ray (2.0-3.0 Å) | `auto` with `--conservative-missing` | Balanced approach |
| Low-res X-ray (> 3.0 Å) | `terminal-only` with low `--max-terminal-residues` | Focus on core structure |
| NMR structures | `terminal-only` | Termini often disordered |
| Homology models | `none` or `non-terminal` | Avoid model artifacts |
| AlphaFold models | `auto` or `all` | Generally high quality |

### **Choose Based on Research Goals:**

| Research Goal | Recommended Strategy | Reasoning |
|---------------|---------------------|-----------|
| Binding site analysis | `non-terminal` | Preserve active site integrity |
| Protein folding | `all` with `--conservative-missing` | Need complete structure |
| Membrane protein dynamics | `terminal-only` with minimal additions | Focus on transmembrane regions |
| Allosteric studies | `auto` | Need complete conformational network |
| Drug screening | `none` or `non-terminal` | Avoid artificial binding sites |

## 🔍 **Monitoring Missing Residue Handling**

EasyMD provides detailed output about missing residue processing:

```
Missing residue strategy: terminal-only
Chain A: keeping 3 terminal missing residues (2 N-terminal, 1 C-terminal) out of 7 total missing
Chain B: keeping all 2 missing residues (within limit)
Will add missing residues:
  Chain A: 3 residues
    <Residue 1: ACE>
    <Residue 2: GLY>
    <Residue 150: NME>
  Chain B: 2 residues
    <Residue 75: NH2>
    <Residue 76: COOH>
Missing residues will be added during addMissingAtoms() step
```

### **Verifying Terminal Residue Limits**

When using `--max-terminal-residues`, you should see output like:

```bash
python -m EasyMD --protein structure.pdb --steps 10000 --solvate \
    --missing-residues terminal-only --max-terminal-residues 3
```

**Expected Output:**
```
Missing residue strategy: terminal-only
Chain A: keeping 3 terminal missing residues (2 N-terminal, 1 C-terminal) out of 8 total missing
Chain B: keeping all 1 missing residues (within limit)
```

This confirms that:
- Chain A had 8 missing residues but only 3 are being added (respecting the limit)
- Chain B had 1 missing residue which is kept (within the limit)
- The residues are distributed between N-terminal and C-terminal ends

## 📝 **Configuration File Example**

```yaml
# missing_residue_config.yml
protein: "structure.pdb"
steps: 50000
solvate: true
temperature: 310

# Missing residue handling
missing_residues: "terminal-only"
max_terminal_residues: 3
terminal_residue_types: ["ACE", "NME", "NH2"]
conservative_missing: true
skip_missing_loops: false

# Other parameters
water_model: "tip3p"
padding: 12.0
interval: 1000
```

```bash
python -m EasyMD --config missing_residue_config.yml
```

## 🚨 **Common Pitfalls and Solutions**

### **Problem**: Simulation crashes due to missing residues
**Solution**: Use `--missing-residues auto` or `all` to complete the structure

### **Problem**: Unrealistic protein conformations
**Solution**: Use `--conservative-missing` and limit `--max-terminal-residues`

### **Problem**: Binding site artifacts
**Solution**: Use `--missing-residues non-terminal` to preserve binding regions

### **Problem**: Long equilibration times
**Solution**: Consider `--missing-residues none` if structure quality allows

## 🎓 **Best Practices**

1. **Always inspect your structure** before choosing a strategy
2. **Start conservative** and add more residues if needed
3. **Use validation tools** to check added residues
4. **Document your choices** for reproducibility
5. **Test different strategies** for critical simulations
6. **Monitor equilibration** - added residues may need longer equilibration

The missing residue handling in EasyMD gives you the control needed to prepare high-quality systems for reliable molecular dynamics simulations!