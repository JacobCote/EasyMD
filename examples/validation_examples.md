# EasyMD Input Validation Examples

This document demonstrates the enhanced input validation system in EasyMD, showing various scenarios and the helpful error messages provided.

## ✅ Valid Command Examples

### Basic Protein Simulation (Explicit Solvent)
```bash
python -m EasyMD --protein 4zgm.pdb --steps 10000 --solvate --temperature 300
```

### Protein-Ligand Complex (Implicit Solvent)
```bash
python -m EasyMD --protein complex.pdb --ligand LIG --steps 50000 --GBIS --temperature 310
```

### Time-based Simulation
```bash
python -m EasyMD --protein protein.pdb --clock 120 --solvate --water-model tip4pew
```

### Using Configuration File
```bash
python -m EasyMD --config simulation.yml
```

### Restart Simulation
```bash
python -m EasyMD --restart out_5 --steps 10000
```

## ❌ Common Validation Errors and Fixes

### 1. Missing Protein File
**Command:**
```bash
python -m EasyMD --steps 1000 --solvate
```

**Error Message:**
```
❌ INPUT VALIDATION FAILED
════════════════════════════════════════════════════════════

1. ❌ Protein PDB file is required. Use --protein <file.pdb>

QUICK FIXES:
• Add protein file: --protein your_protein.pdb
```

**Fix:**
```bash
python -m EasyMD --protein 4zgm.pdb --steps 1000 --solvate
```

### 2. Conflicting Duration Methods
**Command:**
```bash
python -m EasyMD --protein 4zgm.pdb --steps 1000 --clock 60 --solvate
```

**Error Message:**
```
❌ INPUT VALIDATION FAILED
════════════════════════════════════════════════════════════

1. ❌ Cannot specify both --steps and --clock. Choose one simulation duration method

QUICK FIXES:
• Use either --steps 10000 OR --clock 60 (not both)
```

**Fix:**
```bash
# Option 1: Use steps
python -m EasyMD --protein 4zgm.pdb --steps 1000 --solvate

# Option 2: Use time
python -m EasyMD --protein 4zgm.pdb --clock 60 --solvate
```

### 3. No Solvation Method Specified
**Command:**
```bash
python -m EasyMD --protein 4zgm.pdb --steps 1000
```

**Error Message:**
```
❌ INPUT VALIDATION FAILED
════════════════════════════════════════════════════════════

1. ❌ Must choose exactly one solvation method: --solvate OR --GBIS

QUICK FIXES:
• Choose solvation: --solvate (explicit) OR --GBIS (implicit)
```

**Fix:**
```bash
# Explicit solvent (slower, more accurate)
python -m EasyMD --protein 4zgm.pdb --steps 1000 --solvate

# Implicit solvent (faster, less accurate)
python -m EasyMD --protein 4zgm.pdb --steps 1000 --GBIS
```

### 4. Invalid Physical Parameters
**Command:**
```bash
python -m EasyMD --protein 4zgm.pdb --steps 1000 --solvate --temperature -100
```

**Error Message:**
```
❌ INPUT VALIDATION FAILED
════════════════════════════════════════════════════════════

1. ❌ Temperature must be positive, got -100 K
2. ⚠️  Temperature -100 K is outside typical range (250-400 K)
```

**Fix:**
```bash
python -m EasyMD --protein 4zgm.pdb --steps 1000 --solvate --temperature 300
```

### 5. File Not Found
**Command:**
```bash
python -m EasyMD --protein nonexistent.pdb --steps 1000 --solvate
```

**Error Message:**
```
❌ INPUT VALIDATION FAILED
════════════════════════════════════════════════════════════

1. ❌ Protein file 'nonexistent.pdb' does not exist or is not accessible

QUICK FIXES:
• Check protein file path and ensure file exists
```

**Fix:**
```bash
# Make sure the file exists
ls -la *.pdb
python -m EasyMD --protein existing_file.pdb --steps 1000 --solvate
```

## ⚠️ Common Warnings

### 1. Very Short Simulation
**Command:**
```bash
python -m EasyMD --protein 4zgm.pdb --steps 100 --solvate
```

**Warning:**
```
⚠️  1. Warning: Very short simulation (100 steps). Consider at least 10,000 steps for meaningful results
```

### 2. Unusual Force Field
**Command:**
```bash
python -m EasyMD --protein 4zgm.pdb --steps 1000 --solvate --protein-force-field custom.xml
```

**Warning:**
```
⚠️  1. Warning: Uncommon protein force field 'custom.xml'. Common options: amber14-all.xml, amber99sb-ildn.xml, amber03.xml
```

### 3. Long Ligand Name
**Command:**
```bash
python -m EasyMD --protein complex.pdb --ligand VERYLONGNAME --steps 1000 --GBIS
```

**Warning:**
```
⚠️  1. Warning: Ligand name 'VERYLONGNAME' is longer than 4 characters. PDB format typically uses 3-letter codes
```

## ✅ Successful Validation

When all parameters are valid, you'll see:

```
════════════════════════════════════════════════════════════
✅ INPUT VALIDATION SUCCESSFUL
════════════════════════════════════════════════════════════

Simulation Configuration Summary:
----------------------------------------
Mode: Standard MD Simulation
Protein: 4zgm.pdb
Ligand: LIG
Solvation: Explicit (tip3p, padding=10Å)
Duration: 10,000 steps
Temperature: 300 K
Output: auto-generated directory

════════════════════════════════════════════════════════════
Ready to start simulation!
════════════════════════════════════════════════════════════
```

## Configuration File Examples

### Valid Configuration (simulation.yml)
```yaml
# EasyMD Simulation Configuration
protein: "4zgm.pdb"
ligand: "LIG"
steps: 50000
temperature: 300
solvate: true
water_model: "tip3p"
padding: 12.0
ionic_strength: 0.15
interval: 1000
equilibration_steps: 1000
protein_force_field: "amber14-all.xml"
ligand_force_field: "openff-2.2.0"
water_force_field: "amber/tip3p_standard.xml"
outdir: "my_simulation"
remove: ["DMS", "SO4"]
ph: 7.4
```

### Usage:
```bash
python -m EasyMD --config simulation.yml
```

## Best Practices

1. **Always specify a protein file** for new simulations
2. **Choose one solvation method**: `--solvate` for explicit or `--GBIS` for implicit
3. **Set appropriate simulation length**: At least 10,000 steps for meaningful results
4. **Use reasonable temperatures**: 250-400 K for most biological systems
5. **Check file paths**: Ensure all input files exist and are accessible
6. **Use standard force fields**: Stick to well-tested options unless you have specific needs
7. **Validate before long runs**: Test with short simulations first

## Troubleshooting

If you encounter validation errors:

1. **Read the error messages carefully** - they provide specific guidance
2. **Check the QUICK FIXES section** for immediate solutions
3. **Verify file paths** and ensure files exist
4. **Use `--help`** to see all available options
5. **Start with simple configurations** and add complexity gradually
6. **Test with short simulations** before running long productions

## Getting Help

- Use `python -m EasyMD --help` for detailed option descriptions
- Check example configurations in the `examples/` directory
- Review validation messages for specific guidance
- Consult the documentation for advanced usage patterns