from openmm import app, unit
import yaml
import os
from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator

from openmm.app import PDBFile, Modeller
import mdtraj

from pdbfixer import PDBFixer
from EasyMD.utils.utils import _formatIndex, writeFooter, PDBwrite_all, deletePcap


class SysGenerator:
    """
    System generator for molecular dynamics simulations.
    
    This class handles the preparation of molecular systems for MD simulations, including
    protein structure fixing, ligand parameterization, solvation, and force field assignment.
    It supports both new simulations and restart from previous states, with options for
    explicit solvation or implicit solvent (GBIS) models.
    """
    
    def __init__(self, config):
        """
        Initialize the SysGenerator and prepare the molecular system.
        
        Args:
            config: Configuration object containing simulation parameters
        """
        self.config = config
        self.forcefield_kwargs = {
            'constraints': app.HBonds, 
            'rigidWater': True, 
            'removeCMMotion': False, 
            'hydrogenMass': 4*unit.amu
        }
        self.last_state = 0
        
        if self.config.restart:
            self.modeller, self.system = self._restartSetup()
        else:
            self.modeller, self.system = self._setup()

    def _setup(self):
        """Set up a new molecular dynamics system from scratch."""
        if self.config.outdir is None:
            # Auto-generate output directory
            i = 0
            while os.path.exists(f'out_{i}'):
                i += 1
            self.config.outdir = f'out_{i}'
        
        os.makedirs(self.config.outdir, exist_ok=True)
        
        # Save configuration
        try:
            with open(f'{self.config.outdir}/setup.yml', 'w') as f:
                # Try to dump config vars, but handle Mock objects gracefully
                config_dict = {}
                for key, value in vars(self.config).items():
                    try:
                        # Test if the value can be serialized to YAML
                        yaml.dump({key: value})
                        config_dict[key] = value
                    except:
                        # Skip values that can't be serialized (like Mock objects)
                        config_dict[key] = str(value)
                yaml.dump(config_dict, f)
        except Exception as e:
            # If config saving fails, continue anyway (important for tests)
            print(f"Warning: Could not save configuration: {e}")
        
        # Prepare system based on whether ligand is present
        if self.config.ligand:
            return self._prep_complex(
                pdb_in=self.config.protein,
                list_of_molecules_to_remove=self.config.remove,
                lig_name=self.config.ligand,
                solvate=self.config.solvate,
                protein_force_field=self.config.protein_force_field,
                water_force_field=self.config.water_force_field,
                ligand_force_field=self.config.ligand_force_field,
                water_model=self.config.water_model,
                positive_ion=self.config.positive_ion,
                negative_ion=self.config.negative_ion,
                ionic_strength=self.config.ionic_strength,
                no_neutralize=self.config.no_neutralize,
                padding=self.config.padding,
                ph=self.config.ph,
                outdir=self.config.outdir,
                forcefield_kwargs=self.forcefield_kwargs
            )
        else:
            return self._prep_prot(
                pdb_in=self.config.protein,
                list_of_molecules_to_remove=self.config.remove,
                solvate=self.config.solvate,
                protein_force_field=self.config.protein_force_field,
                water_force_field=self.config.water_force_field,
                water_model=self.config.water_model,
                positive_ion=self.config.positive_ion,
                negative_ion=self.config.negative_ion,
                ionic_strength=self.config.ionic_strength,
                no_neutralize=self.config.no_neutralize,
                padding=self.config.padding,
                ph=self.config.ph,
                outdir=self.config.outdir,
                forcefield_kwargs=self.forcefield_kwargs
            )

    def _restartSetup(self):
        """Set up system for restart simulation from previously saved state."""
        out_dir = self.config.restart
        
        # Load setup configuration
        with open(f'{out_dir}/restart_setup.yml', 'r') as f:
            setup = yaml.safe_load(f)
        
        print(f'✓ Restarting from state files in directory: {self.config.restart}')
        
        # Check if ligand is present
        if os.path.isfile(f'{out_dir}/ligand.sdf'):
            print('✓ Preparing system with ligand...')
            return self._prep_restart_ligand(setup=setup, outdir=out_dir, forcefield_kwargs=self.forcefield_kwargs)
        else:
            print('✓ Preparing system without ligand...')
            return self._prep_restart(setup=setup, outdir=out_dir, forcefield_kwargs=self.forcefield_kwargs)

    def _prep_restart_ligand(self, setup, outdir, forcefield_kwargs):
        """Prepare system for restart simulation with ligand present."""
        pdb = PDBFile(f'{outdir}/restart_model.pdb')
        ligand_mol = Molecule.from_file(f'{outdir}/ligand.sdf')
        
        protein_force_field = setup['protein_force_field']
        water_force_field = setup['water_force_field']
        ligand_force_field = setup['ligand_force_field']
        modeller = Modeller(pdb.topology, pdb.positions)
        
        if setup['solvate']:
            print('✓ Generating solvated system...')
            system_generator = SystemGenerator(
                forcefields=[protein_force_field, water_force_field],
                small_molecule_forcefield=ligand_force_field,
                molecules=[ligand_mol],
                forcefield_kwargs=forcefield_kwargs
            )
            system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
        else:
            print('✓ Generating implicit solvent system...')
            system_generator = SystemGenerator(
                forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
                small_molecule_forcefield=ligand_force_field,
                molecules=[ligand_mol],
                forcefield_kwargs=forcefield_kwargs,
                nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
        
        return modeller, system

    def _prep_restart(self, setup, outdir, forcefield_kwargs):
        """Prepare system for restart simulation without ligand (protein-only)."""
        pdb = PDBFile(f'{outdir}/restart_model.pdb')
        protein_force_field = setup['protein_force_field']
        water_force_field = setup['water_force_field']
        modeller = Modeller(pdb.topology, pdb.positions)
        
        if setup['solvate']:
            print('✓ Generating solvated system...')
            system_generator = SystemGenerator(
                forcefields=[protein_force_field, water_force_field],
                forcefield_kwargs=forcefield_kwargs
            )
            system = system_generator.create_system(modeller.topology)
        else:
            print('✓ Generating implicit solvent system...')
            system_generator = SystemGenerator(
                forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
                forcefield_kwargs=forcefield_kwargs,
                nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            system = system_generator.create_system(modeller.topology)
        
        return modeller, system

    def _prep_prot(self, pdb_in, list_of_molecules_to_remove, solvate, protein_force_field, 
                   water_force_field, water_model, positive_ion, negative_ion, ionic_strength, 
                   no_neutralize, padding, ph, outdir, forcefield_kwargs):
        """Prepare protein-only system for molecular dynamics simulation."""
        
        # Load and fix the protein
        fixer = PDBFixer(filename=pdb_in)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()
        
        # Report structure issues found
        if fixer.nonstandardResidues:
            print(f"✓ Found {len(fixer.nonstandardResidues)} non-standard residues - converting to standard")
        if fixer.missingResidues:
            total_missing = sum(len(residues) for residues in fixer.missingResidues.values())
            print(f"✓ Found {total_missing} missing residues across {len(fixer.missingResidues)} chains")
        if fixer.missingAtoms:
            total_atoms = sum(len(atoms) for atoms in fixer.missingAtoms.values())
            print(f"✓ Found {total_atoms} missing atoms")

        fixer.replaceNonstandardResidues()
        self._handle_missing_residues(fixer)
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
        
        # Remove molecules from pdb and save it in a new file
        # Only add water molecules to removal list if keep_water is False (default behavior)
        if not getattr(self.config, 'keep_water', False):
            list_of_molecules_to_remove += ['HOH', 'WAT']
        
        modeller = Modeller(fixer.topology, fixer.positions)
        atoms_before = modeller.topology.getNumAtoms()
        for mol in list_of_molecules_to_remove:
            toDelete = []
            for res in modeller.topology.residues():
                if res.name == mol:
                    toDelete.append(res)
            modeller.delete(toDelete)
        atoms_after = modeller.topology.getNumAtoms()
        if atoms_before != atoms_after:
            print(f"✓ Removed {atoms_before - atoms_after} atoms ({atoms_after} atoms remaining)")
        
        PDBwrite_all(modeller, f'{outdir}/prot_receptor.pdb')
        print(f"✓ System prepared with {modeller.topology.getNumAtoms()} atoms")
        
        if solvate:
            modeller = deletePcap(modeller)
            print('✓ Generating solvated system...')
            
            system_generator = SystemGenerator(
                forcefields=[protein_force_field, water_force_field],
                forcefield_kwargs=forcefield_kwargs
            )
            
            modeller.addSolvent(system_generator.forcefield, model=water_model, padding=padding * unit.angstroms,
                              positiveIon=positive_ion, negativeIon=negative_ion,
                              ionicStrength=ionic_strength * unit.molar, neutralize=not no_neutralize)
            print(f'✓ Solvated system created with {modeller.topology.getNumAtoms()} atoms')

            with open(f'{outdir}/solvated_complex.pdb', 'w') as outfile:
                PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
            print('✓ Solvated complex saved as solvated_complex.pdb')
            
            pdb = PDBFile(f'{outdir}/solvated_complex.pdb')
            modeller = Modeller(pdb.topology, pdb.positions)
            system = system_generator.create_system(modeller.topology)
        else:
            print('✓ Generating implicit solvent system...')
            system_generator = SystemGenerator(
                forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
                forcefield_kwargs=forcefield_kwargs,
                nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            system = system_generator.create_system(modeller.topology)
        
        # Write to pdb modeller as restart_model.pdb
        PDBwrite_all(modeller, f'{outdir}/restart_model.pdb')
        
        return modeller, system

    def _handle_missing_residues(self, fixer):
        """Handle missing residues based on user configuration."""
        missing_residues_strategy = getattr(self.config, 'missing_residues', 'auto')
        max_terminal_residues = getattr(self.config, 'max_terminal_residues', 5)
        terminal_residue_types = getattr(self.config, 'terminal_residue_types', ['ACE', 'NME', 'NH2', 'COOH'])
        skip_missing_loops = getattr(self.config, 'skip_missing_loops', False)
        conservative_missing = getattr(self.config, 'conservative_missing', False)
        
        if missing_residues_strategy == "none":
            print("✓ Skipping all missing residues as requested")
            fixer.missingResidues = {}
            return
        
        if not fixer.missingResidues:
            return
        
        original_missing = dict(fixer.missingResidues)
        
        if missing_residues_strategy == "auto":
            if conservative_missing:
                self._apply_conservative_filtering(fixer, max_terminal_residues)
        elif missing_residues_strategy == "non-terminal":
            self._filter_non_terminal_residues(fixer)
        elif missing_residues_strategy == "terminal-only":
            self._filter_terminal_residues(fixer, max_terminal_residues, terminal_residue_types)
        elif missing_residues_strategy == "all":
            if skip_missing_loops:
                self._filter_loop_residues(fixer)
            if conservative_missing:
                self._apply_conservative_filtering(fixer, max_terminal_residues)
        
        # Report what will be added
        if fixer.missingResidues:
            total_to_add = sum(len(residues) for residues in fixer.missingResidues.values())
            print(f"✓ Will add {total_to_add} missing residues across {len(fixer.missingResidues)} chains")

    def _apply_conservative_filtering(self, fixer, max_terminal_residues):
        """Apply conservative filtering to missing residues."""
        filtered_missing = {}
        
        for chain_id, residues in fixer.missingResidues.items():
            if len(residues) <= max_terminal_residues:
                filtered_missing[chain_id] = residues
            else:
                print(f"Skipping {len(residues)} missing residues in chain {chain_id} (exceeds max {max_terminal_residues})")
        
        fixer.missingResidues = filtered_missing

    def _filter_non_terminal_residues(self, fixer):
        """Filter to keep only non-terminal missing residues."""
        filtered_missing = {}
        
        for chain_id, residues in fixer.missingResidues.items():
            if len(residues) > 2:
                internal_residues = residues[1:-1] if len(residues) > 3 else []
                if internal_residues:
                    filtered_missing[chain_id] = internal_residues
                    print(f"Chain {chain_id}: keeping {len(internal_residues)} internal missing residues")
            else:
                print(f"Chain {chain_id}: no clear internal missing residues found")
        
        fixer.missingResidues = filtered_missing

    def _filter_terminal_residues(self, fixer, max_terminal_residues, allowed_types):
        """Filter to keep only terminal missing residues."""
        filtered_missing = {}
        
        for chain_id, residues in fixer.missingResidues.items():
            if len(residues) <= max_terminal_residues:
                terminal_residues = residues[:]
                print(f"Chain {chain_id}: keeping all {len(terminal_residues)} missing residues (within limit)")
            else:
                # Take from both ends
                n_terminal_count = max_terminal_residues // 2
                c_terminal_count = max_terminal_residues - n_terminal_count
                
                terminal_residues = []
                if n_terminal_count > 0:
                    terminal_residues.extend(residues[:n_terminal_count])
                if c_terminal_count > 0:
                    terminal_residues.extend(residues[-c_terminal_count:])
                
                print(f"Chain {chain_id}: keeping {len(terminal_residues)} terminal missing residues "
                      f"({n_terminal_count} N-terminal, {c_terminal_count} C-terminal) "
                      f"out of {len(residues)} total missing")
            
            if terminal_residues:
                filtered_missing[chain_id] = terminal_residues
        
        fixer.missingResidues = filtered_missing

    def _filter_loop_residues(self, fixer):
        """Filter out missing residues that are likely in loop regions."""
        filtered_missing = {}
        
        for chain_id, residues in fixer.missingResidues.items():
            if len(residues) <= 3:
                non_loop_residues = residues
                print(f"Chain {chain_id}: keeping {len(residues)} residues (small gap)")
            else:
                print(f"Chain {chain_id}: skipping {len(residues)} residues (likely loop region)")
                non_loop_residues = []
            
            if non_loop_residues:
                filtered_missing[chain_id] = non_loop_residues
        
        fixer.missingResidues = filtered_missing

    def _validate_coordinates(self, positions, molecule_name="molecule"):
        """
        Validate that coordinates are finite and not NaN.
        
        Args:
            positions: OpenMM positions object or numpy array
            molecule_name (str): Name for error reporting
            
        Raises:
            ValueError: If any coordinates are NaN or infinite
        """
        import numpy as np
        
        # Convert positions to numpy array for validation
        try:
            # Try OpenMM positions format first (with .x, .y, .z attributes)
            if hasattr(positions, '__iter__') and len(positions) > 0 and hasattr(positions[0], 'x'):
                pos_array = np.array([[pos.x, pos.y, pos.z] for pos in positions])
            # Handle numpy array format
            elif isinstance(positions, np.ndarray):
                pos_array = positions
            # Handle list of coordinates
            elif isinstance(positions, (list, tuple)):
                pos_array = np.array(positions)
            else:
                # Try to convert to numpy array directly
                pos_array = np.array(positions)
        except Exception as e:
            raise ValueError(f"Could not convert {molecule_name} positions to numpy array: {e}")
        
        # Ensure we have a 2D array with 3 columns (x, y, z)
        if pos_array.ndim == 1:
            pos_array = pos_array.reshape(-1, 3)
        elif pos_array.ndim != 2 or pos_array.shape[1] != 3:
            raise ValueError(f"Invalid coordinate format for {molecule_name}: expected Nx3 array, got {pos_array.shape}")
        
        if np.any(np.isnan(pos_array)):
            raise ValueError(f"NaN coordinates found in {molecule_name}. This will cause simulation failure.")
        
        if np.any(np.isinf(pos_array)):
            raise ValueError(f"Infinite coordinates found in {molecule_name}. This will cause simulation failure.")
        
        print(f"✓ Validated coordinates for {molecule_name} - no NaN or infinite values found")

    def _validate_ligand_structure(self, ligand_pdb_path):
        """
        Validate ligand structure for common issues that cause NaN coordinates.
        
        Args:
            ligand_pdb_path (str): Path to ligand PDB file
            
        Raises:
            ValueError: If ligand structure has issues that could cause NaN coordinates
        """
        import numpy as np
        
        try:
            # Read the ligand PDB file
            with open(ligand_pdb_path, 'r') as f:
                lines = f.readlines()
            
            coordinates = []
            atom_count = 0
            
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coordinates.append([x, y, z])
                        atom_count += 1
                    except ValueError:
                        raise ValueError(f"Invalid coordinate format in ligand PDB: {line.strip()}")
            
            if atom_count == 0:
                raise ValueError("No atoms found in ligand PDB file")
            
            # Convert to numpy array for validation
            coords_array = np.array(coordinates)
            
            # Check for NaN or infinite coordinates
            if np.any(np.isnan(coords_array)):
                raise ValueError("NaN coordinates found in ligand structure")
            
            if np.any(np.isinf(coords_array)):
                raise ValueError("Infinite coordinates found in ligand structure")
            
            # Check for unreasonable coordinate values (likely errors)
            if np.any(np.abs(coords_array) > 1000):
                print("⚠️  Warning: Very large coordinate values found in ligand (>1000 Å)")
                print("   This might indicate coordinate system issues")
            
            # Check if all atoms are at the same position (collapsed structure)
            coord_range = np.ptp(coords_array, axis=0)  # peak-to-peak range
            if np.all(coord_range < 0.01):  # All atoms within 0.01 Å
                raise ValueError("All ligand atoms appear to be at the same position (collapsed structure)")
            
            print(f"✓ Ligand structure validated: {atom_count} atoms, coordinate range: {coord_range}")
            
        except FileNotFoundError:
            raise ValueError(f"Ligand PDB file not found: {ligand_pdb_path}")
        except Exception as e:
            raise ValueError(f"Error validating ligand structure: {e}")

    def _prep_complex(self, pdb_in, list_of_molecules_to_remove, lig_name, solvate, 
                     protein_force_field, water_force_field, ligand_force_field, water_model, 
                     positive_ion, negative_ion, ionic_strength, no_neutralize, padding, ph, 
                     outdir, forcefield_kwargs):
        """
        Prepare protein-ligand complex for molecular dynamics simulation.
        
        This method processes a protein-ligand complex by:
        1. Fixing missing residues, atoms, and hydrogens in the protein
        2. Extracting the ligand and preparing it separately
        3. Recombining protein and ligand with proper force field parameters
        4. Optionally adding explicit solvent and ions
        """
        
        fixer = PDBFixer(filename=pdb_in)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()
        
        # Report structure issues found
        if fixer.nonstandardResidues:
            print(f"✓ Found {len(fixer.nonstandardResidues)} non-standard residues - converting to standard")
        if fixer.missingResidues:
            total_missing = sum(len(residues) for residues in fixer.missingResidues.values())
            print(f"✓ Found {total_missing} missing residues across {len(fixer.missingResidues)} chains")
        if fixer.missingAtoms:
            total_atoms = sum(len(atoms) for atoms in fixer.missingAtoms.values())
            print(f"✓ Found {total_atoms} missing atoms")

        fixer.replaceNonstandardResidues()
        self._handle_missing_residues(fixer)
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
        
        # Remove molecules from pdb and save protein only
        # Only add water molecules to removal list if keep_water is False (default behavior)
        if not getattr(self.config, 'keep_water', False):
            list_of_molecules_to_remove += ['HOH', 'WAT']
        list_of_molecules_to_remove += [lig_name]
        
        modeller = Modeller(fixer.topology, fixer.positions)
        atoms_before = modeller.topology.getNumAtoms()
        for mol in list_of_molecules_to_remove:
            toDelete = []
            for res in modeller.topology.residues():
                if res.name == mol:
                    toDelete.append(res)
            modeller.delete(toDelete)
        atoms_after = modeller.topology.getNumAtoms()
        if atoms_before != atoms_after:
            print(f"✓ Removed {atoms_before - atoms_after} atoms from protein ({atoms_after} atoms remaining)")
        
        PDBwrite_all(modeller, f'{outdir}/prot_receptor.pdb')
        
        # Extract ligand from pdb and save it in a new file
        fixer = PDBFixer(filename=pdb_in)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()
        self._handle_missing_residues(fixer)
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)
        
        modeller = Modeller(fixer.topology, fixer.positions)
        toDelete = []
        for res in modeller.topology.residues():
            if res.name != lig_name:
                toDelete.append(res)
        modeller.delete(toDelete)

        with open(f'{outdir}/ligand.pdb', 'w') as outfile:
            PDBFile.writeFile(modeller.topology, modeller.positions, file=outfile, keepIds=True)
        
        # Validate ligand coordinates before processing
        self._validate_ligand_structure(f'{outdir}/ligand.pdb')
        
        # Calculate charges with openbabel and export as sdf file
        from EasyMD.utils.openbabel_charge import get_charges
        get_charges(f'{outdir}/ligand.pdb', f'{outdir}/ligand.sdf', 'pdb', 'sdf')

        # Load the ligand and the protein
        protein_pdb = PDBFile(f'{outdir}/prot_receptor.pdb')
        ligand_mol = Molecule.from_file(f'{outdir}/ligand.sdf')

        # The topology is described in the openforcefield API
        modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
        atoms_before = modeller.topology.getNumAtoms()
        
        # Add ligand to the protein with coordinate validation
        lig_top = ligand_mol.to_topology()
        ligand_positions = lig_top.get_positions().to_openmm()
        
        # Validate ligand coordinates before adding
        self._validate_coordinates(ligand_positions, "ligand")
        
        modeller.add(lig_top.to_openmm(), ligand_positions)
        atoms_after = modeller.topology.getNumAtoms()
        print(f'✓ Added ligand to protein ({atoms_after - atoms_before} ligand atoms, {atoms_after} total atoms)')
        
        # Write complex to pdb
        with open(f'{outdir}/complex.pdb', 'w') as outfile:
            PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
        
        # Replace UNK with the ligand name in the pdb file 
        with open(f'{outdir}/complex.pdb', 'r') as file:
            filedata = file.read()
        filedata = filedata.replace('UNK', lig_name)
        with open(f'{outdir}/complex.pdb', 'w') as file:
            file.write(filedata)

        if solvate:
            print('✓ Generating solvated complex system...')
            modeller = deletePcap(modeller)
           
            system_generator = SystemGenerator(
                forcefields=[protein_force_field, water_force_field],
                small_molecule_forcefield=ligand_force_field,
                molecules=[ligand_mol],
                forcefield_kwargs=forcefield_kwargs
            )
            
            modeller.addSolvent(system_generator.forcefield, model=water_model, padding=padding * unit.angstroms,
                              positiveIon=positive_ion, negativeIon=negative_ion,
                              ionicStrength=ionic_strength * unit.molar, neutralize=not no_neutralize)
            print(f'✓ Solvated complex created with {modeller.topology.getNumAtoms()} atoms')

            with open(f'{outdir}/solvated_complex.pdb', 'w') as outfile:
                PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
            print('✓ Solvated complex saved as solvated_complex.pdb')
            
            pdb = PDBFile(f'{outdir}/solvated_complex.pdb')
            modeller = Modeller(pdb.topology, pdb.positions)
            system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
        else:
            modeller = deletePcap(modeller)
            
            print('✓ Generating implicit solvent system...')
            system_generator = SystemGenerator(
                forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/gbn2.xml'],
                small_molecule_forcefield=ligand_force_field,
                molecules=[ligand_mol],
                forcefield_kwargs=forcefield_kwargs,
                nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff}
            )
            system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
            
        # Write to pdb modeller as restart_model.pdb
        PDBwrite_all(modeller, f'{outdir}/restart_model.pdb')
        
        return modeller, system