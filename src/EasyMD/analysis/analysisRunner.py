"""
Analysis Runner for EasyMD

This module performs trajectory analysis including RMSD, RMSF, distances,
and other structural metrics.
"""

import os
import sys
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import mdtraj as md
from typing import List, Optional, Tuple, Dict, Any


class AnalysisRunner:
    """
    Performs trajectory analysis for molecular dynamics simulations.
    
    This class handles loading trajectory data, performing various analyses
    (RMSD, RMSF, distances, etc.), and generating output plots and data files.
    
    Attributes:
        config: Configuration object containing analysis parameters
        trajectory_dir (Path): Path to trajectory directory
        output_dir (Path): Path to output directory
        topology: MDTraj topology object
        trajectory: Combined MDTraj trajectory object
    """
    
    def __init__(self, config):
        """
        Initialize the AnalysisRunner with configuration.
        
        Args:
            config: Configuration object from AnalysisManager
        """
        self.config = config
        self.trajectory_dir = Path(config.trajectory_directory)
        self.output_dir = Path(config.output_dir)
        self.topology = None
        self.trajectory = None
        
        # Set matplotlib style
        if config.plot_style != "default":
            try:
                mplstyle.use(config.plot_style)
            except OSError:
                print(f"Warning: Style '{config.plot_style}' not available, using default")
    
    def run(self):
        """
        Execute the complete analysis workflow.
        
        This method orchestrates the entire analysis process:
        1. Load trajectory data
        2. Perform selected analyses
        3. Generate plots and save data
        4. Display completion summary
        """
        try:
            print("="*60)
            print("🔬 STARTING TRAJECTORY ANALYSIS")
            print("="*60)
            
            # Load trajectory data
            self._load_trajectory_data()
            
            # Perform analyses
            results = {}
            if self.config.all:
                results.update(self._perform_all_analyses())
            else:
                if self.config.rmsd:
                    results['rmsd'] = self._calculate_rmsd()
                if self.config.rmsf:
                    results['rmsf'] = self._calculate_rmsf()
                if self.config.distances:
                    results['distances'] = self._calculate_distances()
                if self.config.radius_gyration:
                    results['radius_gyration'] = self._calculate_radius_gyration()
                if self.config.secondary_structure:
                    results['secondary_structure'] = self._analyze_secondary_structure()
            
            # Generate outputs
            
            self._generate_outputs(results)
            
            
            # Display completion summary
            self._display_completion_summary(results)
            
            
        except Exception as e:
            print(f"\n❌ Analysis failed with error: {str(e)}")
            print("Check the error message above and ensure all input files are valid.")
            
            sys.exit(1)
    
    def _load_trajectory_data(self):
        """
        Load trajectory files and topology.
        
        This method loads the topology from the pickle file and combines
        all trajectory files into a single MDTraj trajectory object.
        """
        print("📁 Loading trajectory data...")
        
        # Load topology
        topology_file = self.trajectory_dir / 'topology.pkl'
        print(f"   Loading topology from {topology_file}")
        
        try:
            with open(topology_file, 'rb') as f:
                topo_openmm = pickle.load(f)
            self.topology = md.Topology.from_openmm(topo_openmm)
            print(f"   ✅ Topology loaded: {self.topology.n_atoms} atoms, {self.topology.n_residues} residues")
        except Exception as e:
            raise RuntimeError(f"Failed to load topology: {str(e)}")
        
        # Find and load trajectory files
        dcd_files = sorted([f for f in os.listdir(self.trajectory_dir) if f.startswith('output_traj_') and f.endswith('.dcd')])
        
        if not dcd_files:
            raise RuntimeError("No trajectory files found")
        
        print(f"   Found {len(dcd_files)} trajectory files")
        print("   Loading trajectories...")
        
        try:
            trajectories = []
            for i, dcd_file in enumerate(dcd_files):
                traj_path = self.trajectory_dir / dcd_file
                traj = md.load(str(traj_path), top=self.topology)
                trajectories.append(traj)
                print(f"   ✅ Loaded {dcd_file}: {traj.n_frames} frames")
            
            # Combine trajectories
            self.trajectory = md.join(trajectories)
            print(f"   ✅ Combined trajectory: {self.trajectory.n_frames} total frames")
            
            # Apply frame selection if specified
            self._apply_frame_selection()
            
        except Exception as e:
            raise RuntimeError(f"Failed to load trajectories: {str(e)}")
    
    def _apply_frame_selection(self):
        """Apply frame selection based on start, end, and skip parameters."""
        start = self.config.start_frame
        end = self.config.end_frame if self.config.end_frame is not None else self.trajectory.n_frames
        skip = self.config.skip_frames
        
        if start > 0 or end < self.trajectory.n_frames or skip > 1:
            original_frames = self.trajectory.n_frames
            frame_indices = list(range(start, min(end, self.trajectory.n_frames), skip))
            self.trajectory = self.trajectory[frame_indices]
            print(f"   📊 Frame selection applied: {original_frames} → {self.trajectory.n_frames} frames")
    
    def _perform_all_analyses(self) -> Dict[str, Any]:
        """Perform all available analyses."""
        print("\n🔬 Performing all analyses...")
        results = {}
        
        results['rmsd'] = self._calculate_rmsd()
        results['rmsf'] = self._calculate_rmsf()
        results['distances'] = self._calculate_distances()
        results['radius_gyration'] = self._calculate_radius_gyration()
        results['secondary_structure'] = self._analyze_secondary_structure()
        
        return results
    
    def _calculate_rmsd(self) -> Dict[str, Any]:
        """
        Calculate Root Mean Square Deviation (RMSD) for each protein chain separately.
        
        Returns:
            Dictionary containing RMSD data and metadata for each chain
        """
        print("   📈 Calculating RMSD (per chain, excluding water)...")
        
        try:
            reference_frame = min(self.config.reference_frame, self.trajectory.n_frames - 1)
            chains = self._get_protein_chains()
            
            if not chains:
                print("   ⚠️  No protein chains found")
                return None
            
            rmsd_results = {
                'per_chain': {},
                'combined': None,
                'reference_frame': reference_frame,
                'atom_selection': self.config.atom_selection,
                'n_chains': len(chains)
            }
            
            # Calculate RMSD for each chain separately
            all_chain_indices = []
            for chain_id, chain_atom_indices in chains.items():
                if not chain_atom_indices:
                    continue
                
                # Get selection-specific indices for this chain
                selection_indices = self._get_chain_atom_indices(chain_id, self.config.atom_selection)
                
                if not selection_indices:
                    print(f"   ⚠️  No {self.config.atom_selection} atoms found in chain {chain_id}")
                    continue
                
                # Calculate RMSD for this chain
                chain_rmsd = md.rmsd(self.trajectory, self.trajectory, reference_frame, atom_indices=selection_indices)
                chain_rmsd *= 10  # Convert nm to Angstroms
                
                # Get chain name/identifier
                chain_name = f"Chain_{chain_id}"
                try:
                    # Try to get actual chain ID if available
                    first_atom = next(atom for atom in self.trajectory.topology.atoms if atom.index in selection_indices)
                    if hasattr(first_atom.residue.chain, 'id') and first_atom.residue.chain.id:
                        chain_name = f"Chain_{first_atom.residue.chain.id}"
                except:
                    pass
                
                rmsd_results['per_chain'][chain_name] = {
                    'values': chain_rmsd,
                    'frames': np.arange(len(chain_rmsd)),
                    'n_atoms': len(selection_indices),
                    'mean': np.mean(chain_rmsd),
                    'std': np.std(chain_rmsd),
                    'max': np.max(chain_rmsd),
                    'min': np.min(chain_rmsd),
                    'chain_id': chain_id
                }
                
                all_chain_indices.extend(selection_indices)
                
                print(f"   ✅ {chain_name}: {len(selection_indices)} atoms, Mean RMSD: {np.mean(chain_rmsd):.2f} Å")
            
            # Calculate combined RMSD for all protein chains
            if all_chain_indices:
                combined_rmsd = md.rmsd(self.trajectory, self.trajectory, reference_frame, atom_indices=all_chain_indices)
                combined_rmsd *= 10  # Convert nm to Angstroms
                
                rmsd_results['combined'] = {
                    'values': combined_rmsd,
                    'frames': np.arange(len(combined_rmsd)),
                    'n_atoms': len(all_chain_indices),
                    'mean': np.mean(combined_rmsd),
                    'std': np.std(combined_rmsd),
                    'max': np.max(combined_rmsd),
                    'min': np.min(combined_rmsd)
                }
                
                print(f"   ✅ Combined: {len(all_chain_indices)} atoms, Mean RMSD: {np.mean(combined_rmsd):.2f} Å")
            
            return rmsd_results
            
        except Exception as e:
            print(f"   ❌ RMSD calculation failed: {str(e)}")
            return None
    
    def _calculate_rmsf(self) -> Dict[str, Any]:
        """
        Calculate Root Mean Square Fluctuation (RMSF).
        
        Returns:
            Dictionary containing RMSF data and metadata
        """
        print("   📊 Calculating RMSF...")
    
        per_chain = {}

        for at in self.trajectory.topology.atoms:
            if at.name == "CA" :
                continue


        
        try:

            chains = self._get_protein_chains()
            
            if not chains:
                print("   ⚠️  No protein chains found")
                return None
            

            indices = dict()
        
            for i in chains.keys():
                
                indices[f"chain_{i}"] = [atom.index for atom in self.trajectory.topology.atoms if atom.name == 'CA' and atom.residue.chain.index == i]

            
            # For RMSF, typically use CA atoms for proteins
            ca_indices = [atom.index for atom in self.trajectory.topology.atoms if atom.name == 'CA']
            
            
            if not ca_indices:
                print("   ⚠️  No CA atoms found, using all atoms")
                ca_indices = None
            
            # Calculate RMSF
            reference_frame = min(self.config.reference_frame, self.trajectory.n_frames - 1)
            per_chain = {}
            
          
            for chain,indice in indices.items():
                print(reference_frame)
                rmsf_values = md.rmsf(self.trajectory, self.trajectory, reference_frame, atom_indices=indice)
                per_chain[chain] = rmsf_values
       

            rmsf_values = []
           

            
            for chain, values in per_chain.items():

                
                rmsf_values = rmsf_values + list(values)
       
            
            
            #print(per_chain)



            
            # Convert to Angstroms
            rmsf_values *= 10  # nm to Angstroms
            
            print(f"   ✅ RMSF calculated: {len(rmsf_values)} values")
            print(f"      Mean RMSF: {np.mean(rmsf_values):.2f} Å")
            print(f"      Max RMSF: {np.max(rmsf_values):.2f} Å")
            
            # Get residue information for plotting
            
      
            return {
                'values': rmsf_values,
                "per_chain" : per_chain,
                
                'reference_frame': reference_frame,
                'n_residues': len(rmsf_values),
                'mean': np.mean(rmsf_values),
                'std': np.std(rmsf_values),
                'max': np.max(rmsf_values),
                'min': np.min(rmsf_values)
            }
            
        except Exception as e:
            print(f"   ❌ RMSF calculation failed: {str(e)}")
            return None
    
    def _calculate_distances(self) -> Dict[str, Any]:
        """
        Calculate distances between centers of mass of chains.
        
        Returns:
            Dictionary containing distance data and metadata
        """
        print("   📏 Calculating inter-chain distances...")
        
        try:
            # Get chains
            chains = list(self.trajectory.topology.chains)
            
            if len(chains) < 2:
                print("   ⚠️  Less than 2 chains found, skipping distance calculation")
                return None
            
            print(f"   Found {len(chains)} chains")
            
            # Calculate center of mass for each chain
            chain_distances = {}
            
            for i, chain1 in enumerate(chains):
                for j, chain2 in enumerate(chains[i+1:], i+1):
                    # Get atom indices for each chain
                    chain1_atoms = [atom.index for atom in chain1.atoms if atom.name == 'CA']
                    chain2_atoms = [atom.index for atom in chain2.atoms if atom.name == 'CA']
                    
                    # Calculate center of mass distances
                    distances = []
                    for frame in range(self.trajectory.n_frames):
                        com1 = md.compute_center_of_mass(self.trajectory[frame].atom_slice(chain1_atoms))
                        com2 = md.compute_center_of_mass(self.trajectory[frame].atom_slice(chain2_atoms))
                        dist = np.linalg.norm(com1 - com2) * 10  # Convert to Angstroms
                        distances.append(dist)
                    
                    chain_pair = f"Chain_{chain1.index}-Chain_{chain2.index}"
                    chain_distances[chain_pair] = np.array(distances)
                    
                    print(f"   ✅ {chain_pair}: Mean distance = {np.mean(distances):.2f} Å")
            
            return {
                'distances': chain_distances,
                'frames': np.arange(self.trajectory.n_frames),
                'n_chains': len(chains),
                'chain_pairs': list(chain_distances.keys())
            }
            
        except Exception as e:
            print(f"   ❌ Distance calculation failed: {str(e)}")
            return None
    
    def _calculate_radius_gyration(self) -> Dict[str, Any]:
        """
        Calculate radius of gyration.
        
        Returns:
            Dictionary containing radius of gyration data and metadata
        """
        print("   🎯 Calculating radius of gyration...")
        
        try:
            rg_values = md.compute_rg(self.trajectory)
            rg_values *= 10  # Convert to Angstroms
            
            print(f"   ✅ Radius of gyration calculated: {len(rg_values)} values")
            print(f"      Mean Rg: {np.mean(rg_values):.2f} Å")
            print(f"      Std Rg: {np.std(rg_values):.2f} Å")
            
            return {
                'values': rg_values,
                'frames': np.arange(len(rg_values)),
                'mean': np.mean(rg_values),
                'std': np.std(rg_values),
                'max': np.max(rg_values),
                'min': np.min(rg_values)
            }
            
        except Exception as e:
            print(f"   ❌ Radius of gyration calculation failed: {str(e)}")
            return None
    
    def _analyze_secondary_structure(self) -> Dict[str, Any]:
        """
        Analyze secondary structure evolution.
        
        Returns:
            Dictionary containing secondary structure data and metadata
        """
        print("   🧬 Analyzing secondary structure...")
        
        try:
            # Calculate secondary structure using DSSP
            ss = md.compute_dssp(self.trajectory, simplified=True)
            
            # Count secondary structure elements over time
            ss_counts = {}
            ss_types = ['H', 'E', 'C']  # Helix, Sheet, Coil
            
            for ss_type in ss_types:
                counts = np.sum(ss == ss_type, axis=1)
                ss_counts[ss_type] = counts
            
            print(f"   ✅ Secondary structure analyzed: {ss.shape[0]} frames, {ss.shape[1]} residues")
            
            return {
                'ss_matrix': ss,
                'ss_counts': ss_counts,
                'frames': np.arange(ss.shape[0]),
                'n_residues': ss.shape[1],
                'ss_types': ss_types
            }
            
        except Exception as e:
            print(f"   ❌ Secondary structure analysis failed: {str(e)}")
            print("   Note: DSSP analysis requires protein structures")
            return None
    
    def _get_atom_indices(self, selection: str) -> List[int]:
        """
        Get atom indices based on selection string, excluding water molecules.
        
        Args:
            selection: Atom selection type ('all', 'backbone', 'ca', 'heavy')
            
        Returns:
            List of atom indices (water molecules excluded)
        """
        # Get protein atoms only (exclude water and ions)
        protein_atoms = [atom for atom in self.trajectory.topology.atoms 
                        if atom.residue.name not in ['HOH', 'WAT', 'TIP', 'H2O', 'Na+', 'Cl-', 'K+', 'Mg2+', 'Ca2+', 'Zn2+', 'SO4', 'PO4']]
        
        if selection == 'all':
            return [atom.index for atom in protein_atoms]
        elif selection == 'backbone':
            return [atom.index for atom in protein_atoms 
                   if atom.name in ['N', 'CA', 'C', 'O']]
        elif selection == 'ca':
            return [atom.index for atom in protein_atoms 
                   if atom.name == 'CA']
        elif selection == 'heavy':
            return [atom.index for atom in protein_atoms 
                   if atom.element.symbol != 'H']
        else:
            raise ValueError(f"Unknown atom selection: {selection}")
    
    def _get_protein_chains(self) -> Dict[str, List[int]]:
        """
        Get atom indices for each protein chain separately.
        
        Returns:
            Dictionary mapping chain IDs to lists of atom indices
        """
        chains = {}
        
        for atom in self.trajectory.topology.atoms:
            # Skip water molecules and ions
            if not atom.residue.is_protein :
                continue
            
            chain_id = atom.residue.chain.index
            if chain_id not in chains:
                chains[chain_id] = []
            chains[chain_id].append(atom.index)

        
        return chains
    
    def _get_chain_atom_indices(self, chain_id: int, selection: str) -> List[int]:
        """
        Get atom indices for a specific chain based on selection string.
        
        Args:
            chain_id: Chain identifier
            selection: Atom selection type ('all', 'backbone', 'ca', 'heavy')
            
        Returns:
            List of atom indices for the specified chain
        """
        chain_atoms = [atom for atom in self.trajectory.topology.atoms 
                      if atom.residue.chain.index == chain_id and 
                      atom.residue.name not in ['HOH', 'WAT', 'TIP', 'H2O', 'Na+', 'Cl-', 'K+', 'Mg2+', 'Ca2+', 'Zn2+', 'SO4', 'PO4']]
        
        if selection == 'all':
            return [atom.index for atom in chain_atoms]
        elif selection == 'backbone':
            return [atom.index for atom in chain_atoms 
                   if atom.name in ['N', 'CA', 'C', 'O']]
        elif selection == 'ca':
            return [atom.index for atom in chain_atoms 
                   if atom.name == 'CA']
        elif selection == 'heavy':
            return [atom.index for atom in chain_atoms 
                   if atom.element.symbol != 'H']
        else:
            raise ValueError(f"Unknown atom selection: {selection}")
    
    def _generate_outputs(self, results: Dict[str, Any]):
        """
        Generate plots and save data files.
        
        Args:
            results: Dictionary containing analysis results
        """
        print("\n📊 Generating outputs...")
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        for analysis_type, data in results.items():
            if data is None:
                continue
                
            print(f"   📈 Generating {analysis_type} outputs...")
            
            
            if analysis_type == 'rmsd':
                self._plot_rmsd(data)
            
            elif analysis_type == 'rmsf':
                self._plot_rmsf(data)
            elif analysis_type == 'distances':
                self._plot_distances(data)
            elif analysis_type == 'radius_gyration':
                self._plot_radius_gyration(data)
            elif analysis_type == 'secondary_structure':
                self._plot_secondary_structure(data)
            
            
            # Save data if requested
            if self.config.save_data:
                self._save_data(analysis_type, data)
    
    def _plot_rmsd(self, data: Dict[str, Any]):
        """Generate RMSD plot."""
        if self.config.no_plots:
            return
   
        
        plt.figure(figsize=(10, 6))
        for i in data["per_chain"].keys():
            plt.plot(data["per_chain"][i]['frames'], data["per_chain"][i]['values'], linewidth=1.5,label = i)
        plt.title(f'RMSD vs Time ({data["atom_selection"]} atoms)', fontsize=14, fontweight='bold')
        plt.xlabel('Frame', fontsize=12)
        plt.ylabel('RMSD (Å)', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Add statistics text
        #stats_text = f"Mean: {data['mean']:.2f} Å\\nStd: {data['std']:.2f} Å\\nMax: {data['max']:.2f} Å"
        #plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
        #        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        output_file = self.output_dir / f'rmsd.{self.config.output_format}'
        plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
        plt.close()
        print(f"      ✅ RMSD plot saved: {output_file}")
    
    def _plot_rmsf(self, data: Dict[str, Any]):
        """Generate RMSF plot."""
        if self.config.no_plots:
            return
       
            
        plt.figure(figsize=(12, 6))

        for i in data["per_chain"].keys():
            plt.plot(data["per_chain"][i], linewidth=1.5,label = i)
        #plt.plot(data['residue_ids'], data['values'], linewidth=1.5)
        plt.title('RMSF per Residue', fontsize=14, fontweight='bold')
        plt.xlabel('Residue Number', fontsize=12)
        plt.ylabel('RMSF (Å)', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Add statistics text
        stats_text = f"Mean: {data['mean']:.2f} Å\\nStd: {data['std']:.2f} Å\\nMax: {data['max']:.2f} Å"
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        output_file = self.output_dir / f'rmsf.{self.config.output_format}'
        plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
        plt.close()
        print(f"      ✅ RMSF plot saved: {output_file}")
    
    def _plot_distances(self, data: Dict[str, Any]):
        """Generate distance plots."""
        if self.config.no_plots or not data['distances']:
            return
            
        n_pairs = len(data['distances'])
        fig, axes = plt.subplots(n_pairs, 1, figsize=(10, 4*n_pairs), squeeze=False)
        
        for i, (pair_name, distances) in enumerate(data['distances'].items()):
            ax = axes[i, 0]
            ax.plot(data['frames'], distances, linewidth=1.5)
            ax.set_title(f'Distance: {pair_name}', fontsize=12, fontweight='bold')
            ax.set_xlabel('Frame')
            ax.set_ylabel('Distance (Å)')
            ax.grid(True, alpha=0.3)
            
            # Add statistics
            mean_dist = np.mean(distances)
            std_dist = np.std(distances)
            stats_text = f"Mean: {mean_dist:.2f} Å\\nStd: {std_dist:.2f} Å"
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        output_file = self.output_dir / f'distances.{self.config.output_format}'
        plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
        plt.close()
        print(f"      ✅ Distance plots saved: {output_file}")
    
    def _plot_radius_gyration(self, data: Dict[str, Any]):
        """Generate radius of gyration plot."""
        if self.config.no_plots:
            return
            
        plt.figure(figsize=(10, 6))
        plt.plot(data['frames'], data['values'], linewidth=1.5)
        plt.title('Radius of Gyration vs Time', fontsize=14, fontweight='bold')
        plt.xlabel('Frame', fontsize=12)
        plt.ylabel('Radius of Gyration (Å)', fontsize=12)
        plt.grid(True, alpha=0.3)
        
        # Add statistics text
        stats_text = f"Mean: {data['mean']:.2f} Å\\nStd: {data['std']:.2f} Å"
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        output_file = self.output_dir / f'radius_gyration.{self.config.output_format}'
        plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
        plt.close()
        print(f"      ✅ Radius of gyration plot saved: {output_file}")
    
    def _plot_secondary_structure(self, data: Dict[str, Any]):
        """Generate secondary structure plots."""
        if self.config.no_plots:
            return
            
        # Plot secondary structure evolution
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Heatmap of secondary structure
        im = ax1.imshow(data['ss_matrix'].T, aspect='auto', cmap='viridis', origin='lower')
        ax1.set_title('Secondary Structure Evolution', fontsize=14, fontweight='bold')
        ax1.set_xlabel('Frame')
        ax1.set_ylabel('Residue')
        plt.colorbar(im, ax=ax1, label='SS Type')
        
        # Secondary structure counts over time
        for ss_type, counts in data['ss_counts'].items():
            ax2.plot(data['frames'], counts, label=f'{ss_type} ({"Helix" if ss_type=="H" else "Sheet" if ss_type=="E" else "Coil"})', linewidth=2)
        
        ax2.set_title('Secondary Structure Content vs Time', fontsize=14, fontweight='bold')
        ax2.set_xlabel('Frame')
        ax2.set_ylabel('Number of Residues')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_file = self.output_dir / f'secondary_structure.{self.config.output_format}'
        plt.savefig(output_file, dpi=self.config.dpi, bbox_inches='tight')
        plt.close()
        print(f"      ✅ Secondary structure plots saved: {output_file}")
    
    def _save_data(self, analysis_type: str, data: Dict[str, Any]):
        """Save analysis data to CSV files."""
        try:
            if analysis_type == 'rmsd':
                df = pd.DataFrame({
                    'Frame': data['frames'],
                    'RMSD_A': data['values']
                })
                output_file = self.output_dir / 'rmsd_data.csv'
                
            elif analysis_type == 'rmsf':
                df = pd.DataFrame({
                    'Residue_ID': data['residue_ids'],
                    'Residue_Name': data['residue_names'],
                    'RMSF_A': data['values']
                })
                output_file = self.output_dir / 'rmsf_data.csv'
                
            elif analysis_type == 'distances':
                df_data = {'Frame': data['frames']}
                for pair_name, distances in data['distances'].items():
                    df_data[f'{pair_name}_A'] = distances
                df = pd.DataFrame(df_data)
                output_file = self.output_dir / 'distances_data.csv'
                
            elif analysis_type == 'radius_gyration':
                df = pd.DataFrame({
                    'Frame': data['frames'],
                    'Rg_A': data['values']
                })
                output_file = self.output_dir / 'radius_gyration_data.csv'
                
            elif analysis_type == 'secondary_structure':
                df_data = {'Frame': data['frames']}
                for ss_type, counts in data['ss_counts'].items():
                    df_data[f'SS_{ss_type}'] = counts
                df = pd.DataFrame(df_data)
                output_file = self.output_dir / 'secondary_structure_data.csv'
            
            else:
                return
            
            df.to_csv(output_file, index=False)
            print(f"      ✅ Data saved: {output_file}")
            
        except Exception as e:
            print(f"      ❌ Failed to save {analysis_type} data: {str(e)}")
    
    def _display_completion_summary(self, results: Dict[str, Any]):

        """Display analysis completion summary."""
        print("\n" + "="*60)
        print("🎉 ANALYSIS COMPLETED SUCCESSFULLY")
        print("="*60)
        
        completed_analyses = [name for name, data in results.items() if data is not None]
        failed_analyses = [name for name, data in results.items() if data is None]
        
        print(f"\n✅ Completed Analyses ({len(completed_analyses)}):")
        for analysis in completed_analyses:
            print(f"   • {analysis.replace('_', ' ').title()}")
        
        if failed_analyses:
            print(f"\n❌ Failed Analyses ({len(failed_analyses)}):")
            for analysis in failed_analyses:
                print(f"   • {analysis.replace('_', ' ').title()}")
        
        print(f"\n📁 Output Directory: {self.output_dir}")
        
        # List generated files
        output_files = list(self.output_dir.glob('*'))
        if output_files:
            print(f"📄 Generated Files ({len(output_files)}):")
            for file in sorted(output_files):
                print(f"   • {file.name}")
        
        print("\n" + "="*60)
        print("Analysis complete! Check the output directory for results.")
        print("="*60)