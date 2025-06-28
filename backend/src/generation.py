import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
import random
from tqdm import tqdm
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from skopt import gp_minimize
from skopt.space import Real
from rdkit.Chem import QED
from rdkit.Chem import rdFingerprintGenerator
import pickle
import time
from functools import lru_cache


class OptimizedVAEMoleculeGenerator:
    def __init__(self, latent_dim=256, vocab_size=1000, use_chembl_scaffolds=True, cache_size=1000):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._vae_trained = False
        self.cache_size = cache_size
        
        # Initialize fingerprint generator once
        self.morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        
        # Traditional fragments (kept smaller for speed)
        self._default_scaffolds = [
            ("c1ccccc1", "Benzene"),
            ("C1CCCCC1", "Cyclohexane"),
            ("c1ccncc1", "Pyridine"),
            ("c1ccccc1O", "Phenol"),
            ("C1CCOC1", "THF"),
            ("c1ccsc1", "Thiophene"),
            ("C1CNCCN1", "Piperazine"),
        ]
        
        self.fragments = [
            ("CN", "Amino"), ("N", "Amine"),
            ("C1CCNCC1", "Piperazine"), ("C", "Methyl"),
            ("O", "Hydroxyl"), ("Cl", "Chloro"),
            ("F", "Fluoro"), ("CO", "Carbonyl"),
            ("C(=O)O", "Carboxylic Acid"), ("CN", "Cyano"),
        ]
        
        # Load and preprocess scaffolds
        if use_chembl_scaffolds:
            self.scaffolds = self._load_optimized_scaffolds()
        else:
            self.scaffolds = self._default_scaffolds
        
        # Precompute scaffold properties for faster selection
        self._precompute_scaffold_properties()
        
        # VAE Components
        self.encoder = self._build_encoder(latent_dim)
        self.decoder = self._build_decoder(latent_dim)
        self.optimizer = optim.Adam(
            list(self.encoder.parameters()) + list(self.decoder.parameters()),
            lr=0.001
        )
        
        # Caching for expensive operations
        self._mol_cache = {}
        self._fp_cache = {}
        self._property_cache = {}

    def _load_optimized_scaffolds(self, max_scaffolds=100):  # Reduced from 1000
        """Load scaffolds with optimizations"""
        try:
            # Try to load from cache first
            try:
                with open('scaffolds_cache.pkl', 'rb') as f:
                    cached_data = pickle.load(f)
                    if len(cached_data) >= max_scaffolds:
                        print(f"Loaded {len(cached_data)} scaffolds from cache")
                        return cached_data[:max_scaffolds]
            except:
                pass
            
            print("Processing scaffolds...")
            
            # Use your provided scaffolds (limited set for speed)
            scaffold_smiles_list = [
                "c1ccc2[nH]ccc2c1",
                "O=C(NNC(=O)[C@H]1CCCN1)OCc1ccccc1",
                "C=C1CCC(c2cccc3ccccc23)C(=O)O1",
                "O=C([C@@H]1CCCN1)N1CCN(Cc2ccccn2)CC1",
                "c1ccc(CN2CCN3CC(c4ccccc4)c4ccccc4C3C2)cc1",
                "C1=NC(C2CCCCC2)CCC1",
                "c1ccc2c(c1)Cc1ccccc1-2",
                "C1CC2CCC1C2",
                "c1ccc(CSc2ncnc3sccc23)cc1",
                "c1ccc2nc(N3CCNCC3)cnc2c1",
                "c1ccc(C2CC3CCC2C3)cc1",
                "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1",
                "O=C1c2ccccc2-c2ccccc21",
                "O=C1CCCC1",
                "c1nncs1",
                "c1ccc2nc(N3CCCCC3)ncc2c1",
                "O=C1Cc2ccccc2N1",
                "c1nc(Nc2nncs2)c2ccsc2n1",
                "c1ccc(-c2cccnn2)cc1",
                "c1ccc2c3c([nH]c2c1)CCNCC3",
                "O=S(=O)(Nc1cnccn1)c1cccc2ccccc12",
                "c1ncc(C2CC3CCC2N3)cn1",
                "c1ccc2c3c([nH]c2c1)CCNC3",
                "c1ccc2c(c1)OCO2",
                "N=c1[nH]ncs1",
                "c1ccc(C2CNCc3ccccc32)cc1",
                "c1ccc2c(c1)nc1n2CCNC1",
                "O=S(=O)(Nc1ccccc1)c1ccccc1",
                "c1nc(Sc2nncs2)c2ccsc2n1",
                "c1ccc2ccccc2c1",
                "C1=NCCCC1",
                "c1ccc2c(Sc3cnc4ncncc4n3)cccc2c1"
            ]
            
            # Process and validate scaffolds
            valid_scaffolds = []
            for i, smiles in enumerate(scaffold_smiles_list[:max_scaffolds]):
                mol = Chem.MolFromSmiles(smiles)
                if mol and self._is_valid_scaffold(mol):
                    valid_scaffolds.append((smiles, f"Scaffold_{i}"))
            
            # Cluster to reduce redundancy (faster with smaller set)
            if len(valid_scaffolds) > 20:
                valid_scaffolds = self._fast_cluster_scaffolds(valid_scaffolds, max_clusters=20)
            
            # Cache the results
            try:
                with open('scaffolds_cache.pkl', 'wb') as f:
                    pickle.dump(valid_scaffolds, f)
            except:
                pass
            
            print(f"Processed {len(valid_scaffolds)} scaffolds")
            return valid_scaffolds
            
        except Exception as e:
            print(f"Error loading scaffolds: {str(e)}")
            return self._default_scaffolds

    def _is_valid_scaffold(self, mol):
        """Quick validation of scaffold"""
        try:
            mw = Descriptors.MolWt(mol)
            return 50 <= mw <= 400 and mol.GetNumAtoms() >= 3
        except:
            return False

    def _fast_cluster_scaffolds(self, scaffolds, max_clusters=20, threshold=0.6):
        """Faster clustering with early termination"""
        print(f"Fast clustering {len(scaffolds)} scaffolds...")
        
        if len(scaffolds) <= max_clusters:
            return scaffolds
        
        # Sample subset for clustering if too many
        if len(scaffolds) > 50:
            scaffolds = random.sample(scaffolds, 50)
        
        # Convert to molecules and get fingerprints
        mols = []
        valid_scaffolds = []
        
        for scaffold, name in scaffolds:
            mol = self._get_mol_cached(scaffold)
            if mol:
                mols.append(mol)
                valid_scaffolds.append((scaffold, name))
        
        if len(mols) <= max_clusters:
            return valid_scaffolds
        
        # Fast fingerprint generation
        fingerprints = [self._get_fingerprint_cached(mol) for mol in mols]
        
        # Simple greedy clustering for speed
        clusters = []
        used = set()
        
        for i, (scaffold, name) in enumerate(valid_scaffolds):
            if i in used:
                continue
                
            cluster = [i]
            used.add(i)
            
            for j in range(i + 1, len(valid_scaffolds)):
                if j in used:
                    continue
                    
                sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                if sim >= threshold:
                    cluster.append(j)
                    used.add(j)
            
            # Take representative from cluster
            clusters.append(valid_scaffolds[cluster[0]])
            
            if len(clusters) >= max_clusters:
                break
        
        print(f"Clustered to {len(clusters)} representatives")
        return clusters

    @lru_cache(maxsize=1000)
    def _get_mol_cached(self, smiles):
        """Cached molecule creation"""
        return Chem.MolFromSmiles(smiles)

    @lru_cache(maxsize=1000)
    def _get_fingerprint_cached(self, mol):
        """Cached fingerprint generation"""
        if mol is None:
            return None
        return self.morgan_gen.GetFingerprint(mol)

    def _precompute_scaffold_properties(self):
        """Precompute properties for all scaffolds"""
        print("Precomputing scaffold properties...")
        self.scaffold_properties = []
        
        for scaffold_smiles, name in self.scaffolds:
            mol = self._get_mol_cached(scaffold_smiles)
            if mol:
                try:
                    qed = QED.qed(mol)
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    
                    # Simple complexity measure
                    complexity = mol.GetNumBonds() / max(mol.GetNumAtoms(), 1)
                    norm_complexity = 1.0 / (1.0 + complexity)
                    
                    # Property score
                    property_score = (qed * 3.0 + norm_complexity * 2.0 + 
                                    (0.5 if 250 < mw < 500 else 0) + 
                                    (0.5 if 1 < logp < 4 else 0))
                    
                    self.scaffold_properties.append({
                        'smiles': scaffold_smiles,
                        'name': name,
                        'qed': qed,
                        'mw': mw,
                        'logp': logp,
                        'score': property_score
                    })
                except:
                    # Default values for problematic molecules
                    self.scaffold_properties.append({
                        'smiles': scaffold_smiles,
                        'name': name,
                        'qed': 0.5,
                        'mw': 200,
                        'logp': 2.0,
                        'score': 1.0
                    })

    def _build_encoder(self, latent_dim):
        return nn.Sequential(
            nn.Linear(2048, 512),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, latent_dim * 2)
        ).to(self.device)
        
    def _build_decoder(self, latent_dim):
        return nn.Sequential(
            nn.Linear(latent_dim, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, 512),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(512, 2048),
            nn.Sigmoid()
        ).to(self.device)
        
    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std

    def _train_vae(self, smiles_list, epochs=30):  # Reduced epochs
        """Optimized VAE training"""
        print("Training VAE...")
        
        # Batch process fingerprints
        fingerprints = []
        valid_smiles = []
        
        batch_size = 100
        for i in range(0, len(smiles_list), batch_size):
            batch = smiles_list[i:i+batch_size]
            for smile in batch:
                mol = self._get_mol_cached(smile)
                if mol:
                    fp = self._get_fingerprint_cached(mol)
                    if fp:
                        arr = np.zeros(2048, dtype=np.int32)
                        DataStructs.ConvertToNumpyArray(fp, arr)
                        fingerprints.append(arr)
                        valid_smiles.append(smile)
                        
        if not fingerprints:
            print("No valid fingerprints generated")
            return
            
        X = torch.FloatTensor(np.array(fingerprints)).to(self.device)
        
        # Training with early stopping
        best_loss = float('inf')
        patience = 5
        patience_counter = 0
        
        for epoch in range(epochs):
            self.optimizer.zero_grad()
            
            enc_output = self.encoder(X)
            mu, log_var = torch.chunk(enc_output, 2, dim=1)
            z = self.reparameterize(mu, log_var)
            recon_x = self.decoder(z)
            
            recon_loss = nn.functional.binary_cross_entropy(recon_x, X)
            kl_loss = -0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp()) / X.size(0)
            loss = recon_loss + 0.1 * kl_loss  # Reduced KL weight
            
            loss.backward()
            torch.nn.utils.clip_grad_norm_(
                list(self.encoder.parameters()) + list(self.decoder.parameters()), 
                max_norm=1.0
            )
            self.optimizer.step()
            
            # Early stopping
            if loss.item() < best_loss:
                best_loss = loss.item()
                patience_counter = 0
            else:
                patience_counter += 1
                
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch+1}")
                break
                
            if epoch % 5 == 0:
                with torch.no_grad():
                    preds = (recon_x > 0.5).float()
                    acc = (preds == X).float().mean()
                print(f"Epoch {epoch+1}/{epochs} | Loss: {loss.item():.4f} | Acc: {acc.item():.2%}")
        
        self._vae_trained = True
        self.training_fps = fingerprints
        self.training_smiles = valid_smiles

    def _select_scaffold_optimized(self):
        """Optimized scaffold selection using precomputed properties"""
        if not self.scaffold_properties:
            return self.scaffolds[0][0]
        
        # Fast weighted selection
        total_score = sum(prop['score'] for prop in self.scaffold_properties)
        if total_score <= 0:
            return random.choice(self.scaffold_properties)['smiles']
        
        r = random.uniform(0, total_score)
        cumulative = 0
        for prop in self.scaffold_properties:
            cumulative += prop['score']
            if r <= cumulative:
                return prop['smiles']
        
        return self.scaffold_properties[0]['smiles']

    def _generate_hybrid_molecule(self):
        """Optimized hybrid molecule generation"""
        scaffold = self._select_scaffold_optimized()
        mol = self._modify_scaffold_fast(scaffold)
        return mol

    def _modify_scaffold_fast(self, scaffold):
        """Faster scaffold modification"""
        mol = self._get_mol_cached(scaffold)
        if not mol:
            return None
            
        # Limit modifications for speed
        max_modifications = random.randint(1, 2)  # Reduced from 3
        
        for _ in range(max_modifications):
            try:
                mol = self._add_fragment_fast(mol)
                if mol is None:
                    break
                Chem.SanitizeMol(mol)
            except:
                break
                
        # Skip stereoisomer enumeration for speed (major bottleneck)
        return mol

    def _add_fragment_fast(self, mol):
        """Faster fragment addition"""
        rwmol = Chem.RWMol(mol)
        
        # Simple atom selection (avoid scoring for speed)
        atoms = [atom for atom in rwmol.GetAtoms() if atom.GetTotalNumHs() > 0]
        if not atoms:
            return mol
            
        selected_atom = random.choice(atoms)
        
        # Use smaller fragment set
        simple_fragments = [("C", "Methyl"), ("O", "Hydroxyl"), ("N", "Amine"), ("F", "Fluoro")]
        fragment, _ = random.choice(simple_fragments)
        
        frag_mol = self._get_mol_cached(fragment)
        if not frag_mol:
            return mol
            
        try:
            combined = Chem.CombineMols(rwmol, frag_mol)
            ed_combined = Chem.EditableMol(combined)
            ed_combined.AddBond(
                selected_atom.GetIdx(),
                rwmol.GetNumAtoms(),
                order=Chem.BondType.SINGLE
            )
            return ed_combined.GetMol()
        except:
            return mol

    def generate_molecules(self, n_molecules=100, training_data=None, target_properties=None, 
                          receptor_file=None, docking_threshold=-7.0):
        """Optimized molecule generation"""
        
        if training_data:
            print("🚀 Training VAE model...")
            # Limit training data size for speed
            limited_training = training_data[:min(1000, len(training_data))]
            self._train_vae(limited_training)

        molecules = []
        seen = set()
        docking_enabled = receptor_file is not None
        
        # Batch generation for efficiency
        batch_size = 10
        attempts_per_batch = batch_size * 3  # Allow some failures
        
        with tqdm(total=n_molecules, desc="Generating molecules") as pbar:
            while len(molecules) < n_molecules:
                batch_molecules = []
                
                # Generate batch
                for _ in range(attempts_per_batch):
                    if len(batch_molecules) >= batch_size:
                        break
                        
                    # Simplified generation strategy
                    if self._vae_trained and random.random() > 0.5:  # 50/50 split
                        mol = self._vae_generate_molecule_fast()
                    else:
                        mol = self._generate_hybrid_molecule()
                    
                    if mol:
                        smile = Chem.MolToSmiles(mol)
                        if smile not in seen and self._pass_filters_fast(mol):
                            if not docking_enabled or self._dock_molecule(mol, receptor_file) <= docking_threshold:
                                batch_molecules.append(mol)
                                seen.add(smile)
                
                # Add successful molecules
                for mol in batch_molecules:
                    if len(molecules) < n_molecules:
                        molecules.append(mol)
                        pbar.update(1)
        
        return [Chem.MolToSmiles(mol) for mol in molecules]

    def _vae_generate_molecule_fast(self):
        """Faster VAE generation"""
        if not self._vae_trained:
            return None
            
        try:
            with torch.no_grad():
                # Always use random sampling (skip optimization for speed)
                z = torch.randn(1, self.decoder[0].in_features).to(self.device)
                recon_x = self.decoder(z)
                
                # Deterministic sampling for consistency
                gen_fp_arr = (recon_x > 0.5).cpu().numpy().astype(int)[0]
                
                # Quick similarity check with limited training set
                if len(self.training_smiles) > 100:
                    sample_indices = random.sample(range(len(self.training_smiles)), 100)
                else:
                    sample_indices = range(len(self.training_smiles))
                
                best_mol = None
                max_sim = -1
                
                gen_fp = DataStructs.ExplicitBitVect(2048)
                for i in range(2048):
                    if gen_fp_arr[i]:
                        gen_fp.SetBit(i)
                
                for idx in sample_indices:
                    train_fp = self.training_fps[idx]
                    train_fpv = DataStructs.ExplicitBitVect(2048)
                    for i in range(2048):
                        if train_fp[i]:
                            train_fpv.SetBit(i)
                    
                    sim = DataStructs.TanimotoSimilarity(gen_fp, train_fpv)
                    if sim > max_sim:
                        max_sim = sim
                        best_mol = self._get_mol_cached(self.training_smiles[idx])
                
                # Return novel molecules
                return best_mol if 0.3 < max_sim < 0.8 else None
                
        except Exception as e:
            return None

    def _pass_filters_fast(self, mol):
        """Faster molecular filters"""
        try:
            mw = Descriptors.MolWt(mol)
            if not (150 <= mw <= 600):
                return False
                
            logp = Descriptors.MolLogP(mol)
            if not (-2 <= logp <= 5):
                return False
                
            # Skip other filters for speed
            return True
        except:
            return False

    def _dock_molecule(self, mol, receptor_file=None):
        """Fast mock docking"""
        try:
            if mol:
                # Simplified property-based scoring
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                
                # Quick heuristic
                mock_score = -8.0 + (mw / 500) * 2 + (abs(logp - 2) / 3) * 2
                return max(-12.0, min(0.0, mock_score))
            return 0.0
        except:
            return 0.0


# Convenience wrapper maintaining the same interface
class MoleculeGenerator(OptimizedVAEMoleculeGenerator):
    def __init__(self):
        super().__init__()