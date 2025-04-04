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

class VAEMoleculeGenerator:
    def __init__(self, latent_dim=256, vocab_size=1000, use_chembl_scaffolds=True):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._vae_trained = False
        
        # Traditional fragments
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
        
        # Load ChEMBL scaffolds 
        if use_chembl_scaffolds:
            self.scaffolds = self._load_chembl_scaffolds()
            # Cluster scaffolds to reduce redundancy
            self.scaffolds = self._cluster_scaffolds(self.scaffolds)
        else:
            self.scaffolds = self._default_scaffolds
        
        # VAE Components
        self.encoder = self._build_encoder(latent_dim)
        self.decoder = self._build_decoder(latent_dim)
        self.optimizer = optim.Adam(
            list(self.encoder.parameters()) +
            list(self.decoder.parameters()),
            lr=0.001
        )


    def _load_chembl_scaffolds(self, max_scaffolds=1000):
        try:
            print("Loading ChEMBL scaffolds...")
            
            # Load from a local file if available
            try:
                scaffolds_df = pd.read_csv('chembl_scaffolds.csv')
                scaffolds = scaffolds_df['scaffold_smiles'].tolist()[:max_scaffolds]
                print(f"Loaded {len(scaffolds)} scaffolds from local file")
                return [(s, f"ChEMBL_Scaffold_{i}") for i, s in enumerate(scaffolds)]
            except:
                pass
            
            # Extracting from ChEMBL API
            try:
                print("Fetching compounds from ChEMBL API...")
                from chembl_webresource_client.new_client import new_client
                
                compounds_api = new_client.molecule
                compounds = compounds_api.filter(
                    molecule_structures__canonical_smiles__isnull='false'
                ).only('molecule_structures')[:1000]
                
                # Extracting SMILES
                smiles_list = [
                    compound['molecule_structures']['canonical_smiles'] 
                    for compound in compounds 
                    if compound['molecule_structures']
                ]
                
                # Murcko scaffolds
                from rdkit.Chem.Scaffolds import MurckoScaffold
                
                scaffolds = set()
                for smiles in smiles_list:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                        scaffold_smiles = Chem.MolToSmiles(scaffold)
                        if scaffold_smiles and len(scaffold_smiles) > 5:  #tiny scaffolds
                            scaffolds.add(scaffold_smiles)
                
                scaffolds = list(scaffolds)[:max_scaffolds]
                print(f"Extracted {len(scaffolds)} scaffolds from ChEMBL API")
                
                # Save for future use
                import pandas as pd
                pd.DataFrame(scaffolds, columns=['scaffold_smiles']).to_csv('chembl_scaffolds.csv', index=False)
                
                return [(s, f"ChEMBL_Scaffold_{i}") for i, s in enumerate(scaffolds)]
            except Exception as e:
                print(f"Error fetching from ChEMBL API: {str(e)}")
            
            # default scaffolds
            print("Using default scaffolds as fallback")
            return self._default_scaffolds
            
        except Exception as e:
            print(f"Error loading ChEMBL scaffolds: {str(e)}")
            return self._default_scaffolds




    def _cluster_scaffolds(self, scaffolds, threshold=0.7):
        print("Clustering scaffolds...")

        # scaffolds to molecules and fingerprints
        mols = []
        valid_scaffolds = []

        for scaffold, name in scaffolds:
            try:
                mol = Chem.MolFromSmiles(scaffold)
                if mol:
                    mols.append(mol)
                    valid_scaffolds.append((scaffold, name))
            except:
                continue
            
        if len(mols) <= 1:
            return valid_scaffolds

        fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for mol in mols]

        # distance matrix (1 - similarity)
        n = len(fingerprints)
        distance_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i+1, n):
                sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                distance_matrix[i, j] = 1.0 - sim
                distance_matrix[j, i] = 1.0 - sim

        # hierarchical clustering
        Z = linkage(distance_matrix, method='ward')

        # Form flat clusters
        clusters = fcluster(Z, t=1 - threshold, criterion='distance')

        # representatives from each cluster
        clustered_scaffolds = []
        for cluster_id in set(clusters):
            cluster_indices = np.where(clusters == cluster_id)[0]
            # Use the first scaffold in each cluster as representative
            clustered_scaffolds.append(valid_scaffolds[cluster_indices[0]])

        print(f"Reduced {len(valid_scaffolds)} scaffolds to {len(clustered_scaffolds)} clusters")
        return clustered_scaffolds




    def _build_encoder(self, latent_dim):
        return nn.Sequential(
            nn.Linear(2048, 512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Linear(256, latent_dim * 2)
        ).to(self.device)
        
    def _build_decoder(self, latent_dim):
        return nn.Sequential(
            nn.Linear(latent_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 512),
            nn.ReLU(),
            nn.Linear(512, 2048),
            nn.Sigmoid()
        ).to(self.device)
        
    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std
        
    def _train_vae(self, smiles_list, epochs=50):
        fingerprints = []
        valid_smiles = []
        
        for smile in smiles_list:
            mol = Chem.MolFromSmiles(smile)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                arr = np.zeros(2048, dtype=np.int32)
                AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
                fingerprints.append(arr)
                valid_smiles.append(smile)
                
        if not fingerprints:
            return
            
        # X = torch.FloatTensor(fingerprints).to(self.device)
        X = torch.FloatTensor(np.array(fingerprints)).to(self.device)

        best_loss = float('inf')
        
        for epoch in range(epochs):
            self.optimizer.zero_grad()
            
            enc_output = self.encoder(X)
            mu, log_var = torch.chunk(enc_output, 2, dim=1)
            z = self.reparameterize(mu, log_var)
            recon_x = self.decoder(z)
            
            recon_loss = nn.functional.binary_cross_entropy(recon_x, X)
            kl_loss = -0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())
            loss = recon_loss + kl_loss
            
            loss.backward()
            self.optimizer.step()
            
            # reconstruction accuracy
            with torch.no_grad():
                preds = (recon_x > 0.5).float()
                acc = (preds == X).float().mean()
                
            if loss.item() < best_loss:
                best_loss = loss.item()
                
            print(f"Epoch {epoch+1}/{epochs} | Loss: {loss.item():.4f} | Acc: {acc.item():.2%}")
            
        self._vae_trained = True
        self.training_fps = fingerprints
        self.training_smiles = valid_smiles
        
    def generate_molecules(self, n_molecules=100, training_data=None, target_properties=None, 
                      receptor_file=None, docking_threshold=-7.0):
        
        if training_data:
            print("🚀 Training VAE model...")
            self._train_vae(training_data)

        molecules = []
        seen = set()
        docking_enabled = receptor_file is not None

        if docking_enabled:
            print("Docking filter enabled - molecules will be filtered by docking score")

        with tqdm(total=n_molecules, desc="Generating molecules") as pbar:
            while len(molecules) < n_molecules:
                # Hybrid generation strategy
                if self._vae_trained and random.random() > 0.3:  # 70% VAE, 30% traditional
                    mol = self._vae_generate_molecule()
                    if not mol:
                        mol = self._generate_hybrid_molecule()
                else:
                    mol = self._generate_hybrid_molecule()

                if mol:
                    smile = Chem.MolToSmiles(mol)

                    passes_filters = smile not in seen and self._pass_filters(mol)

                    # docking filter if enabled
                    if passes_filters and docking_enabled:
                        docking_score = self._dock_molecule(mol, receptor_file)
                        passes_filters = docking_score <= docking_threshold

                    if passes_filters:
                        seen.add(smile)
                        molecules.append(mol)
                        pbar.update(1)

        return [Chem.MolToSmiles(mol) for mol in molecules]
  
    def _vae_generate_molecule(self):
        try:
            with torch.no_grad():
                # Using optimization with 20% probability, random sampling otherwise
                if random.random() < 0.2:
                    z = self._optimize_latent_space(n_calls=10)
                else:
                    z = torch.randn(1, self.decoder[0].in_features).to(self.device)

                if z is None:
                    return None

                recon_x = self.decoder(z)

                # Probabilistic sampling
                gen_fp_arr = (torch.rand_like(recon_x) < recon_x).cpu().numpy().astype(int)[0]

                # Create fingerprint object
                gen_fp = DataStructs.ExplicitBitVect(2048)
                for i in range(2048):
                    if gen_fp_arr[i]:
                        gen_fp.SetBit(i)

                # nearest neighbor
                max_sim = -1
                best_mol = None

                for train_fp, smile in zip(self.training_fps, self.training_smiles):
                    train_fpv = DataStructs.ExplicitBitVect(2048)
                    for i in range(2048):
                        if train_fp[i]:
                            train_fpv.SetBit(i)

                    sim = DataStructs.TanimotoSimilarity(gen_fp, train_fpv)
                    if sim > max_sim:
                        max_sim = sim
                        best_mol = Chem.MolFromSmiles(smile)

                # Only return novel molecules
                if max_sim < 0.6:
                    return self._optimize_molecule(gen_fp_arr)

                return best_mol if max_sim < 0.9 else None

        except Exception as e:
            print(f"VAE generation failed: {e}")
            return None
      

    def _dock_molecule(self, mol, receptor_file=None):
        # This is a placeholder for integrating with docking software
        # In a real implementation, this would call AutoDock Vina or similar
        try:
            # Mock implementation - returns random score between -12 (good) and 0 (poor)
            if mol:
                # In a real implementation:
                # 1. Convert mol to 3D
                # 2. Prepare for docking (add hydrogens, charges)
                # 3. Call docking software
                # 4. Parse and return results
                
                # For now, just return a property-based estimate
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                rotatable_bonds = Descriptors.NumRotatableBonds(mol)
                
                # Simple heuristic: smaller, more rigid molecules tend to dock better
                mock_score = -10.0 + (mw / 500) + (logp / 5) + (rotatable_bonds / 10)
                return max(-12.0, min(0.0, mock_score))
            return 0.0
        except:
            return 0.0

        
    def _generate_hybrid_molecule(self):
        scaffold, _ = random.choice(self.scaffolds)
        mol = self._modify_scaffold(scaffold)
        return mol
        
    def _modify_scaffold(self, scaffold):
        mol = Chem.MolFromSmiles(scaffold)
        if not mol:
            return None
            
        for _ in range(random.randint(1, 3)):
            try:
                mol = self._add_fragment(mol)
                Chem.SanitizeMol(mol)
            except:
                return None
                
        stereoisomers = list(EnumerateStereoisomers(mol))
        return random.choice(stereoisomers) if stereoisomers else mol
        
    def _predict_fragment(self, mol):
        # Count functional groups
        has_aromatic = False
        has_hydroxyl = False
        has_amine = False
        has_carbonyl = False

        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                has_aromatic = True
            if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() > 0:
                has_hydroxyl = True
            if atom.GetSymbol() == 'N':
                has_amine = True

        # carbonyl groups
        patt = Chem.MolFromSmarts('C=O')
        if mol.HasSubstructMatch(patt):
            has_carbonyl = True

        # Molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)

        # Rules for fragment selection
        if has_aromatic and not has_hydroxyl:
            return ("O", "Hydroxyl")  # Add hydroxyl to aromatic rings
        elif has_aromatic and not has_amine:
            return ("N", "Amine")  # Add amine to aromatic rings
        elif has_hydroxyl and not has_carbonyl:
            return ("CO", "Carbonyl")  # Add carbonyl near hydroxyl
        elif mw < 200:
            return ("c1ccccc1", "Benzene")  # Add aromatic ring to small molecules
        elif logp > 3:
            return ("O", "Hydroxyl")  # Add hydroxyl to lipophilic molecules
        else:
            # Default: random selection
            return random.choice(self.fragments)

    def _add_fragment(self, mol):
        rwmol = Chem.RWMol(mol)
        atoms = [atom for atom in rwmol.GetAtoms() if atom.GetTotalNumHs() > 0]

        if not atoms:
            return mol

        atom = random.choice(atoms)
        fragment, _ = self._predict_fragment(mol)
        frag_mol = Chem.MolFromSmiles(fragment)

        if not frag_mol:
            fragment, _ = random.choice(self.fragments)
            frag_mol = Chem.MolFromSmiles(fragment)

        combined = Chem.CombineMols(rwmol, frag_mol)
        ed_combined = Chem.EditableMol(combined)

        ed_combined.AddBond(
            atom.GetIdx(),
            rwmol.GetNumAtoms(),
            order=Chem.BondType.SINGLE
        )

        return ed_combined.GetMol()
    

    def _optimize_latent_space(self, n_calls=10, property_func=None):
        try:
            if not self._vae_trained:
                print("VAE not trained, cannot optimize latent space")
                return None

            latent_dim = self.decoder[0].in_features

            # maximize QED (drug-likeness)
            if property_func is None:
                from rdkit.Chem import QED

                def property_func(mol):
                    if mol:
                        return -QED.qed(mol)  
                    return 0

            # optimizing
            def objective(x):

                z = torch.FloatTensor(x).reshape(1, -1).to(self.device)

                with torch.no_grad():
                    recon_x = self.decoder(z)
                    gen_fp_arr = (torch.rand_like(recon_x) < recon_x).cpu().numpy().astype(int)[0]

                    #nearest molecule in training
                    gen_fp = DataStructs.ExplicitBitVect(2048)
                    for i in range(2048):
                        if gen_fp_arr[i]:
                            gen_fp.SetBit(i)

                    max_sim = -1
                    best_mol = None

                    for train_fp, smile in zip(self.training_fps, self.training_smiles):
                        train_fpv = DataStructs.ExplicitBitVect(2048)
                        for i in range(2048):
                            if train_fp[i]:
                                train_fpv.SetBit(i)

                        sim = DataStructs.TanimotoSimilarity(gen_fp, train_fpv)
                        if sim > max_sim:
                            max_sim = sim
                            best_mol = Chem.MolFromSmiles(smile)

                    return property_func(best_mol)

            #  search space
            space = [Real(-3.0, 3.0, name=f'z_{i}') for i in range(latent_dim)]

            #optimization
            result = gp_minimize(objective, space, n_calls=n_calls, random_state=42)

            # best latent vector
            best_z = result.x

            return torch.FloatTensor(best_z).reshape(1, -1).to(self.device)

        except ImportError:
            print("scikit-optimize not installed, using random sampling instead")
            return torch.randn(1, self.decoder[0].in_features).to(self.device)
        except Exception as e:
            print(f"Error in latent space optimization: {e}")
            return torch.randn(1, self.decoder[0].in_features).to(self.device)


        
    def _pass_filters(self, mol):
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        return (150 <= mw <= 600) and (-2 <= logp <= 5) and \
               (Descriptors.NumHDonors(mol) <= 5) and \
               (Descriptors.NumHAcceptors(mol) <= 10)

class MoleculeGenerator(VAEMoleculeGenerator):
    def __init__(self):
        super().__init__()
