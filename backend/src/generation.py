import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
import random
from tqdm import tqdm

class VAEMoleculeGenerator:
    def __init__(self, latent_dim=256, vocab_size=1000):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._vae_trained = False
        
        # Traditional fragments
        self.scaffolds = [
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
        
        # VAE 
        self.encoder = self._build_encoder(latent_dim)
        self.decoder = self._build_decoder(latent_dim)
        self.optimizer = optim.Adam(
            list(self.encoder.parameters()) +
            list(self.decoder.parameters()),
            lr=0.001
        )
        
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
        
    def generate_molecules(self, n_molecules=100, training_data=None, target_properties=None):
        if training_data:
            print("🚀 Training VAE model...")
            self._train_vae(training_data)
            
        molecules = []
        seen = set()
        
        with tqdm(total=n_molecules, desc="Generating molecules") as pbar:
            while len(molecules) < n_molecules:
                # Hybrid generation strategy
                if self._vae_trained and random.random() > 0.3: # 70% VAE, 30% traditional
                    mol = self._vae_generate_molecule()
                    if not mol:
                        mol = self._generate_hybrid_molecule()
                else:
                    mol = self._generate_hybrid_molecule()
                    
                if mol:
                    smile = Chem.MolToSmiles(mol)
                    if smile not in seen and self._pass_filters(mol):
                        seen.add(smile)
                        molecules.append(mol)
                        pbar.update(1)
                        
        return [Chem.MolToSmiles(mol) for mol in molecules]
        
    def _vae_generate_molecule(self):
        """Generate novel molecules through VAE latent space"""
        try:
            with torch.no_grad():
                # Generating from latent space
                z = torch.randn(1, self.decoder[0].in_features).to(self.device)
                recon_x = self.decoder(z)
                
                # Probabilistic sampling
                gen_fp_arr = (torch.rand_like(recon_x) < recon_x).cpu().numpy().astype(int)[0]
                
                # fingerprint object
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
                    
                if max_sim < 0.6:
                    return self._optimize_molecule(gen_fp_arr)
                    
                return best_mol if max_sim < 0.9 else None
                
        except Exception as e:
            print(f"VAE generation failed: {e}")
            return None
            
    # def _optimize_molecule(self, fp_arr):
    #     # Placeholder for genetic optimization
    #     return None 
        
    def _generate_hybrid_molecule(self):
        """Traditional fragment-based generation"""
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
        
    def _add_fragment(self, mol):
        rwmol = Chem.RWMol(mol)
        atoms = [atom for atom in rwmol.GetAtoms() if atom.GetTotalNumHs() > 0]
        
        if not atoms:
            return mol
            
        atom = random.choice(atoms)
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
        
    def _pass_filters(self, mol):
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        return (150 <= mw <= 600) and (-2 <= logp <= 5) and \
               (Descriptors.NumHDonors(mol) <= 5) and \
               (Descriptors.NumHAcceptors(mol) <= 10)

class MoleculeGenerator(VAEMoleculeGenerator):
    def __init__(self):
        super().__init__()
