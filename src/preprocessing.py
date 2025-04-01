from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.SaltRemover import SaltRemover
import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)
remover = SaltRemover()

class Preprocessor:
    def __init__(self):
        self.descriptors = [
            'MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors',
            'TPSA', 'NumRotatableBonds', 'RingCount'
        ]

    def process_smiles(self, smiles_series):
        """Convert SMILES to cleaned molecules and fingerprints"""
        valid_mols = []
        valid_smiles = []

        for smile in smiles_series:
            try:
                if pd.isna(smile):
                    continue

                # Sanitize and remove salts
                mol = Chem.MolFromSmiles(smile)
                if mol:
                    mol = remover.StripMol(mol)
                    Chem.SanitizeMol(mol)
                    valid_mols.append(mol)
                    valid_smiles.append(Chem.MolToSmiles(mol))
            except Exception as e:
                logger.warning(f"Invalid SMILES {smile}: {str(e)}")
                continue

        if not valid_mols:
            raise ValueError("No valid molecules remaining after preprocessing")

        # Generate fingerprints
        fingerprints = []
        for mol in valid_mols:
            from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
            fp = GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            arr = np.zeros((1,))
            AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
            fingerprints.append(arr)

        # Calculate descriptors
        descriptor_data = [
            self._get_descriptors(mol) for mol in valid_mols
        ]

        # Combine features
        X = np.hstack([np.array(fingerprints), np.array(descriptor_data)])
        return X, valid_smiles

    def _get_descriptors(self, mol):
        """Calculate molecular descriptors"""
        return [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.RingCount(mol)
        ]
