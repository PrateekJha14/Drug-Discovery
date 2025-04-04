from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import rdFingerprintGenerator
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

        # fingerprints
        fingerprints = []
        for mol in valid_mols:
            morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            fp = morgan_gen.GetFingerprint(mol)
            arr = np.zeros((1,))
            AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
            fingerprints.append(arr)

        # Calculating descriptors
        descriptor_data = [
            self._get_descriptors(mol) for mol in valid_mols
        ]

        X = np.hstack([np.array(fingerprints), np.array(descriptor_data)])
        return X, valid_smiles

    def _get_descriptors(self, mol):
        return [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.RingCount(mol)
        ]

    def process_single_molecule(self, mol):
        if not mol:
            return None, None

        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fp = morgan_gen.GetFingerprint(mol)
        arr = np.zeros((1,))
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
        
        descriptors = self._get_descriptors(mol)
        
        X = np.hstack([arr, np.array(descriptors)])
        
        return X, Chem.MolToSmiles(mol)
