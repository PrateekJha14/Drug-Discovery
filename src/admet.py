import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging

logger = logging.getLogger(__name__)

class ADMETPredictor:
    def __init__(self):
        self.properties = [
            'MolWt', 'LogP', 'TPSA', 'NumHDonors', 'NumHAcceptors',
            'NumRotatableBonds', 'RingCount', 'HBD_Lipinski', 'HBA_Lipinski'
        ]

    def calculate_admet(self, smiles):
        """Calculate ADMET properties for a given molecule."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None

        return {
            'SMILES': smiles,
            'MolWt': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'TPSA': Descriptors.TPSA(mol),
            'NumHDonors': Descriptors.NumHDonors(mol),
            'NumHAcceptors': Descriptors.NumHAcceptors(mol),
            'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
            'RingCount': Descriptors.RingCount(mol),
            'HBD_Lipinski': int(Descriptors.NumHDonors(mol) <= 5),
            'HBA_Lipinski': int(Descriptors.NumHAcceptors(mol) <= 10)
        }

    def analyze_candidates(self, candidates):
        """Compute ADMET properties for top drug candidates."""
        results = []
        for smiles, ic50 in candidates:
            admet_props = self.calculate_admet(smiles)
            if admet_props:
                admet_props['Predicted_IC50'] = ic50
                results.append(admet_props)
        return results
