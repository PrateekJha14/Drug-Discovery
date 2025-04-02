import numpy as np
import logging
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, MolSurf
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Fragments import fr_benzene, fr_halogen

logger = logging.getLogger(__name__)

class EnhancedADMETPredictor:
    def __init__(self):
        # Basic properties
        self.basic_properties = [
            'MolWt', 'LogP', 'TPSA', 'NumHDonors', 'NumHAcceptors',
            'NumRotatableBonds', 'RingCount', 'AromaticRings'
        ]
        
        # ADME properties
        self.adme_properties = [
            'HIA', 'Caco2', 'BBB_Penetration', 'Pgp_Substrate', 'Pgp_Inhibitor',
            'WaterSolubility', 'PlasmaProteinBinding'
        ]
        
        # Metabolism properties
        self.metabolism_properties = [
            'CYP1A2_Inhibitor', 'CYP2C9_Inhibitor', 'CYP2C19_Inhibitor', 
            'CYP2D6_Inhibitor', 'CYP3A4_Inhibitor', 'CYP_Promiscuity'
        ]
        
        # Toxicity properties
        self.toxicity_properties = [
            'AMES_Toxicity', 'hERG_Inhibition', 'Carcinogenicity', 
            'AcuteOralToxicity', 'HepatotoxicityRisk'
        ]
        
        # Weights for ADMET scoring (based on admetSAR research)
        self.property_weights = {
            'AMES_Toxicity': 0.6021,
            'AcuteOralToxicity': 0.5867,
            'Caco2': 0.3074,
            'Carcinogenicity': 0.6857,
            'CYP1A2_Inhibitor': 0.2379,
            'CYP2C19_Inhibitor': 0.2685,
            'CYP2C9_Inhibitor': 0.2966,
            'CYP2D6_Inhibitor': 0.2829,
            'CYP3A4_Inhibitor': 0.3421,
            'CYP_Promiscuity': 0.4808,
            'hERG_Inhibition': 0.4059,
            'HIA': 0.7548,
            'Pgp_Inhibitor': 0.3735,
            'Pgp_Substrate': 0.3464
        }

    def calculate_admet(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None
        

        properties = self._calculate_basic_properties(mol)
        properties.update(self._calculate_adme_properties(mol))
        properties.update(self._calculate_metabolism_properties(mol))
        properties.update(self._calculate_toxicity_properties(mol))
        
        # Lipinski Rule of 5 compliance
        properties.update(self._calculate_druglikeness(mol))
        
        properties['SMILES'] = smiles
        
        # Calculating ADMET score
        properties['ADMET_Score'] = self.calculate_admet_score(properties)
        
        return properties
    
    def _calculate_basic_properties(self, mol):

        return {
            'MolWt': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'TPSA': Descriptors.TPSA(mol),
            'NumHDonors': Descriptors.NumHDonors(mol),
            'NumHAcceptors': Descriptors.NumHAcceptors(mol),
            'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
            'RingCount': Descriptors.RingCount(mol),
            'AromaticRings': Lipinski.NumAromaticRings(mol),
            'QED': QED.qed(mol)  
        }
    
    def _calculate_adme_properties(self, mol):
        """Calculate ADME properties using simple rule-based models."""
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        mw = Descriptors.MolWt(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        
        # Human Intestinal Absorption (HIA)
        hia_score = 1 if (tpsa < 140 and h_donors <= 5) else 0
        
        # Caco-2 permeability
        caco2_score = 1 if (tpsa < 100 and logp > 0 and logp < 5) else 0
        
        # Blood-Brain Barrier penetration
        bbb_score = 1 if (logp > 0 and logp < 3 and tpsa < 90 and mw < 450) else 0
        
        # P-glycoprotein substrate
        pgp_substrate = 1 if (mw > 400 and logp > 4) else 0
        
        # P-glycoprotein inhibitor
        pgp_inhibitor = 1 if (logp > 4.5 and mw > 400) else 0
        
        # Water solubility
        if logp < 0:
            water_solubility = "High"
        elif logp < 3:
            water_solubility = "Moderate"
        else:
            water_solubility = "Low"
        
        # Plasma protein binding
        ppb = 1 if (logp > 3 and mw > 400) else 0
        
        return {
            'HIA': hia_score,
            'Caco2': caco2_score,
            'BBB_Penetration': bbb_score,
            'Pgp_Substrate': pgp_substrate,
            'Pgp_Inhibitor': pgp_inhibitor,
            'WaterSolubility': water_solubility,
            'PlasmaProteinBinding': ppb
        }
    
    def _calculate_metabolism_properties(self, mol):
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        
        # CYP inhibition models (simplified)
        cyp1a2_inhibitor = 1 if (aromatic_rings >= 2 and logp > 2.5) else 0
        cyp2c9_inhibitor = 1 if (logp > 3.5 and mw > 220) else 0
        cyp2c19_inhibitor = 1 if (logp > 3 and mw > 300) else 0
        cyp2d6_inhibitor = 1 if (Descriptors.NumHDonors(mol) >= 1 and logp > 2) else 0
        cyp3a4_inhibitor = 1 if (logp > 3.5 and mw > 300) else 0
        
        # CYP promiscuity (inhibits multiple CYP isoforms)
        inhibition_count = sum([cyp1a2_inhibitor, cyp2c9_inhibitor, 
                               cyp2c19_inhibitor, cyp2d6_inhibitor, 
                               cyp3a4_inhibitor])
        cyp_promiscuity = 1 if inhibition_count >= 3 else 0
        
        return {
            'CYP1A2_Inhibitor': cyp1a2_inhibitor,
            'CYP2C9_Inhibitor': cyp2c9_inhibitor,
            'CYP2C19_Inhibitor': cyp2c19_inhibitor,
            'CYP2D6_Inhibitor': cyp2d6_inhibitor,
            'CYP3A4_Inhibitor': cyp3a4_inhibitor,
            'CYP_Promiscuity': cyp_promiscuity
        }
    
    def _calculate_toxicity_properties(self, mol):
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        
        # Check for structural alerts
        has_halogen = sum([atom.GetAtomicNum() in [9, 17, 35, 53] for atom in mol.GetAtoms()]) > 0
        
        # AMES mutagenicity 
        ames_toxic = 1 if (has_halogen and logp > 3 and mw > 250) else 0
        
        # hERG inhibition (cardiotoxicity)
        herg_inhibition = 1 if (logp > 3.7 and Descriptors.NumHDonors(mol) > 0) else 0
        
        # Carcinogenicity 
        carcinogenic = 1 if (ames_toxic == 1 and mw > 300) else 0
        
        # Acute oral toxicity 
        if logp > 5 or mw > 600:
            acute_oral_toxicity = "High"
        elif logp > 3 or mw > 400:
            acute_oral_toxicity = "Moderate"
        else:
            acute_oral_toxicity = "Low"
            
        # Hepatotoxicity risk
        hepatotoxicity = 1 if (logp > 3 and mw > 300) else 0
        
        return {
            'AMES_Toxicity': ames_toxic,
            'hERG_Inhibition': herg_inhibition,
            'Carcinogenicity': carcinogenic,
            'AcuteOralToxicity': acute_oral_toxicity,
            'HepatotoxicityRisk': hepatotoxicity
        }
    
    def _calculate_druglikeness(self, mol):
        """Calculate Lipinski's Rule of 5 and other druglikeness parameters."""
        # Lipinski's Rule of 5
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if h_donors > 5: violations += 1
        if h_acceptors > 10: violations += 1
        
        # Veber rules
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)
        veber_violations = 0
        if rotatable_bonds > 10: veber_violations += 1
        if tpsa > 140: veber_violations += 1
        
        return {
            'Lipinski_Violations': violations,
            'Veber_Violations': veber_violations,
            'Druglikeness': "High" if violations == 0 else "Medium" if violations == 1 else "Low"
        }
    
    def calculate_admet_score(self, properties):
        """Calculate a comprehensive ADMET score based on weighted properties."""
        score = 0
        total_weight = 0
        
        # Transforming categorical values to binary
        binary_props = properties.copy()
        if 'AcuteOralToxicity' in binary_props:
            binary_props['AcuteOralToxicity'] = 1 if binary_props['AcuteOralToxicity'] == 'Low' else 0
        if 'WaterSolubility' in binary_props:
            binary_props['WaterSolubility'] = 1 if binary_props['WaterSolubility'] in ['Moderate', 'High'] else 0
            
        # Calculating weighted score
        for prop, weight in self.property_weights.items():
            if prop in binary_props:
                # For toxicity endpoints, we want 0 (non-toxic) to contribute positively
                if prop in ['AMES_Toxicity', 'Carcinogenicity', 'hERG_Inhibition', 
                           'CYP1A2_Inhibitor', 'CYP2C9_Inhibitor', 'CYP2C19_Inhibitor', 
                           'CYP2D6_Inhibitor', 'CYP3A4_Inhibitor', 'CYP_Promiscuity',
                           'Pgp_Inhibitor']:
                    score += (1 - binary_props[prop]) * weight
                else:
                    score += binary_props[prop] * weight
                total_weight += weight
        
        # Normalizing score between 0 and 1
        if total_weight > 0:
            normalized_score = score / total_weight
        else:
            normalized_score = 0
            
        # Adjust based on Lipinski violations
        if 'Lipinski_Violations' in properties:
            if properties['Lipinski_Violations'] > 1:
                normalized_score *= (1 - 0.1 * properties['Lipinski_Violations'])
        
        return max(0, min(1, normalized_score))
    
    def analyze_candidates(self, candidates):
        """Computing ADMET properties for top drug candidates."""
        results = []
        for smiles, ic50 in candidates:
            admet_props = self.calculate_admet(smiles)
            if admet_props:
                admet_props['Predicted_IC50'] = ic50
                results.append(admet_props)
        

        if results:
            results.sort(key=lambda x: x.get('ADMET_Score', 0), reverse=True)
            
        return results
