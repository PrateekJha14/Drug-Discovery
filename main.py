
import sys

from backend.src.pipeline import DrugDiscoveryPipeline

def format_admet(properties):
    basic = f"""
- Molecular Weight: {properties.get('MolWt', 'N/A'):.2f}
- LogP: {properties.get('LogP', 'N/A'):.2f}
- TPSA: {properties.get('TPSA', 'N/A'):.2f}
- H-Bond Donors: {properties.get('NumHDonors', 'N/A')}
- H-Bond Acceptors: {properties.get('NumHAcceptors', 'N/A')}
- Rotatable Bonds: {properties.get('NumRotatableBonds', 'N/A')}
- Rings: {properties.get('RingCount', 'N/A')}
- QED: {properties.get('QED', 'N/A'):.2f}
- Lipinski Violations: {properties.get('Lipinski_Violations', 'N/A')}
- Veber Violations: {properties.get('Veber_Violations', 'N/A')}
- Druglikeness: {properties.get('Druglikeness', 'N/A')}"""

    adme = f"""
- Water Solubility: {properties.get('WaterSolubility', 'N/A')}
- BBB Penetration: {'Yes' if properties.get('BBB_Penetration') else 'No'}
- Human Intestinal Absorption: {'High' if properties.get('HIA') else 'Low'}
- Caco-2 Permeability: {'High' if properties.get('Caco2') else 'Low'}
- P-gp Substrate: {'Yes' if properties.get('Pgp_Substrate') else 'No'}
- P-gp Inhibitor: {'Yes' if properties.get('Pgp_Inhibitor') else 'No'}
- Plasma Protein Binding: {'High' if properties.get('PlasmaProteinBinding') else 'Low'}"""

    metabolism = f"""
- CYP1A2 Inhibitor: {'Yes' if properties.get('CYP1A2_Inhibitor') else 'No'}
- CYP2C9 Inhibitor: {'Yes' if properties.get('CYP2C9_Inhibitor') else 'No'}
- CYP2C19 Inhibitor: {'Yes' if properties.get('CYP2C19_Inhibitor') else 'No'}
- CYP2D6 Inhibitor: {'Yes' if properties.get('CYP2D6_Inhibitor') else 'No'}
- CYP3A4 Inhibitor: {'Yes' if properties.get('CYP3A4_Inhibitor') else 'No'}
- CYP Promiscuity: {'High' if properties.get('CYP_Promiscuity') else 'Low'}"""

    toxicity = f"""
- AMES Toxicity: {'Yes' if properties.get('AMES_Toxicity') else 'No'}
- hERG Inhibition: {'Yes' if properties.get('hERG_Inhibition') else 'No'}
- Carcinogenicity: {'Yes' if properties.get('Carcinogenicity') else 'No'}
- Acute Oral Toxicity: {properties.get('AcuteOralToxicity', 'N/A')}
- Hepatotoxicity Risk: {'Yes' if properties.get('HepatotoxicityRisk') else 'No'}"""

    overall = f"""
- ADMET Score: {properties.get('ADMET_Score', 'N/A'):.2f}"""

    return f"{basic}\n\nADME Properties:{adme}\n\nMetabolism:{metabolism}\n\nToxicity:{toxicity}\n\nOverall:{overall}"

def format_multi_target_results(candidate):
    """Format multi-target prediction results"""
    targets_info = ""
    if 'targets' in candidate:
        targets_info = "\n Target Specificity:"
        for target, weight in zip(candidate.get('targets', []), candidate.get('target_weights', [])):
            targets_info += f"\n  - {target}: {weight:.2f}"
    
    return targets_info

if __name__ == "__main__":
    print("\n" + "="*40)
    print("MULTI-TARGET DRUG DISCOVERY WORKFLOW")
    print("="*40)

    protein_name = input("\nEnter target protein name (e.g., 'EGFR', 'DRD2'): ").strip()

    try:
        pipeline = DrugDiscoveryPipeline()
        candidates = pipeline.run_pipeline(protein_name)

        if candidates:
            print("\nTOP CANDIDATES WITH MULTI-TARGET PROPERTIES:")
            for i, candidate in enumerate(candidates, 1):
                print(f"\n{i}. {candidate['SMILES']}")
                print(f" Predicted IC50: {candidate['Predicted_IC50']:.2f} nM")
                print(" ADMET Properties:" + format_admet(candidate))
                print(format_multi_target_results(candidate))
        else:
            print("\nNo valid candidates generated")

    except Exception as e:
        print(f"\nERROR: {str(e)}")
        print("\nTroubleshooting:")
        print("- Try different protein names (e.g., 'DRD2', 'ACHE')")
        print("- Check internet connection for ChEMBL and STRING API access")
        print("- Ensure RDKit and other dependencies are properly installed")
        sys.exit(1)
