from src.pipeline import DrugDiscoveryPipeline
import sys

def format_admet(properties):
    return f"""
- Molecular Weight: {properties.get('MolWt', 'N/A'):.2f}
- LogP: {properties.get('LogP', 'N/A'):.2f}
- TPSA: {properties.get('TPSA', 'N/A'):.2f}
- H-Bond Donors: {properties.get('NumHDonors', 'N/A')}
- H-Bond Acceptors: {properties.get('NumHAcceptors', 'N/A')}
- Rotatable Bonds: {properties.get('NumRotatableBonds', 'N/A')}
- Rings: {properties.get('RingCount', 'N/A')}
- Lipinski HBD: {'Pass' if properties.get('HBD_Lipinski') else 'Fail'}
- Lipinski HBA: {'Pass' if properties.get('HBA_Lipinski') else 'Fail'}"""

if __name__ == "__main__":
    print("\n" + "="*40)
    print("DRUG DISCOVERY WORKFLOW")
    print("="*40)

    protein_name = input("\nEnter target protein name (e.g., 'EGFR', 'DRD2'): ").strip()

    try:
        pipeline = DrugDiscoveryPipeline()
        candidates = pipeline.run_pipeline(protein_name)

        if candidates:
            print("\nTOP CANDIDATES WITH ADMET PROPERTIES:")
            for i, candidate in enumerate(candidates, 1):
                print(f"\n{i}. {candidate['SMILES']}")
                print(f" Predicted IC50: {candidate['Predicted_IC50']:.2f} nM")
                print(" ADMET Properties:" + format_admet(candidate))
        else:
            print("\nNo valid candidates generated")

    except Exception as e:
        print(f"\nERROR: {str(e)}")
        print("\nTroubleshooting:")
        print("- Try different protein names (e.g., 'DRD2', 'ACHE')")
        print("- Check internet connection for ChEMBL API access")
        print("- Ensure RDKit and other dependencies are properly installed")
        sys.exit(1)
