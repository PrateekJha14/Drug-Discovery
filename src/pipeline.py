import pandas as pd
from rdkit import Chem
from src.admet import  EnhancedADMETPredictor
from src.data_retrieval import AdvancedDataRetriever
from src.preprocessing import Preprocessor
from src.modeling import DrugModel
from src.generation import MoleculeGenerator
import logging

logging.basicConfig(level=logging.INFO)

class DrugDiscoveryPipeline:
    def __init__(self):
        self.model = DrugModel()
        self.admet_predictor = EnhancedADMETPredictor()
        self.retriever = AdvancedDataRetriever()
        self.preprocessor = Preprocessor()
        self.generator = MoleculeGenerator()

    def retrieve_data(self, protein_name):
        # return self.retriever.get_target_data(protein_name)
        return self.retriever.run_integrated_search(protein_name)

    def preprocess_data(self, df):
        X, valid_smiles = self.preprocessor.process_smiles(df['canonical_smiles'])
        y = df['standard_value'].values
        return X, y, valid_smiles

    def train_model(self, X, y):
        self.model.train(X, y)

    def generate_molecules(self, num_samples=100, training_data=None):
        return self.generator.generate_molecules(n_molecules=num_samples, training_data=training_data)

    def filter_top_candidates(self, predictions, threshold=1.0):
        return [p for p in predictions if p[1] < threshold]

    def analyze_candidates(self, candidates):
        """Analyze ADMET properties of selected candidates and return full data"""
        analyzed = self.admet_predictor.analyze_candidates(candidates)
        return [c for c in analyzed if c is not None]

    def run_pipeline(self, protein_name):
        print(f"\n{'='*30}\nDRUG DISCOVERY PIPELINE: {protein_name.upper()}\n{'='*30}")

        try:
            print("\n[1/5] Retrieving bioactivity data...")
            df = self.retrieve_data(protein_name)
            print("="*40)
            print(df)
            print("="*40)

            if df.empty:
                raise ValueError(f"No bioactivity data found for {protein_name}")
            print(f"Retrieved {len(df)} records")

            # required_columns = ['canonical_smiles', 'standard_value']
            # if not all(col in df.columns for col in required_columns):
            #     raise ValueError("Missing required columns in bioactivity data")

            print("\n[2/5] Preprocessing molecules...")
            X, y, valid_smiles = self.preprocess_data(df)

            if X.shape[0] < 1:
                raise ValueError(f"Only {X.shape[0]} samples available - need at least 50 for modeling")

            print("\n[3/5] Training prediction model...")
            self.train_model(X, y)

            print("\n[4/5] Generating novel molecules...")
            generated_smiles = self.generate_molecules(
                num_samples=100,
                training_data=valid_smiles # Pass preprocessed SMILES from your dataset
            )

            print("\n[5/5] Screening candidates...")
            X_new, valid_new = self.preprocessor.process_smiles(generated_smiles)

            if X_new.shape[0] == 0:
                raise ValueError("No valid molecules generated")

            predictions = []
            for smile, features in zip(valid_new, X_new):
                try:
                    ic50 = self.model.predict(features.reshape(1, -1))[0]
                    predictions.append((smile, ic50))
                except Exception as e:
                    logging.warning(f"Prediction failed for {smile}: {str(e)}")

            results = sorted(predictions, key=lambda x: x[1])
            n_top = max(5, len(results) // 20)
            top_candidates = results[:n_top]

            return self.analyze_candidates(top_candidates)

        except Exception as e:
            logging.error(f"Pipeline failed: {str(e)}")
            raise
