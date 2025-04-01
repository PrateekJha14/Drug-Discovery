import logging
from chembl_webresource_client.new_client import new_client
import pandas as pd
import requests
import json

from tqdm import tqdm

class AdvancedDataRetriever:
    def __init__(self, confidence_threshold=0.7, data_quality_filter=True):
        self.client = new_client
        self.logger = logging.getLogger(__name__)
        self.confidence_threshold = confidence_threshold
        self.data_quality_filter = data_quality_filter
        self.string_api_url = "https://string-db.org/api"
        self.string_output_format = "json"
        self.string_version = "11.5"
        
    def get_protein_interactions(self, protein_name, species=9606, score_threshold=700):
        """Query STRING database for protein-protein interactions"""
        params = {
            "identifiers": protein_name,
            "species": species,
            "caller_identity": "drug_discovery_pipeline",
            "network_type": "functional",
            "required_score": score_threshold,
            "add_nodes": "all" 

        }

        # Use the correct endpoint for JSON format
        url = f"{self.string_api_url}/json/network"

        try:
            response = requests.post(url, data=params)
            # print('=/'*30)
            # print(response)
            # print('=/'*30)
            response.raise_for_status()

            # Parse JSON response
            interactions = response.json()

            # Process the JSON response into the expected format
            processed_interactions = []
            for item in interactions:
                    # Extract gene names instead of string IDs
                interactor_id = item.get("stringId_B", "")
                interactor_name = item.get("preferredName_B", interactor_id.split(".")[-1])
                processed_interactions.append({
                    "queryItem": item.get("stringId_A", ""),
                    "interactor": interactor_name, 
                    "score": item.get("score", 0)
                })

            # print('=/'*30)
            # print(processed_interactions)
            # print('=/'*30)


            return processed_interactions
        except Exception as e:
            self.logger.error(f"STRING API error: {str(e)}")
            return []

            
    def get_target_network(self, protein_name):
        """Get protein interaction network and identify potential drug targets"""
        interactions = self.get_protein_interactions(protein_name)
        
        # Sort interactions by confidence score
        interactions = sorted(interactions, key=lambda x: x["score"], reverse=True)
        
        # Get top interactors as potential targets
        potential_targets = [item["interactor"] for item in interactions[:10]]
        
        print(f"Found {len(potential_targets)} potential targets from STRING database")
        return potential_targets
        
    # def run_integrated_search(self, protein_name):
        """Integrated search using both STRING and ChEMBL"""
        # First get protein interaction network from STRING
        potential_targets = self.get_target_network(protein_name)

        
        # Add the original protein to the list
        all_targets = [protein_name] + potential_targets

        print("="*40)
        print(all_targets)
        print("="*40)
        
        # Try each target with ChEMBL
        i=0
        for target in all_targets:
            try:
                print("="*40)
                print(i,target)
                print("="*40)
                print(f"Searching ChEMBL for target: {target} ////////////////////////")
                data = self.get_target_data(target)
                if len(data) > 0:
                    print(f"Found {len(data)} bioactivity records for {target}############")
                    return data
            except Exception as e:
                print("="*40)
                print(i,target)
                print("="*40)
                self.logger.warning(f"Failed to retrieve data for {target}: {str(e)}]]]]]]]]]]]")
                continue
                
        # If we get here, no data was found for any target
        raise ValueError(f"No bioactivity data found for {protein_name} or its interactors")

    def run_integrated_search(self, protein_name):
        """Integrated search using both STRING and ChEMBL to collect data for multiple targets"""
        # Get protein interaction network from STRING
        potential_targets = self.get_target_network(protein_name)

        # Add the original protein to the list
        all_targets = [protein_name] + potential_targets

        print(f"Found {len(potential_targets)} potential targets from STRING database")
        print("="*40)
        print(all_targets)
        print("="*40)

        # Store interaction scores for weighting
        target_weights = {protein_name: 1.0}  # Primary target has weight 1.0

        # Get interaction scores for all targets
        interactions = self.get_protein_interactions(protein_name)
        for interaction in interactions:
            target = interaction["interactor"]
            # Normalize score to 0-1 range (STRING scores are 0-1000)
            score = interaction["score"] / 1000.0
            if target in all_targets:
                target_weights[target] = score

        # Collect data for all targets
        target_data = {}
        for target in all_targets:
            try:
                print(f"Searching ChEMBL for target: {target}")
                data = self.get_target_data(target)
                if len(data) > 0:
                    print(f"Found {len(data)} bioactivity records for {target}")
                    # Add target information to the dataframe
                    data['target_name'] = target
                    data['target_weight'] = target_weights.get(target, 0.5)
                    target_data[target] = data
            except Exception as e:
                self.logger.warning(f"Failed to retrieve data for {target}: {str(e)}")
                continue
            
        if not target_data:
            raise ValueError(f"No bioactivity data found for {protein_name} or its interactors")

        # Combine all target data
        combined_df = pd.concat(target_data.values(), ignore_index=True)
        return combined_df, target_weights

    def get_target_data(self, protein_name):
        """Fetch ChEMBL ID and bioactivity data with more flexible retrieval"""
        print(f"\n{'='*30}\nSearching for target: {protein_name} PPPPPPPPPPPPPPPPP]]]]]]]][[[]]]")
        
        # Expanded search strategies
        search_strategies = [
            # 1. Exact match with human targets
            lambda: self.client.target.filter(
                target_synonym__iexact=protein_name,
                organism='Homo sapiens'
            ),
            # 2. Broader match with human targets
            lambda: self.client.target.filter(
                pref_name__icontains=protein_name,
                organism='Homo sapiens'
            ),
            # 3. Any match across all organisms
            lambda: self.client.target.filter(pref_name__icontains=protein_name),
            # 4. Very broad match
            lambda: self.client.target.filter(target_synonym__icontains=protein_name)
        ]
        
        targets = None
        for strategy in search_strategies:
            targets = strategy()
            if targets:
                break
                
        if not targets:
            raise ValueError(f"No targets found for {protein_name}")
            
        # Sort targets by available data
        targets = sorted(targets,
                        key=lambda x: len(x.get('target_components', [])),
                        reverse=True)
        
        # Try multiple targets if first one fails
        for target in targets[:3]: # Try top 3 targets
            target_id = target['target_chembl_id']
            print(f"Attempting target: {target.get('pref_name')} ({target_id})")
            
            # Expand activity type search
            activity_types = ["IC50", "Ki", "Kd", "EC50"]
            for activity_type in activity_types:
                print(f"Searching with activity type: {activity_type}")
                activities = list(tqdm(
                    self.client.activity.filter(
                        target_chembl_id=target_id,
                        standard_type=activity_type,
                        relation="=",
                        assay_type="B"
                    ).only(
                        'molecule_chembl_id', 'canonical_smiles', 'standard_value',
                        'standard_type', 'standard_units', 'pchembl_value'
                    ),
                    desc="Downloading records"
                ))
                
                if activities:
                    df = pd.DataFrame(activities)
                    
                    # Clean data
                    df = df[df.standard_value.notna() & df.canonical_smiles.notna()]
                    df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
                    df = df.dropna(subset=['standard_value'])
                    
                    # Convert to nM units
                    df = self._convert_to_nM(df)
                    
                    print("[]"*40)
                    print(df)
                    print("[]"*40)

                    if len(df) > 0:
                        print(f"\nRetrieved {len(df)} valid bioactivity records")
                        return df

            
                        
        raise ValueError(f"No bioactivity data found for {protein_name}")

    def _convert_to_nM(self, df):
        """Robust unit conversion to nM"""
        unit_conversion = {
            'μM': 1000,
            'nM': 1,
            'mM': 1000000,
            'pM': 0.001
        }

        df['standard_value'] = df.apply(
            lambda row: row['standard_value'] * unit_conversion.get(row['standard_units'], 1),
            axis=1
        )

        df['standard_units'] = 'nM'
        return df

    def ensembl_to_gene_symbol(self, ensembl_id):
        """Convert Ensembl protein ID to gene symbol"""
        # Extract just the ID part without species prefix
        if "." in ensembl_id:
            ensembl_id = ensembl_id.split(".")[-1]
        
        # For IDs that are already gene symbols
        if not ensembl_id.startswith("ENSP"):
            return ensembl_id
            
        # This would ideally use a proper mapping service
        # For demonstration, return a simplified ID
        return ensembl_id.replace("ENSP", "GENE")
