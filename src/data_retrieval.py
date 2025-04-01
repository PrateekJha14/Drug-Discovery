from chembl_webresource_client.new_client import new_client
import pandas as pd
from tqdm import tqdm
import logging

class AdvancedDataRetriever:
    def __init__(self, confidence_threshold=0.7, data_quality_filter=True):
        """
        Advanced data retrieval with enhanced target and activity selection
        
        Args:
            confidence_threshold (float): Confidence score for target selection
            data_quality_filter (bool): Apply additional data quality filters
        """
        self.client = new_client
        self.logger = logging.getLogger(__name__)
        self.confidence_threshold = confidence_threshold
        self.data_quality_filter = data_quality_filter
        
    def get_target_data(self, protein_name):
        """Fetch ChEMBL ID and bioactivity data with more flexible retrieval"""
        print(f"\n{'='*30}\nSearching for target: {protein_name}")
        
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
                    
                    if len(df) > 0:
                        print(f"\nRetrieved {len(df)} valid bioactivity records")
                        return df
                        
        raise ValueError(f"No bioactivity data found for {protein_name}")
    
    def _find_targets(self, protein_name):
        """Find targets with multiple search strategies"""
        search_strategies = [
            {'target_synonym__iexact': protein_name, 'organism': 'Homo sapiens'},
            {'pref_name__icontains': protein_name, 'organism': 'Homo sapiens'},
            {'pref_name__icontains': protein_name}
        ]
        
        for strategy in search_strategies:
            targets = self.client.target.filter(**strategy)
            if targets:
                return targets
                
        return []
        
    def _select_best_target(self, targets):
        """Select best target based on multiple criteria"""
        targets = sorted(targets, key=lambda x: (
            len(x.get('target_components', [])), # Component count
            x.get('confidence_score', 0), # Confidence score
            len(str(x.get('pref_name', ''))) # Name length
        ), reverse=True)
        
        best_target = targets[0]
        self.logger.info(f"Selected target: {best_target.get('pref_name')}")
        return best_target
        
    def _fetch_bioactivities(self, target):
        """Fetch bioactivity data with advanced filtering"""
        target_id = target['target_chembl_id']
        activities = list(tqdm(
            self.client.activity.filter(
                target_chembl_id=target_id,
                standard_type__in=["IC50", "Ki", "Kd"],
                relation="=",
                assay_type="B",
                standard_units__in=["nM", "μM", "pM", "mM"]
            ).only(
                'molecule_chembl_id', 'canonical_smiles',
                'standard_value', 'standard_type',
                'standard_units', 'pchembl_value',
                'confidence_score'
            ),
            desc="Downloading records"
        ))
        
        return activities
        
    def _process_activities(self, activities):
        """Process and filter bioactivity data"""
        df = pd.DataFrame(activities)
        
        # Comprehensive data cleaning
        df = df[
            df['standard_value'].notna() &
            df['canonical_smiles'].notna() &
            (df['confidence_score'] >= self.confidence_threshold if self.data_quality_filter else True)
        ]
        
        df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
        df = df.dropna(subset=['standard_value'])
        
        # Convert to nM
        df = self._convert_to_nM(df)
        
        return df
        
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

# Add compatibility layer for existing code
class DataRetriever(AdvancedDataRetriever):
    def __init__(self):
        super().__init__()
