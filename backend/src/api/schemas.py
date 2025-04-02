from pydantic import BaseModel, Field, validator
from typing import List, Dict, Optional, Any, Union

class ProteinRequest(BaseModel):
    protein_name: str = Field(..., description="Name of the protein target (e.g., 'DRD2', 'EGFR')")
    
    @validator('protein_name')
    def protein_name_must_not_be_empty(cls, v):
        if not v or not v.strip():
            raise ValueError('Protein name cannot be empty')
        return v.strip()

class MoleculeRequest(BaseModel):
    smiles: str = Field(..., description="SMILES string of the molecule")
    
    @validator('smiles')
    def smiles_must_not_be_empty(cls, v):
        if not v or not v.strip():
            raise ValueError('SMILES string cannot be empty')
        return v.strip()

class ADMETProperty(BaseModel):
    name: str
    value: Union[str, float, int, bool]
    category: str

class ADMETResponse(BaseModel):
    smiles: str
    predicted_ic50: float
    admet_score: float
    properties: Dict[str, Any]

class CandidateResponse(BaseModel):
    smiles: str
    predicted_ic50: float
    admet_score: float
    properties: Dict[str, Any]
    targets: Optional[List[str]] = None
    target_weights: Optional[List[float]] = None

class DiscoveryResponse(BaseModel):
    candidates: List[CandidateResponse]
    count: int
    query: str
    processing_time: float
