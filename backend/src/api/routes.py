from fastapi import APIRouter, HTTPException, BackgroundTasks
from ..pipeline import DrugDiscoveryPipeline
from ..admet import EnhancedADMETPredictor
from .schemas import ProteinRequest, MoleculeRequest, ADMETResponse, DiscoveryResponse, CandidateResponse
import time
import logging
from typing import List, Dict, Any

router = APIRouter()
logger = logging.getLogger(__name__)

# Store recent results in memory cache
results_cache = {}

@router.post("/discover", response_model=DiscoveryResponse)
async def discover_compounds(request: ProteinRequest):
    """Run the drug discovery pipeline for a protein target"""
    start_time = time.time()
    protein_name = request.protein_name
    
    # Check cache
    if protein_name in results_cache:
        cached_result = results_cache[protein_name]
        # Update processing time
        cached_result["processing_time"] = time.time() - start_time
        return cached_result
    
    try:
        pipeline = DrugDiscoveryPipeline()
        
        candidates = pipeline.run_pipeline(protein_name)
        
        
        candidate_responses = []
        for candidate in candidates:
            candidate_responses.append(CandidateResponse(
                smiles=candidate["SMILES"],
                predicted_ic50=candidate["Predicted_IC50"],
                admet_score=candidate.get("ADMET_Score", 0.0),
                properties=candidate,
                targets=candidate.get("targets"),
                target_weights=candidate.get("target_weights")
            ))
        
        response = DiscoveryResponse(
            candidates=candidate_responses,
            count=len(candidate_responses),
            query=protein_name,
            processing_time=time.time() - start_time
        )
        
        # Cache result
        results_cache[protein_name] = response
        
        return response
        
    except Exception as e:
        logger.error(f"Error in discovery pipeline: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Pipeline error: {str(e)}")

@router.post("/admet", response_model=ADMETResponse)
async def predict_admet(request: MoleculeRequest):
    """Predict ADMET properties for a molecule"""
    try:
        predictor = EnhancedADMETPredictor()
        
        properties = predictor.calculate_admet(request.smiles)
        
        if not properties:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        # response
        response = ADMETResponse(
            smiles=request.smiles,
            predicted_ic50=0.0,  
            admet_score=properties.get("ADMET_Score", 0.0),
            properties=properties
        )
        
        return response
        
    except Exception as e:
        logger.error(f"Error in ADMET prediction: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Prediction error: {str(e)}")

@router.get("/status")
async def get_status():
    """Check API status"""
    return {"status": "online", "version": "1.0.0"}
