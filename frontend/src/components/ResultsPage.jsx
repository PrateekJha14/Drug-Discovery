import { useLocation } from "react-router-dom";

export default function ResultsPage() {
  const { state } = useLocation();
  const results = state?.results;

  return (
    <div className="results-container">
      <h2>Results for {results?.query}</h2>
      <div className="summary">
        <p>
          Found {results?.count} compounds in{" "}
          {results?.processing_time.toFixed(2)}s
        </p>
      </div>

      <div className="candidates-grid">
        {console.log(results?.candidates)}
        {results?.candidates.map((candidate, index) => (
          <div key={index} className="candidate-card">
            <h3>Candidate #{index + 1}</h3>
            <div className="card-section">
              <h4>Basic Info</h4>
              <div className="card-section">
                <h4>Molecular Structure</h4>
                <img
                  src={`https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(
                    candidate.smiles
                  )}/image`}
                  alt={`Molecular structure of ${candidate.smiles}`}
                  className="molecule-image"
                />
              </div>
              <p>
                <strong>SMILES:</strong> {candidate.smiles}
              </p>
              <p>
                <strong>Predicted IC50:</strong>{" "}
                {candidate.predicted_ic50.toFixed(2)}
              </p>
              <p>
                <strong>ADMET Score:</strong> {candidate.admet_score.toFixed(2)}
              </p>
            </div>

            <div className="card-section">
              <h4>Key Properties</h4>
              <div className="property-grid">
                <div className="property">
                  <span>Molecular Weight</span>
                  <span>{candidate.properties.MolWt.toFixed(2)}</span>
                </div>
                <div className="property">
                  <span>LogP</span>
                  <span>{candidate.properties.LogP.toFixed(2)}</span>
                </div>
                <div className="property">
                  <span>Water Solubility</span>
                  <span>{candidate.properties.WaterSolubility}</span>
                </div>
                <div className="property">
                  <span>Druglikeness</span>
                  <span>{candidate.properties.Druglikeness}</span>
                </div>
              </div>
            </div>

            <div className="card-section">
              <h4>Safety Profile</h4>
              <div className="safety-grid">
                <div
                  className={`safety-item ${
                    candidate.properties.AMES_Toxicity ? "bad" : "good"
                  }`}
                >
                  {candidate.properties.AMES_Toxicity ? "Toxic" : "Non-toxic"}
                </div>
                <div
                  className={`safety-item ${
                    candidate.properties.hERG_Inhibition ? "bad" : "good"
                  }`}
                >
                  {candidate.properties.hERG_Inhibition
                    ? "hERG Risk"
                    : "hERG Safe"}
                </div>
                <div
                  className={`safety-item ${
                    candidate.properties.HepatotoxicityRisk ? "bad" : "good"
                  }`}
                >
                  {candidate.properties.HepatotoxicityRisk
                    ? "Liver Risk"
                    : "Liver Safe"}
                </div>
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}
