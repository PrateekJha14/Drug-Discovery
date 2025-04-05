// src/components/ResultsPage.jsx
import { useLocation, useNavigate } from "react-router-dom";
import { useContext } from "react";
import { useState } from "react";
import { ThemeContext } from "../context/ThemeContext";


export default function ResultsPage() {
  const { state } = useLocation();
  const navigate = useNavigate();
  const results = state?.results;
  const [sortBy, setSortBy] = useState("ic50");
  const [filterSafety, setFilterSafety] = useState("all");
  const { theme, toggleTheme } = useContext(ThemeContext);


  if (!results) {
    return (
      <div className="app">
        <div className="spinner-container">
          <p>No results found. Please return to the home page.</p>
          <button 
            className="action-button primary"
            onClick={() => navigate("/")}
          >
            Go Home
          </button>
        </div>
      </div>
    );
  }

  // Sort candidates
  const sortedCandidates = [...results.candidates].sort((a, b) => {
    if (sortBy === "ic50") {
      return a.predicted_ic50 - b.predicted_ic50;
    } else if (sortBy === "admet") {
      return b.admet_score - a.admet_score;
    } 
    return 0;
  });

  // Filter candidates
  const filteredCandidates = sortedCandidates.filter(candidate => {
    if (filterSafety === "all") return true;
    if (filterSafety === "safe") {
      return !candidate.properties.AMES_Toxicity && 
             !candidate.properties.hERG_Inhibition && 
             !candidate.properties.HepatotoxicityRisk;
    }
    if (filterSafety === "risk") {
      return candidate.properties.AMES_Toxicity || 
             candidate.properties.hERG_Inhibition || 
             candidate.properties.HepatotoxicityRisk;
    }
    return true;
  });

  return (
    <div className="app">
      <div className="results-container">
      <div className="results-header">
  <h2>Results for {results.query}</h2>
  <div className="header-actions">
    <div className="theme-toggle" style={{margin : 25, padding : 0}} onClick={toggleTheme}>
      <div className={`theme-option ${theme === 'light' ? 'active' : ''}`}>
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
          <circle cx="12" cy="12" r="5"></circle>
          <line x1="12" y1="1" x2="12" y2="3"></line>
          <line x1="12" y1="21" x2="12" y2="23"></line>
          <line x1="4.22" y1="4.22" x2="5.64" y2="5.64"></line>
          <line x1="18.36" y1="18.36" x2="19.78" y2="19.78"></line>
          <line x1="1" y1="12" x2="3" y2="12"></line>
          <line x1="21" y1="12" x2="23" y2="12"></line>
          <line x1="4.22" y1="19.78" x2="5.64" y2="18.36"></line>
          <line x1="18.36" y1="5.64" x2="19.78" y2="4.22"></line>
        </svg>
      </div>
      <div className={`theme-option ${theme === 'dark' ? 'active' : ''}`}>
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
          <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"></path>
        </svg>
      </div>
    </div>
    <button className="back-button" onClick={() => navigate("/")}>
      <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
        <line x1="19" y1="12" x2="5" y2="12"></line>
        <polyline points="12 19 5 12 12 5"></polyline>
      </svg>
      Back to Search
    </button>
  </div>
</div>


        <div className="summary">
          <p>
            <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <path d="M22 12h-4l-3 9L9 3l-3 9H2"></path>
            </svg>
            Found <strong>{results.count}</strong> potential compounds
          </p>
          <p>
            <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <circle cx="12" cy="12" r="10"></circle>
              <polyline points="12 6 12 12 16 14"></polyline>
            </svg>
            Processing time: <strong>{results.processing_time.toFixed(2)}s</strong>
          </p>
        </div>

        <div className="filtering">
          <div className="filter-group">
            <span>Safety:</span>
            <button 
              className={`filter-button ${filterSafety === 'all' ? 'active' : ''}`}
              onClick={() => setFilterSafety('all')}
            >
              All
            </button>
            <button 
              className={`filter-button ${filterSafety === 'safe' ? 'active' : ''}`}
              onClick={() => setFilterSafety('safe')}
            >
              Safe Only
            </button>
            <button 
              className={`filter-button ${filterSafety === 'risk' ? 'active' : ''}`}
              onClick={() => setFilterSafety('risk')}
            >
              Has Risks
            </button>
          </div>
          
          <div className="filter-group">
            <span>Sort by:</span>
            <select 
              className="sort-select"
              value={sortBy}
              onChange={(e) => setSortBy(e.target.value)}
            >
              <option value="ic50">Predicted IC50</option>
              <option value="admet">ADMET Score</option>
              {/* <option value="mw">Molecular Weight</option> */}
            </select>
          </div>
        </div>

        <div className="candidates-grid">
          {filteredCandidates.map((candidate, index) => (
            <div key={index} className="candidate-card">
              <div className="card-header">
                <h3>Candidate #{index + 1}</h3>
                <span className="card-badge">
                  {candidate.admet_score > 0.7 ? "High" : candidate.admet_score > 0.5 ? "Medium" : "Low"} Potential
                </span>
              </div>
              
              <div className="card-body">
                <div className="molecule-container">
                  <img
                    src={`https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(
                      candidate.smiles
                    )}/image`}
                    alt={"unable to generate molecule image"}
                    className="molecule-image"
                  />
                </div>
                
                <div className="smiles-code">{candidate.smiles}</div>
                
                <div className="prediction-score">
                  <div className="score-item">
                    <span className="score-value">{candidate.predicted_ic50.toFixed(2)}</span>
                    <span className="score-label">Predicted IC50</span>
                  </div>
                  <div className="score-item">
                    <span className="score-value">{candidate.admet_score.toFixed(2)}</span>
                    <span className="score-label">ADMET Score</span>
                  </div>
                </div>
                
                <div className="property-section">
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
                
                <div className="property-section safety-section">
                  <h4>Safety Profile</h4>
                  <div className="safety-grid">
                    <div
                      className={`safety-item ${
                        candidate.properties.AMES_Toxicity ? "bad" : "good"
                      }`}
                    >
                      {candidate.properties.AMES_Toxicity ? "AMES Toxic" : "AMES Safe"}
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
                
                {/* <div className="card-actions">
                  <button className="action-button secondary">
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                      <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
                      <polyline points="7 10 12 15 17 10"></polyline>
                      <line x1="12" y1="15" x2="12" y2="3"></line>
                    </svg>
                    Export Data
                  </button>
                  <button className="action-button primary">
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                      <circle cx="11" cy="11" r="8"></circle>
                      <line x1="21" y1="21" x2="16.65" y2="16.65"></line>
                    </svg>
                    Analyze Further
                  </button>
                </div> */}
              </div>
            </div>
          ))}
        </div>
        
        {filteredCandidates.length === 0 && (
          <div style={{ textAlign: 'center', padding: '3rem 0' }}>
            <p>No compounds match the current filters.</p>
            <button 
              className="action-button primary" 
              style={{ margin: '1rem auto' }}
              onClick={() => setFilterSafety('all')}
            >
              Reset Filters
            </button>
          </div>
        )}
      </div>
    </div>
  );
}