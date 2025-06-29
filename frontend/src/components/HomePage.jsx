// src/components/HomePage.jsx
import { useState, useContext } from 'react';
import { useNavigate } from 'react-router-dom';
import axios from 'axios';
import BeatLoader from 'react-spinners/BeatLoader';
import { ThemeContext } from '../context/ThemeContext';

export default function HomePage() {
  const [protein, setProtein] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const navigate = useNavigate();
  const { theme, toggleTheme } = useContext(ThemeContext);

  const handleDiscover = async (e) => {
    e.preventDefault();
    setLoading(true);
      try {
    const response = await axios.post(
      'https://drug-discovery-ai-9vxo.onrender.com/api/discover',
      { protein_name: protein },
      { withCredentials: true }  // Add this line
    );
    navigate('/results', { state: { results: response.data } });
  } catch (err) {
      setError(err.response?.data?.detail || 'Discovery failed');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="app">
      <header className="app-header">
        <div className="app-logo">
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
            <path d="M4.5 16.5c-1.5 1.26-2 5-2 5s3.74-.5 5-2c.71-.84.7-2.13-.09-2.91a2.18 2.18 0 0 0-2.91-.09z"></path>
            <path d="M12 15l-3-3a22 22 0 0 1 2-3.95A12.88 12.88 0 0 1 22 2c0 5.5-5.27 10-10.5 13-1.05.6-2 1.37-3 2z"></path>
          </svg>
          <h1>MoleculeAI</h1>
        </div>
        
        <div className="theme-toggle" onClick={toggleTheme}>
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
      </header>

      <div className="home-container">
        <h1>AI-Powered Drug Discovery</h1>
        <p className="home-subtitle">
          Enter a protein target to discover potential drug compounds with optimized properties and safety profiles
        </p>
        
        <form className="search-form" onSubmit={handleDiscover}>
          <input
            type="text"
            value={protein}
            onChange={(e) => setProtein(e.target.value)}
            placeholder="Enter protein target (e.g., DRD2, EGFR, BRAF)"
            required
          />
          <button type="submit" disabled={loading}>
            {loading ? (
              <BeatLoader size={8} color="#fff" />
            ) : (
              <>
                <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                  <circle cx="11" cy="11" r="8"></circle>
                  <line x1="21" y1="21" x2="16.65" y2="16.65"></line>
                </svg>
                Discover Compounds
              </>
            )}
          </button>
        </form>
        <p className="home-info">
          
        </p>
        
        {error && <div className="error">{error}</div>}
        
        {/* How It Works Section */}
        <div className="how-it-works-section">
          <h2>How It Works</h2>
          <div className="workflow-steps">
            <div className="workflow-step">
              <div className="step-number">1</div>
              <div className="step-content">
                <h3>Data Retrieval</h3>
                <p>Fetches bioactivity data for the target protein from the ChEMBL database</p>
              </div>
            </div>
            
            <div className="workflow-step">
              <div className="step-number">2</div>
              <div className="step-content">
                <h3>Preprocessing</h3>
                <p>Cleans molecules and generates fingerprints and descriptors</p>
              </div>
            </div>
            
            <div className="workflow-step">
              <div className="step-number">3</div>
              <div className="step-content">
                <h3>Model Training</h3>
                <p>Trains a machine learning model to predict molecule activity</p>
              </div>
            </div>
            
            <div className="workflow-step">
              <div className="step-number">4</div>
              <div className="step-content">
                <h3>Molecule Generation</h3>
                <p>Creates novel drug-like molecules using VAE and fragment-based approaches</p>
              </div>
            </div>
            
            <div className="workflow-step">
              <div className="step-number">5</div>
              <div className="step-content">
                <h3>Candidate Analysis</h3>
                <p>Evaluates and ranks molecules based on predicted activity and ADMET properties</p>
              </div>
            </div>
          </div>
        </div>
        
        {/* <div className="features-grid">
          <div className="feature-card">
            <svg className="feature-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"></path>
              <line x1="12" y1="9" x2="12" y2="13"></line>
              <line x1="12" y1="17" x2="12.01" y2="17"></line>
            </svg>
            <h3 className="feature-title">Advanced AI Models</h3>
            <p>State-of-the-art deep learning models for accurate prediction of binding affinity</p>
          </div>
          
          <div className="feature-card">
            <svg className="feature-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <path d="M9 18l6-6-6-6"></path>
            </svg>
            <h3 className="feature-title">ADMET Scoring</h3>
            <p>Comprehensive evaluation of absorption, distribution, metabolism, excretion, and toxicity</p>
          </div>
          
          <div className="feature-card">
            <svg className="feature-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <circle cx="12" cy="12" r="10"></circle>
              <line x1="2" y1="12" x2="22" y2="12"></line>
              <path d="M12 2a15.3 15.3 0 0 1 4 10 15.3 15.3 0 0 1-4 10 15.3 15.3 0 0 1-4-10 15.3 15.3 0 0 1 4-10z"></path>
            </svg>
            <h3 className="feature-title">Molecular Visualization</h3>
            <p>Interactive structural visualization of potential drug candidates</p>
          </div>
        </div> */}
      </div>
    </div>
  );
}
