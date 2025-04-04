// src/components/HomePage.jsx
import { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import axios from 'axios';
import BeatLoader from 'react-spinners/BeatLoader';

export default function HomePage() {
  const [protein, setProtein] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const navigate = useNavigate();

  const handleDiscover = async (e) => {
    e.preventDefault();
    setLoading(true);
    try {
      const response = await axios.post('http://localhost:8000/api/discover', {
        protein_name: protein
      });
      navigate('/results', { state: { results: response.data } });
    } catch (err) {
      setError(err.response?.data?.detail || 'Discovery failed');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="home-container">
      <h1>Drug Discovery Platform</h1>
      <form onSubmit={handleDiscover}>
        <input
          type="text"
          value={protein}
          onChange={(e) => setProtein(e.target.value)}
          placeholder="Enter protein target (e.g., DRD2, EGFR)"
          required
        />
        <button type="submit" disabled={loading}>
          {loading ? <BeatLoader size={8} color="#fff" /> : 'Discover Compounds'}
        </button>
      </form>
      {error && <div className="error">{error}</div>}
    </div>
  );
}