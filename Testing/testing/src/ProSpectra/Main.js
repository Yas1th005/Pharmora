// src/App.js
import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { ToastContainer, toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';
import MoleculeViewer from './Molecule';
import DescriptorsTable from './Description';
import PharmacologicalProperties from './Properties';
import ExampleMolecules from './Example';

function Main() {
  const [smiles, setSmiles] = useState('');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [exampleMolecules, setExampleMolecules] = useState({});

  useEffect(() => {
    // Fetch example molecules on component mount
    axios.get('http://127.0.0.1:5000/api/examples')
      .then(response => {
        setExampleMolecules(response.data);
      })
      .catch(error => {
        console.error('Error fetching example molecules:', error);
        toast.error('Failed to load example molecules');
      });

      console.log(exampleMolecules)
  }, []);

  const handleSubmit = async (e) => {
    e.preventDefault();
    
    if (!smiles.trim()) {
      toast.error('Please enter a SMILES string');
      return;
    }
    
    setLoading(true);
    
    try {
      const response = await axios.post('http://127.0.0.1:5000/api/analyze', { smiles });
      setResult(response.data);
      toast.success('Molecule analysis complete!');
    } catch (error) {
      console.error('Error analyzing molecule:', error);
      toast.error(error.response?.data?.error || 'Failed to analyze molecule');
      setResult(null);
    } finally {
      setLoading(false);
    }
  };

  const handleExampleSelect = (selectedSmiles) => {
    setSmiles(selectedSmiles);
  };

  return (
    <div className="min-h-screen bg-gray-100">
      <ToastContainer position="top-right" autoClose={3000} />
      
      {/* Header */}
      <header className="bg-blue-700 text-white shadow-lg py-4">
        <div className="container mx-auto px-4">
          <h1 className="text-3xl font-bold">Molecular Analysis Tool</h1>
          <p className="text-blue-100">Analyze molecules using SMILES notation</p>
        </div>
      </header>
      
      {/* Main Content */}
      <main className="container mx-auto px-4 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Left Column - Input Form */}
          <div className="lg:col-span-1">
            <div className="bg-white rounded-lg shadow-md p-6 mb-6">
              <h2 className="text-xl font-semibold mb-4">Enter Molecule</h2>
              <form onSubmit={handleSubmit}>
                <div className="mb-4">
                  <label htmlFor="smiles" className="block text-gray-700 mb-2">SMILES String</label>
                  <textarea
                    id="smiles"
                    className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                    rows="3"
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O for Aspirin"
                  />
                </div>
                <button
                  type="submit"
                  disabled={loading}
                  className="w-full bg-blue-600 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded-md transition duration-300 disabled:opacity-50"
                >
                  {loading ? 'Analyzing...' : 'Analyze Molecule'}
                </button>
              </form>
            </div>
            
            <ExampleMolecules 
              examples={exampleMolecules} 
              onSelect={handleExampleSelect} 
            />
          </div>
          
          {/* Right Column - Results */}
          <div className="lg:col-span-2">
            {loading ? (
              <div className="flex justify-center items-center h-64">
                <div className="animate-spin rounded-full h-16 w-16 border-t-2 border-b-2 border-blue-500"></div>
              </div>
            ) : result ? (
              <div className="space-y-6">
                <MoleculeViewer result={result} />
                <DescriptorsTable descriptors={result.descriptors} />
                <PharmacologicalProperties properties={result.pharmacological_properties} />
              </div>
            ) : (
              <div className="bg-white rounded-lg shadow-md p-6 text-center">
                <p className="text-gray-500">Enter a SMILES string and click "Analyze Molecule" to see results</p>
              </div>
            )}
          </div>
        </div>
      </main>
      
      {/* Footer */}
      <footer className="bg-gray-800 text-white py-6">
        <div className="container mx-auto px-4 text-center">
          <p>Molecular Analysis Tool &copy; {new Date().getFullYear()}</p>
        </div>
      </footer>
    </div>
  );
}

export default Main;