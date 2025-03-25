import React, { useState, useEffect } from 'react';
import axios from 'axios';
import MoleculeViewer from './Molecular';

function DrugVista() {
  const [targetName, setTargetName] = useState('');
  const [targets, setTargets] = useState([]);
  const [selectedTargetIndex, setSelectedTargetIndex] = useState(0);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [results, setResults] = useState(null);
  const [step, setStep] = useState('search'); // search, select, results
  const [viewSmiles, setViewSmiles] = useState(null);
  
  const handleSearch = async () => {
    if (!targetName) return;
    
    setLoading(true);
    setError(null);
    
    try {
      const response = await axios.post('http://localhost:5000/api/search_target', {
        targetname: targetName
      });
      
      setTargets(response.data.targets);
      if (response.data.targets.length > 0) {
        setStep('select');
      } else {
        setError('No targets found. Please try a different search term.');
      }
    } catch (err) {
      setError('Error searching for targets: ' + (err.response?.data?.error || err.message));
    } finally {
      setLoading(false);
    }
  };
  
  const handleProcessTarget = async () => {
    setLoading(true);
    setError(null);
    
    try {
      const response = await axios.post('http://localhost:5000/api/process_target', {
        targetname: targetName,
        target_index: selectedTargetIndex
      });
      
      setResults(response.data);
      setStep('results');
    } catch (err) {
      setError('Error processing target: ' + (err.response?.data?.error || err.message));
    } finally {
      setLoading(false);
    }
  };
  
  const handleDownloadData = () => {
    window.location.href = 'http://localhost:5000/api/download_data';
  };
  
  return (
    <div className="min-h-screen bg-gray-50">
      <header className="bg-blue-600 shadow-lg">
        <div className="max-w-7xl mx-auto px-4 py-6 sm:px-6 lg:px-8">
          <h1 className="text-3xl font-bold text-white">Bioactivity Analysis Tool</h1>
        </div>
      </header>
      
      <main className="max-w-7xl mx-auto px-4 py-6 sm:px-6 lg:px-8">
        <div className="bg-white shadow-md rounded-lg p-6 mb-6">
          {step === 'search' && (
            <div>
              <h2 className="text-xl font-semibold mb-4">Search for Target</h2>
              <div className="flex flex-col md:flex-row gap-4">
                <input
                  type="text"
                  value={targetName}
                  onChange={(e) => setTargetName(e.target.value)}
                  placeholder="Enter target name (e.g., 'Acetylcholinesterase', 'EGFR')"
                  className="flex-1 p-2 border border-gray-300 rounded-md focus:ring-blue-500 focus:border-blue-500"
                />
                <button
                  onClick={handleSearch}
                  disabled={loading || !targetName}
                  className="px-4 py-2 bg-blue-500 text-white rounded-md hover:bg-blue-600 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 disabled:opacity-50"
                >
                  {loading ? 'Searching...' : 'Search'}
                </button>
              </div>
              {error && <p className="mt-2 text-red-600">{error}</p>}
            </div>
          )}
          
          {step === 'select' && (
            <div>
              <div className="flex justify-between items-center mb-4">
                <h2 className="text-xl font-semibold">Select Target</h2>
                <button
                  onClick={() => setStep('search')}
                  className="text-blue-500 hover:text-blue-700"
                >
                  Back to Search
                </button>
              </div>
              
              <div className="mb-4">
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Found {targets.length} targets matching "{targetName}"
                </label>
                <select
                  value={selectedTargetIndex}
                  onChange={(e) => setSelectedTargetIndex(parseInt(e.target.value))}
                  className="w-full p-2 border border-gray-300 rounded-md focus:ring-blue-500 focus:border-blue-500"
                >
                  {targets.map((target, index) => (
                    <option key={target.target_chembl_id} value={index}>
                      {target.pref_name} ({target.target_chembl_id})
                    </option>
                  ))}
                </select>
              </div>
              
              <div className="mt-4">
                <button
                  onClick={handleProcessTarget}
                  disabled={loading}
                  className="px-4 py-2 bg-blue-500 text-white rounded-md hover:bg-blue-600 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 disabled:opacity-50"
                >
                  {loading ? 'Processing...' : 'Process Target'}
                </button>
                {error && <p className="mt-2 text-red-600">{error}</p>}
              </div>
            </div>
          )}
          
          {step === 'results' && results && (
            <div>
              <div className="flex justify-between items-center mb-6">
                <h2 className="text-xl font-semibold">Results for {results.selected_target}</h2>
                <button
                  onClick={() => setStep('search')}
                  className="text-blue-500 hover:text-blue-700"
                >
                  New Search
                </button>
              </div>
              
              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                <div className="bg-gray-50 rounded-lg p-4">
                  <h3 className="text-lg font-medium mb-3">Data Distribution Plots</h3>
                  {results.plot_data && (
                    <img 
                      src={`data:image/png;base64,${results.plot_data}`} 
                      alt="Bioactivity Data Plots" 
                      className="w-full rounded-md shadow-sm"
                    />
                  )}
                </div>
                
                <div className="bg-gray-50 rounded-lg p-4">
                  <h3 className="text-lg font-medium mb-3">Model Performance</h3>
                  {results.scatter_plot && (
                    <div className="mb-4">
                      <img 
                        src={`data:image/png;base64,${results.scatter_plot}`} 
                        alt="Model Predictions" 
                        className="w-full rounded-md shadow-sm"
                      />
                    </div>
                  )}
                  
                  {results.model_results?.model_stats && (
                    <div className="grid grid-cols-2 gap-4 mb-4">
                      <div className="bg-white p-3 rounded-md shadow-sm">
                        <p className="text-sm text-gray-500">Mean Squared Error</p>
                        <p className="text-xl font-semibold">{results.model_results.model_stats.mse.toFixed(4)}</p>
                      </div>
                      <div className="bg-white p-3 rounded-md shadow-sm">
                        <p className="text-sm text-gray-500">R² Score</p>
                        <p className="text-xl font-semibold">{results.model_results.model_stats.r2.toFixed(4)}</p>
                      </div>
                    </div>
                  )}
                </div>
              </div>
              
              <div className="mt-6 grid grid-cols-1 lg:grid-cols-2 gap-6">
                <div className="bg-gray-50 rounded-lg p-4">
                  <h3 className="text-lg font-medium mb-3">Top Important Features</h3>
                  <div className="overflow-x-auto">
                    <table className="min-w-full divide-y divide-gray-200">
                      <thead className="bg-gray-100">
                        <tr>
                          <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Feature</th>
                          <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Importance</th>
                          <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Bit Substructure</th>
                        </tr>
                      </thead>
                      <tbody className="bg-white divide-y divide-gray-200">
                        {results.model_results?.top_features?.map((feature, index) => (
                          <tr key={index} className={index % 2 === 0 ? 'bg-white' : 'bg-gray-50'}>
                            <td className="px-4 py-2 text-sm text-gray-900">{feature.Feature}</td>
                            <td className="px-4 py-2 text-sm text-gray-900">{feature.Importance.toFixed(4)}</td>
                            <td className="px-4 py-2 text-sm text-gray-900">{feature['Bit Substructure'] || 'N/A'}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
                
                <div className="bg-gray-50 rounded-lg p-4">
                    <h3 className="text-lg font-medium mb-3">Generated Molecules</h3>
                    <div className="grid grid-cols-1 gap-4">
                        {results.model_results?.generated_molecules?.map((molecule, index) => (
                        <div key={index} className="bg-white p-3 rounded-md shadow-sm">
                            <p className="font-medium">{molecule.description}</p>
                            <p className="text-sm font-mono bg-gray-100 p-2 mt-1 rounded whitespace-pre-wrap break-all">
                            {molecule.smiles}
                            </p>
                            
                            {/* Add the MoleculeViewer component here */}
                            <div className="mt-3">
                            <h4 className="text-sm font-medium text-gray-700 mb-2">3D Structure</h4>
                            <MoleculeViewer smiles={molecule.smiles} height="250px" />
                            </div>
                            
                            <div className="mt-2">
                            {molecule.smiles && (
                                <a 
                                href={`https://pubchem.ncbi.nlm.nih.gov/compound/${encodeURIComponent(molecule.smiles)}`}
                                target="_blank"
                                rel="noopener noreferrer"
                                className="text-blue-500 hover:text-blue-700 text-sm"
                                >
                                View in PubChem
                                </a>
                            )}
                            </div>
                        </div>
                        ))}
                    </div>
                </div>
              </div>
              
              <div className="mt-6">
                <h3 className="text-lg font-medium mb-3">Bioactivity Data Sample</h3>
                <div className="overflow-x-auto bg-white rounded-md shadow-sm">
                  <table className="min-w-full divide-y divide-gray-200">
                    <thead className="bg-gray-100">
                      <tr>
                        <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Molecule ID</th>
                        <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">SMILES</th>
                        <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Bioactivity Class</th>
                        <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">pIC50</th>
                      </tr>
                    </thead>
                    <tbody className="bg-white divide-y divide-gray-200">
                    {results.bioactivity_data?.slice(0, 10).map((item, index) => (
                        <tr key={index} className={index % 2 === 0 ? 'bg-white' : 'bg-gray-50'}>
                        <td className="px-4 py-2 text-sm text-gray-900">{item.molecule_chembl_id}</td>
                        <td className="px-4 py-2 text-sm font-mono text-gray-900 whitespace-pre-wrap break-all max-w-xs">
                            {item.canonical_smiles?.length > 40 
                            ? item.canonical_smiles.substring(0, 40) + '...' 
                            : item.canonical_smiles}
                            <button 
                            className="ml-2 text-xs text-blue-500 hover:text-blue-700"
                            onClick={() => setViewSmiles(item.canonical_smiles)}
                            >
                            View 3D
                            </button>
                        </td>
                        {/* Rest of the row content */}
                        </tr>
                    ))}
                    </tbody>
                  </table>

                  {viewSmiles && (
                    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
                        <div className="bg-white p-6 rounded-lg max-w-2xl w-full">
                        <div className="flex justify-between items-center mb-4">
                            <h3 className="text-lg font-medium">3D Molecule Viewer</h3>
                            <button 
                            onClick={() => setViewSmiles(null)}
                            className="text-gray-500 hover:text-gray-700"
                            >
                            Close
                            </button>
                        </div>
                        <MoleculeViewer smiles={viewSmiles} height="400px" />
                        </div>
                    </div>
                    )}
                </div>
                
                <div className="mt-4">
                  <button
                    onClick={handleDownloadData}
                    className="px-4 py-2 bg-green-500 text-white rounded-md hover:bg-green-600 focus:outline-none focus:ring-2 focus:ring-green-500 focus:ring-offset-2"
                  >
                    Download Complete Dataset
                  </button>
                </div>
              </div>
            </div>
          )}
        </div>
      </main>
      
      <footer className="bg-gray-100 border-t border-gray-200">
        <div className="max-w-7xl mx-auto px-4 py-6 sm:px-6 lg:px-8">
          <p className="text-center text-gray-500 text-sm">
            Bioactivity Analysis Tool — Powered by RDKit, scikit-learn and ChemGPT
          </p>
        </div>
      </footer>
    </div>
  );
}

export default DrugVista;