// App.jsx
import React, { useState, useEffect, useRef } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';
import Plot from 'react-plotly.js';

function MolViewer({ pdbData }) {
  const viewerRef = useRef(null);
  
  useEffect(() => {
    if (pdbData && viewerRef.current) {
      // Initialize 3Dmol viewer
      const viewer = window.$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: 'white',
      });
      
      // Add the PDB data
      viewer.addModel(pdbData, "pdb");
      
      // Style the protein
      viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
      
      // Zoom to fit the protein
      viewer.zoomTo();
      
      // Start spinning
      viewer.spin(true);
      
      // Render the scene
      viewer.render();
      
      return () => {
        viewer.clear();
      };
    }
  }, [pdbData]);

  return (
    <div 
      ref={viewerRef} 
      style={{ height: '400px', width: '100%', position: 'relative' }}
      className="border rounded"
    />
  );
}

function AminoAcidComposition({ sequence }) {
  const aminoAcids = [...new Set(sequence)].sort();
  const counts = aminoAcids.map(aa => ({
    aa,
    count: (sequence.match(new RegExp(aa, 'g')) || []).length,
    percentage: ((sequence.match(new RegExp(aa, 'g')) || []).length / sequence.length) * 100
  }));

  return (
    <div>
      {/* Pie Chart */}
      <Plot
        data={[{
          values: counts.map(c => c.count),
          labels: counts.map(c => c.aa),
          type: 'pie',
          hole: 0.4,
          marker: {
            colors: counts.map(() => 
              `hsl(${Math.random() * 360}, 70%, 50%)`
            )
          }
        }]}
        layout={{
          height: 400,
          title: 'Amino Acid Distribution',
          showlegend: true,
          legend: {
            orientation: 'h',
            y: -0.2
          }
        }}
        config={{ responsive: true }}
      />

      {/* Detailed Table */}
      <div className="mt-4 overflow-x-auto">
        <table className="min-w-full divide-y divide-gray-200">
          <thead>
            <tr>
              <th className="px-6 py-3 bg-gray-50 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Amino Acid
              </th>
              <th className="px-6 py-3 bg-gray-50 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Count
              </th>
              <th className="px-6 py-3 bg-gray-50 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Percentage
              </th>
            </tr>
          </thead>
          <tbody className="bg-white divide-y divide-gray-200">
            {counts.map(({ aa, count, percentage }) => (
              <tr key={aa}>
                <td className="px-6 py-4 whitespace-nowrap">{aa}</td>
                <td className="px-6 py-4 whitespace-nowrap">{count}</td>
                <td className="px-6 py-4 whitespace-nowrap">{percentage.toFixed(2)}%</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function ProteinModel() {
  const [sequence, setSequence] = useState('');
  const [defaultSequences, setDefaultSequences] = useState({});
  const [selectedSequence, setSelectedSequence] = useState('Custom');
  const [analysisOptions, setAnalysisOptions] = useState([
    'Structure Prediction',
    'Sequence Composition',
    'Hydrophobicity'
  ]);
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  useEffect(() => {
    fetch('http://localhost:5000/api/sequences')
      .then(res => res.json())
      .then(data => setDefaultSequences(data))
      .catch(err => setError('Failed to load default sequences'));
  }, []);

  const handleAnalyze = async () => {
    setLoading(true);
    setError(null);
    try {
      const response = await fetch('http://localhost:5000/api/analyze', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          sequence,
          analysis_options: analysisOptions
        }),
      });
      const data = await response.json();
      if (!response.ok) throw new Error(data.error);
      setResults(data);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-black py-8">
      <div className="max-w-7xl mx-auto px-4">
      <h1 className="text-3xl font-bold mb-8 text-white font-monument">🧬 BioFold: Protein Analysis Platform</h1>
        
        {/* Input Section */}
        <div className="bg-white rounded-lg shadow p-6 mb-8">
          <select 
            className="w-full p-2 mb-4 border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-blue-500"
            value={selectedSequence}
            onChange={(e) => {
              setSelectedSequence(e.target.value);
              if (e.target.value !== 'Custom') {
                setSequence(defaultSequences[e.target.value]);
              }
            }}
          >
            <option value="Custom">Custom Sequence</option>
            {Object.keys(defaultSequences).map(key => (
              <option key={key} value={key}>{key}</option>
            ))}
          </select>

          <textarea
            className="w-full p-2 mb-4 border border-gray-300 rounded h-32 focus:outline-none focus:ring-2 focus:ring-blue-500"
            value={sequence}
            onChange={(e) => setSequence(e.target.value)}
            placeholder="Enter protein sequence..."
          />

          <div className="mb-4">
            <h3 className="font-semibold mb-2">Analysis Options:</h3>
            <div className="space-y-2">
              {['Structure Prediction', 'Sequence Composition', 'Hydrophobicity'].map(option => (
                <label key={option} className="flex items-center space-x-2">
                  <input
                    type="checkbox"
                    className="form-checkbox"
                    checked={analysisOptions.includes(option)}
                    onChange={(e) => {
                      if (e.target.checked) {
                        setAnalysisOptions([...analysisOptions, option]);
                      } else {
                        setAnalysisOptions(analysisOptions.filter(o => o !== option));
                      }
                    }}
                  />
                  <span>{option}</span>
                </label>
              ))}
            </div>
          </div>

          <button
            className={`w-full py-2 px-4 rounded font-semibold text-white ${
              loading || !sequence 
                ? 'bg-gray-400 cursor-not-allowed' 
                : 'bg-blue-500 hover:bg-blue-600'
            }`}
            onClick={handleAnalyze}
            disabled={loading || !sequence}
          >
            {loading ? 'Analyzing...' : 'Analyze Protein'}
          </button>
        </div>

        {/* Error Display */}
        {error && (
          <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-4">
            {error}
          </div>
        )}


        {/* Results Section */}
        {results && (
          <div className="space-y-8">
            {/* Molecular Structure Viewer */}
            {results.structure && !results.structure.error && (
              <div className="bg-white rounded-lg shadow p-6">
                <h2 className="text-xl font-bold mb-4">Protein Structure</h2>
                <MolViewer pdbData={results.structure.pdb_string} />
                <div className="mt-4">
                  <p className="font-medium">Confidence Score (pLDDT):</p>
                  <p className="text-2xl">{results.structure.b_value.toFixed(2)}/100</p>
                </div>
                <button
                  className="mt-4 bg-green-500 hover:bg-green-600 text-white py-2 px-4 rounded"
                  onClick={() => {
                    const blob = new Blob([results.structure.pdb_string], { type: 'text/plain' });
                    const url = window.URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'predicted_structure.pdb';
                    a.click();
                  }}
                >
                  Download PDB File
                </button>
              </div>
            )}

            {/* Amino Acid Composition */}
            <div className="bg-white rounded-lg shadow p-6">
              <h2 className="text-xl font-bold mb-4">Amino Acid Composition</h2>
              <AminoAcidComposition sequence={sequence} />
            </div>

                        {/* Molecular Properties */}
                        <div className="bg-white rounded-lg shadow p-6">
              <h2 className="text-xl font-bold mb-4">Molecular Properties</h2>
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
                {Object.entries(results.molecular_properties).map(([key, value]) => (
                  <div key={key} className="bg-gray-50 p-4 rounded">
                    <h3 className="font-semibold capitalize">{key.replace(/_/g, ' ')}</h3>
                    <p className="text-2xl">{typeof value === 'number' ? value.toFixed(2) : value}</p>
                  </div>
                ))}
              </div>
            </div>

            {/* Hydrophobicity Plot */}
            {results.hydrophobicity && (
              <div className="bg-white rounded-lg shadow p-6">
                <h2 className="text-xl font-bold mb-4">Hydrophobicity Profile</h2>
                <div className="h-64">
                  <ResponsiveContainer width="100%" height="100%">
                    <LineChart
                      data={results.hydrophobicity.map((value, index) => ({
                        position: index,
                        value
                      }))}
                    >
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis dataKey="position" />
                      <YAxis />
                      <Tooltip />
                      <Line type="monotone" dataKey="value" stroke="#3B82F6" />
                    </LineChart>
                  </ResponsiveContainer>
                </div>
              </div>
            )}

            {/* Advanced Analysis */}
            {results.advanced_analysis && (
              <div className="bg-white rounded-lg shadow p-6">
                <h2 className="text-xl font-bold mb-4">Advanced Analysis</h2>
                
                {/* Motifs */}
                <div className="mb-6">
                  <h3 className="font-semibold mb-2">Detected Motifs:</h3>
                  {Object.entries(results.advanced_analysis.motifs).map(([motif, positions]) => (
                    <div key={motif} className="mb-4 bg-gray-50 p-3 rounded">
                      <p className="font-medium">{motif}:</p>
                      <p className="text-gray-600">
                        {positions.length > 0 ? 
                          `Found at positions: ${positions.join(', ')}` : 
                          'None found'}
                      </p>
                    </div>
                  ))}
                </div>

                {/* Secondary Structure */}
                <div>
                  <h3 className="font-semibold mb-2">Secondary Structure Composition:</h3>
                  {Object.entries(results.advanced_analysis.secondary_structure).map(([structure, value]) => (
                    <div key={structure} className="mb-3">
                      <div className="flex items-center justify-between mb-1">
                        <span className="font-medium">{structure}</span>
                        <span>{(value * 100).toFixed(1)}%</span>
                      </div>
                      <div className="w-full bg-gray-200 rounded-full h-2">
                        <div
                          className="bg-blue-500 rounded-full h-2"
                          style={{ width: `${value * 100}%` }}
                        />
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Structure Prediction Results */}
            {results.structure && !results.structure.error && (
              <div className="bg-white rounded-lg shadow p-6">
                <h2 className="text-xl font-bold mb-4">Structure Prediction</h2>
                <div className="mb-4">
                  <p className="font-medium">Confidence Score (pLDDT):</p>
                  <p className="text-2xl">{results.structure.b_value.toFixed(2)}/100</p>
                </div>
                <button
                  className="bg-green-500 hover:bg-green-600 text-white py-2 px-4 rounded"
                  onClick={() => {
                    const blob = new Blob([results.structure.pdb_string], { type: 'text/plain' });
                    const url = window.URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = 'predicted_structure.pdb';
                    a.click();
                  }}
                >
                  Download PDB File
                </button>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

export default ProteinModel;