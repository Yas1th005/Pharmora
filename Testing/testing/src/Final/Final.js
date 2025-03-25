import React, { useState } from 'react';

function Final() {
  const [substructureInput, setSubstructureInput] = useState('');
  const [drugData, setDrugData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  // Sample data
  const sampleData = `
Importance    Bit substructure
0.0416       >= 5 unsaturated non-aromatic heteroatom-containing ring size 6
0.0303       >= 5 unsaturated non-aromatic heteroatom-containing ring size 5
0.0301       >= 2 saturated or aromatic heteroatom-containing ring size 4
0.0190       >= 1 Bi
0.0184       >= 3 saturated or aromatic nitrogen-containing ring size 6
0.0152       >= 1 Ni
  `;

  const handleGenerateDrug = async () => {
    setLoading(true);
    setError('');
    try {
      const response = await fetch('http://localhost:5000/generate-drug', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ substructure_input: substructureInput }),
      });
      const data = await response.json();
      if (data.success) {
        setDrugData(data);
      } else {
        setError(data.error || 'Failed to generate drug compound');
      }
    } catch (err) {
      setError('An error occurred while generating the drug compound');
    } finally {
      setLoading(false);
      
    }
  };

  const handleUseSampleData = () => {
    setSubstructureInput(sampleData.trim());
  };
  console.log(drugData)
  return (
    <div className="min-h-screen bg-gray-100 py-8">
      <div className="max-w-4xl mx-auto px-4">
        <h1 className="text-3xl font-bold text-center mb-8">AI Drug Designer</h1>
        
        <div className="bg-white shadow-md rounded-lg p-6">
          <h2 className="text-xl font-semibold mb-4">Enter Substructure Data</h2>
          <textarea
            className="w-full p-2 border border-gray-300 rounded mb-4"
            rows="6"
            placeholder="Enter substructure data..."
            value={substructureInput}
            onChange={(e) => setSubstructureInput(e.target.value)}
          ></textarea>
          
          <div className="flex gap-2">
            <button
              className="bg-blue-500 text-white px-6 py-2 rounded hover:bg-blue-600"
              onClick={handleGenerateDrug}
              disabled={loading}
            >
              {loading ? 'Generating...' : 'Generate Drug Compound'}
            </button>
            <button
              className="bg-gray-500 text-white px-6 py-2 rounded hover:bg-gray-600"
              onClick={handleUseSampleData}
            >
              Use Sample Data
            </button>
          </div>
        </div>

        {error && (
          <div className="mt-6 bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded">
            {error}
          </div>
        )}

        {drugData && (
          <div className="mt-6 bg-white shadow-md rounded-lg p-6">
            <h2 className="text-xl font-semibold mb-4">Generated Drug Compound</h2>
            
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <div>
                <h3 className="text-lg font-semibold mb-2">Molecular Structure</h3>
                <img
                  src={`data:image/png;base64,${drugData.image}`}
                  alt="Molecular Structure"
                  className="w-full rounded"
                />
              </div>
              
              <div>
                <h3 className="text-lg font-semibold mb-2">Drug Information</h3>
                <p><strong>Name:</strong> {drugData.drug_name}</p>
                <p><strong>SMILES Notation:</strong></p>
                <pre className="bg-gray-100 p-2 rounded">{drugData.smiles}</pre>
                
                <h3 className="text-lg font-semibold mt-4 mb-2">Molecular Properties</h3>
                <ul>
                  {Object.entries(drugData.properties).map(([key, value]) => (
                    <li key={key}><strong>{key}:</strong> {value}</li>
                  ))}
                </ul>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default Final;