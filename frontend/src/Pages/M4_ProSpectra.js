// src/App.js
import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { ToastContainer, toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';
import MoleculeViewer from '../Pages_Components/Model4_ProSpectra/Molecule';
import DescriptorsTable from '../Pages_Components/Model4_ProSpectra/Description';
import PharmacologicalProperties from '../Pages_Components/Model4_ProSpectra/Properties';
import ExampleMolecules from '../Pages_Components/Model4_ProSpectra/Example';
import Loading from "../Pages_Components/Model1_Protein/Loading"

function ProSpectra() {
  const [smiles, setSmiles] = useState('');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [exampleMolecules, setExampleMolecules] = useState({});

  useEffect(() => {
    // Fetch example molecules on component mount
    axios.get('http://127.0.0.1:5003/api/examples')
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
      const response = await axios.post('http://127.0.0.1:5003/api/analyze', { smiles });
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
    <div className="min-h-screen bg-black">
      <ToastContainer position="top-right" autoClose={3000} />
      
      {/* Header */}
      <header className="bg-black text-white shadow-lg py-4">
        <div className="container mx-auto px-4">
          <h1 className="text-3xl font-bold text-center font-poppins mt-4">Molecular Analysis Tool</h1>
          {/* <p className="text-white text-center font-roboto">Analyze molecules using SMILES notation</p> */}
        </div>
      </header>
      
      {/* Main Content */}
      <main className="container mx-auto px-4 py-8">
        <div className="">
          {/* Left Column - Input Form */}
          <div className="grid grid-cols-2 gap-4">
            <div className="bg-black text-white rounded-lg shadow-md mb-6 border border-white h-[300px] w-[100%] col-span-1 flex flex-col justify-center items-center">
              <h2 className="text-xl font-semibold mb-4 text-center font-poppins">Enter SMILES String</h2>
              <form onSubmit={handleSubmit}>
                <div className="mb-4">
                  {/* <label htmlFor="smiles" className="block text-white mb-2 text-center font-roboto">SMILES String</label> */}
                  <textarea
                    id="smiles"
                    className="w-[650px] px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500 text-white bg-black"
                    rows="3"
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O for Aspirin"
                  />
                </div>
                <button
                  type="submit"
                  disabled={loading}
                  className="w-full py-2 font-semibold rounded-[10px] bg-violet-400 text-gray-900 font-poppins"
                >
                  Analyze Molecule
                </button>
              </form>
            </div>
            
            <ExampleMolecules 
              examples={exampleMolecules} 
              onSelect={handleExampleSelect} 
            />
          </div>
          
          {/* Right Column - Results */}
          <div className="">
            {loading ? (
              <div className="flex justify-center items-center h-screen">
                <Loading/>
              </div>
            ) : result ? (
              <div className="space-y-6">
                <MoleculeViewer result={result} />
                <DescriptorsTable descriptors={result.descriptors} />
                <PharmacologicalProperties properties={result.pharmacological_properties} />
              </div>
            ) : (
              <div className="bg-black rounded-lg shadow-md p-6 text-center">
                <p className="text-gray-500">Enter a SMILES string and click "Analyze Molecule" to see results</p>
              </div>
            )}
          </div>
        </div>
      </main>
  
    </div>
  );
}

export default ProSpectra;