import React, { useState,useEffect } from 'react';
import { useLocation } from "react-router-dom";
import Loading from "../Pages_Components/Model1_Protein/Loading"

function Final() {
  const [substructureInput, setSubstructureInput] = useState('');
  const [drugData, setDrugData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const location = useLocation();
  const data = location.state;

  
  // Sample data in dictionary format
  const sampleData = {
    "Importance": data.name,
    "Bit substructure": data.imp
  };

  useEffect(() => {
    if (Object.keys(sampleData).length > 0) {
      handleGenerateDrug();
    }
  }, []); // Trigger when sampleData changes

  const handleGenerateDrug = async () => {
    setLoading(true);
    setError('');
    try {
      // Convert the dictionary to the desired POST format
      const postData = {
        substructure_input: Object.entries(sampleData).map(([key, value]) => 
          value.map((val, index) => `${key}: ${val}`).join('\n')
        ).join('\n')
      };

      const response = await fetch('http://localhost:5005/generate-drug', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(postData),
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
    // Convert the dictionary to the desired format for display
    const formattedSampleData = Object.entries(sampleData).map(([key, value]) => 
      value.map((val, index) => `${key}: ${val}`).join('\n')
    ).join('\n');
    setSubstructureInput(formattedSampleData);
  };

  return (
    <div className='bg-black'>
      {loading?<div className='flex justify-center items-center h-screen'>
      <Loading/>
      </div>
      :
      <>
{drugData && (
          <div className="bg-black shadow-md rounded-lg p-6 h-fit text-white w-[80vw] mx-auto pt-[50px] pb-[50px]">
            <h2 className="text-3xl font-semibold mb-9 text-center font-poppins">Generated Drug Compound</h2>
            
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <div>
                <img
                  src={`data:image/png;base64,${drugData.image}`}
                  alt="Molecular Structure"
                  className="w-full rounded-[20px]"
                />
              </div>
              
              <div className=' p-2 rounded-[20px] flex justify-center items-center flex-col'>
                {/* <h3 className="text-lg font-semibold mb-4 font-poppins text-center">Drug Information</h3> */}
                {/* <p className='mb-4'><strong>Name:</strong> {drugData.drug_name}</p> */}
                <p className='mb-2 text-lg'><strong>SMILES Notation:</strong></p>
                <pre className="bg-black border border-white p-4 rounded text-lg">{drugData.smiles}</pre>            
              </div>
            </div>

            <div className='border border-white rounded-[20px] mt-8'>
              <h3 className="text-xl font-semibold mt-4 font-poppins text-center">Molecular Properties</h3>
                <div className='grid grid-cols-3 gap-4 p-6'>
                  {Object.entries(drugData.properties).map(([key, value]) => (
                    <div className='border border-white p-4 rounded-[10px]'><p key={key}><strong>{key}:</strong> {value}</p></div>
                  ))}
                </div>
            </div>
          </div>
        )}
      </>}
      

        {error && (
          <div className="mt-6 bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded">
            {error}
          </div>
        )}

        
      </div>
  );
}

export default Final;