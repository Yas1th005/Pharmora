import React, { useState, useEffect } from 'react';
import axios from 'axios';
import MoleculeViewer from '../Pages_Components/Model6_DrugVista/MolecularViewer';
import { Link } from "react-router-dom";
import Loading from "../Pages_Components/Model1_Protein/Loading";
import one from "../imgaes/Icons/1.png" 
import two from "../imgaes/Icons/2.png" 
import three from "../imgaes/Icons/3.png" 

function DrugVista() {
  const [targetName, setTargetName] = useState('');
  const [targets, setTargets] = useState([]);
  const [selectedTargetIndex, setSelectedTargetIndex] = useState(0);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [results, setResults] = useState(null);
  const [step, setStep] = useState('search'); // search, select, results
  const [viewSmiles, setViewSmiles] = useState(null);
  const Array1=[];
  const Array2=[];

  const send={imp:Array1,name:Array2};

  console.log(results)

  const serviceList = [
    {
      icon: one,
      title: "Bioactivity Profiling and Classification",
      description:
        "Analyzes ChEMBL bioactivity data to classify the compounds as active or inactive against targets while calculating Lipinski descriptors and pIC50 values for structure-activity visualization.",
    },
    {
      icon: two,
      title: "AI-Driven Molecular Feature Analysis",
      description:
        "Employs Random Forest models to identify critical molecular descriptors that contribute to target binding and also revealing key pharmacophore elements through variance threshold filtering and statistical validation.",
    },
    {
      icon: three,
      title: "Generative Chemistry for Novel Drug Discovery",
      description:
        "Integrates language models to transform high-ranking molecular fingerprints into novel SMILES strings, bridging the gap between feature identification and practical drug candidate generation.",
    },
  ];

  const ServiceItem = ({ service }) => (
    <div className="bg-black border border-white shadow-xl rounded-xl h-full">
      <div className="p-6 md:p-12">
        <div className="w-[75px] h-[75px] rounded-full text-[26px] shadow-xl flex justify-center items-center mb-6 text-white">
          <img src={service.icon}></img>
        </div>
        <h4 className="text-2xl mb-6 font-bold text-white">{service.title}</h4>
        <p className="opacity-70 leading-[1.8] text-white text-justify">{service.description}</p>
      </div>
    </div>
  );
  
  const handleSearch = async () => {
    if (!targetName) return;
    
    setLoading(true);
    setError(null);
    
    try {
      const response = await axios.post('http://localhost:5004/api/search_target', {
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
      const response = await axios.post('http://localhost:5004/api/process_target', {
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
    window.location.href = 'http://localhost:5004/api/download_data';
  };
  
  console.log(results)

  return (
    <div className="min-h-screen bg-black">
      {/* <header className="bg-black shadow-lg">
        <div className="max-w-7xl mx-auto px-4 py-6 sm:px-6 lg:px-8">
          <h1 className="text-3xl font-bold text-white text-center font-poppins">Bioactivity Analysis Tool</h1>
        </div>
      </header> */}
      
      <main className="w-[100%] mx-auto p-2 overflow-x-hidden">
        <div className="bg-black">
          {step === 'search' && (
            // <div>
            //   <h2 className="text-xl font-semibold mb-4">Search for Target</h2>
            //   <div className="flex flex-col md:flex-row gap-4">
            //     <input
            //       type="text"
            //       value={targetName}
            //       onChange={(e) => setTargetName(e.target.value)}
            //       placeholder="Enter target name (e.g., 'Acetylcholinesterase', 'EGFR')"
            //       className="flex-1 p-2 border border-gray-300 rounded-md focus:ring-blue-500 focus:border-blue-500"
            //     />
            //     <button
            //       onClick={handleSearch}
            //       disabled={loading || !targetName}
            //       className="px-4 py-2 bg-blue-500 text-white rounded-md hover:bg-blue-600 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 disabled:opacity-50"
            //     >
            //       {loading ? 'Searching...' : 'Search'}
            //     </button>
            //   </div>
            //   {error && <p className="mt-2 text-red-600">{error}</p>}
            // </div>

            <section className="light bg-black dark:text-white relative z-[1] rounded-[20px]">
              <div className="absolute bottom-0 left-0 right-0 h-1/2 w-full-z-[1]" />
              <div className="container px-4 mx-auto">
                <div className="">
                  <div className="w-[80vw] mx-auto">
                    <h2 className="text-3xl font-bold text-white font-poppins mb-5 text-center mt-7">
                    Drug Discovery Toolkit
                    </h2>
                    <div>
                      <h2 className="text-xl font-semibold mb-3 text-white font-poppins">Enter Target Name</h2>
                      <div className="flex flex-col md:flex-row gap-4">
                        <input
                          type="text"
                          value={targetName}
                          onChange={(e) => setTargetName(e.target.value)}
                          placeholder="e.g. Acetylcholinesterase, EGFR"
                          className="flex-1 p-4 border border-white bg-black rounded-md focus:ring-blue-500 focus:border-blue-500 w-full text-white"
                        />
                        <button
                          onClick={handleSearch}
                          disabled={loading || !targetName}
                          className="w-[200px] py-2 font-semibold rounded-[10px] bg-violet-400 text-gray-900 font-poppins hover:cursor-pointer"
                        >
                          Search
                        </button>
                      </div>
                      {error && <p className="mt-2 text-red-600">{error}</p>}
                    </div>
                  </div>
                </div>
                <div className="grid grid-cols-6 gap-6 w-[80vw] mx-auto mt-12">
                  {serviceList.map((service, i) => (
                    <div className="col-span-6 md:col-span-3 lg:col-span-2" key={i}>
                      <ServiceItem service={service} />
                    </div>
                  ))}
                </div>
              </div>
            </section>
          )}

          {!loading?<>
            {step === 'select' && (
            <div>
                <div className="text-white w-[80vw] mx-auto mb-4">
                <h2 className="text-2xl font-semibold font-poppins text-center mt-6">Select Target</h2>
                </div>
                {/* <button
                  onClick={() => setStep('search')}
                  className="w-[200px] py-2 font-semibold rounded-[10px] bg-violet-400 text-gray-900 font-poppins hover:cursor-pointer"
                >
                  Back to Search
                </button> */}
                
                <div className="mb-4">
                <p className="text-base font-medium text-white font-poppins mb-2 text-center mb-4 mt-[-10px]">
                    Found {targets.length} targets matching "{targetName}"
                </p>
                
                <div className="overflow-auto h-fit border border-white rounded-md w-[90vw] mx-auto overflow-x-hidden">
                    <table className="min-w-full divide-y divide-gray-200">
                    <thead className="bg-violet-400 text-black font-poppins">
                        <tr>
                        <th scope="col" className="px-6 py-3 text-left text-base font-medium font-poppins text-black uppercase tracking-wider">
                            Organism
                        </th>
                        <th scope="col" className="px-6 py-3 text-left text-base font-medium font-poppins text-black uppercase tracking-wider">
                            Pref_Name
                        </th>
                        <th scope="col" className="px-6 py-3 text-left text-base font-medium font-poppins text-black uppercase tracking-wider">
                            Score
                        </th>
                        <th scope="col" className="px-6 py-3 text-left text-base font-medium font-poppins text-black uppercase tracking-wider">
                            ChEMBL ID
                        </th>
                        <th scope="col" className="px-6 py-3 text-left text-base font-medium font-poppins text-black uppercase tracking-wider">
                            Type
                        </th>
                        <th scope="col" className="px-6 py-3 text-left text-base font-medium font-poppins text-black uppercase tracking-wider">
                            Tax ID
                        </th>
                        </tr>
                    </thead>
                    <tbody className="bg-black divide-y divide-white">
                        {targets.map((target, index) => (
                        <tr 
                            key={target.target_chembl_id}
                            onClick={() => setSelectedTargetIndex(index)}
                            onDoubleClick={() => {
                            setSelectedTargetIndex(index);
                            handleProcessTarget();
                            }}
                            className={`cursor-pointer hover:bg-white/10 ${selectedTargetIndex === index ? 'bg-white/30' : ''}`}
                        >
                            <td className="px-6 py-4 whitespace-nowrap text-sm font-medium text-white">
                            {target.organism}
                            </td>
                            <td className="px-6 py-4 whitespace-nowrap text-sm text-white">
                            {target.pref_name}
                            </td>
                            <td className="px-6 py-4 whitespace-nowrap text-sm text-white">
                            {target.score}
                            </td>
                            <td className="px-6 py-4 whitespace-nowrap text-sm text-white">
                            {target.target_chembl_id}
                            </td>
                            <td className="px-6 py-4 whitespace-nowrap text-sm text-white">
                            {target.target_type}
                            </td>
                            <td className="px-6 py-4 whitespace-nowrap text-sm text-white">
                            {target.tax_id}
                            </td>
                        </tr>
                        ))}
                    </tbody>
                    </table>
                </div>
                
                {/* {loading && (
                    <div className="mt-4 flex items-center text-sm text-gray-600">
                    <svg className="animate-spin h-5 w-5 mr-2 text-blue-500" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                    </svg>
                    Processing...
                    </div>
                )} */}
                
                {error && <p className="mt-2 text-red-600">{error}</p>}
                </div>
            </div>
            )}
          </>:<div className='flex justify-center items-center h-screen bg-black overflow-x-hidden'>
          <Loading/>
          </div>}
          
          
          {step === 'results' && results && (
            <div className='w-[80vw] flex flex-col justify-center items-center mx-auto'>
              <div className="mt-2">
                <h2 className="text-2xl font-semibold text-white text-center font-poppins">Detailed Analysis for {results.selected_target}</h2>
              </div>


              <div className="mt-4">
                <h3 className="text-lg font-medium mb-3 text-white">Bioactivity Data Sample</h3>
                <div className="overflow-x-auto bg-white rounded-md shadow-sm border border-white">
                  <table className="min-w-full divide-y divide-gray-200 w-[80vw]">
                    <thead className="bg-violet-400">
                      <tr>
                        <th className="px-4 py-2 text-left text-base font-medium text-black font-poppins uppercase tracking-wider">Molecule ID</th>
                        <th className="px-4 py-2 text-left text-base font-medium text-black font-poppins uppercase tracking-wider">SMILES</th>
                        <th className="px-4 py-2 text-left text-base font-medium text-black font-poppins uppercase tracking-wider">Bioactivity Class</th>
                        <th className="px-4 py-2 text-left text-base font-medium text-black font-poppins uppercase tracking-wider">pIC50</th>
                      </tr>
                    </thead>
                    <tbody className="bg-black divide-y divide-gray-200">
                    {results.bioactivity_data?.slice(0, 10).map((item, index) => (
                        <tr key={index} className={index % 2 === 0 ? 'bg-black' : 'bg-white/10'}>
                        <td className="px-4 py-2 text-sm text-white">{item.molecule_chembl_id}</td>
                        <td className="px-4 py-2 text-sm font-mono text-white whitespace-pre-wrap break-all max-w-xs">
                            {item.canonical_smiles?.length > 40 
                            ? item.canonical_smiles.substring(0, 40) + '...' 
                            : item.canonical_smiles}
                        </td>
                        <td className="px-4 py-2 text-sm text-white capitalize">{item.bioactivity_class}</td>
                        <td className="px-4 py-2 text-sm text-white">{-Math.log10(item.standard_value).toFixed(3)}</td>

                        {/* Rest of the row content */}
                        </tr>
                    ))}
                    </tbody>
                  </table>

                </div>
          
              </div>
              
              <div className="mt-[30px]">
                <div className="bg-black border border-white rounded-lg px-6 py-3 text-white h-fit mb-5">
                  <h3 className="text-lg font-medium mb-3 text-center font-poppins">Drug Distribution Plots</h3>
                  {results.plot_data && (
                    <img 
                      src={`data:image/png;base64,${results.plot_data}`} 
                      alt="Bioactivity Data Plots" 
                      className="w-full rounded-md shadow-sm"
                    />
                  )}
                </div>

                <DatasetTable dataset={results.dataset} />
                
                <div className='flex justify-center items-center w-[80vw] gap-5 mt-7'>
                  
                  <div className="">
                <div className="bg-black text-white border border-white rounded-lg p-4 h-[700px]">
                  <h3 className="text-lg font-medium mb-3 text-center font-poppins">Most Important Features</h3>
                  <div className="overflow-x-auto border border-white rounded-[10px]">
                    <table className="min-w-full divide-y divide-gray-200">
                      <thead className="bg-violet-400">
                        <tr className='border border-white'>
                          <th className="px-4 py-2 text-left text-xs font-medium text-black font-poppins uppercase tracking-wider">Feature</th>
                          <th className="px-4 py-2 text-left text-xs font-medium text-black font-poppins uppercase tracking-wider">Importance</th>
                          <th className="px-4 py-2 text-left text-xs font-medium text-black font-poppins uppercase tracking-wider">Bit Substructure</th>
                        </tr>
                      </thead>
                      <tbody className="bg-black border border-white divide-y divide-white">
                        {results.model_results?.top_features?.map((feature, index) => {
                          if (index <= 4) {
                            Array1.push(feature['Bit Substructure'] || 'N/A');
                            Array2.push(Number(feature.Importance.toFixed(4)));
                          }

                          return (
                            <tr key={index} className={index % 2 === 0 ? 'bg-black' : 'bg-white/10'}>
                              <td className="px-4 py-2 text-sm text-white">{feature.Feature}</td>
                              <td className="px-4 py-2 text-sm text-white">{feature.Importance.toFixed(4)}</td>
                              <td className="px-4 py-2 text-sm text-white">{feature['Bit Substructure'] || 'N/A'}</td>
                            </tr>
                          );
                        })}
                      </tbody>
                    </table>
                  </div>
                </div>
              </div>

              <div className="bg-black border border-white rounded-lg p-4 h-[700px] flex flex-col justify-center items-center">
                    <h3 className="text-lg font-medium mb-3 text-white">Model Performance</h3>
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
                      <div className="grid grid-cols-2 gap-4 mb-4 w-full mt-4">
                        <div className="bg-black border border-white p-3 rounded-md shadow-sm">
                          <p className="text-sm text-white text-center font-poppins">Mean Squared Error</p>
                          <p className="text-xl font-semibold text-white text-center">{results.model_results.model_stats.mse.toFixed(4)}</p>
                        </div>
                        <div className="bg-black p-3 rounded-md shadow-sm border border-white">
                          <p className="text-sm text-white text-center font-poppins">R² Score</p>
                          <p className="text-xl font-semibold text-white text-center">{results.model_results.model_stats.r2.toFixed(4)}</p>
                        </div>
                      </div>
                    )}
                  </div>

              </div>
              
              
                
                

              </div>

              

              <div className="mt-6 grid grid-cols-4 gap-4 w-[80vw] mb-[50px]">
                
      
                <div className="col-span-1">
                  <button
                    onClick={handleDownloadData}
                    className="px-4 py-4 bg-transparent border border-white text-white rounded-lg hover:bg-green-600 focus:outline-none focus:ring-2 focus:ring-green-500 focus:ring-offset-2"
                  >
                    Download Complete Report
                  </button>
                </div>
                <div className="col-span-1">
                  <button
                    onClick={handleDownloadData}
                    className="px-4 py-4 bg-transparent border border-white text-white rounded-lg hover:bg-green-600 focus:outline-none focus:ring-2 focus:ring-green-500 focus:ring-offset-2"
                  >
                    Download PubChem Doc
                  </button>
                </div>

                <Link className='col-span-2 py-4 font-semibold rounded-[10px] bg-violet-400 text-gray-900 font-poppins hover:cursor-pointer text-center' to={"/Final"} state={send}>Generate New Drug using Top Features</Link>
              </div>
            </div>
          )}
        </div>
      </main>
    </div>
  );
}

export default DrugVista;



const DatasetTable = ({ dataset }) => {
  if (!dataset || dataset.length === 0) {
    return <p className="text-white">No dataset available.</p>;
  }

  return (
    <div className="mt-6 w-[80vw]">
      <h3 className="text-lg font-medium mb-3 text-white font-poppins">Molecular Fingerprints</h3>
      <div className="overflow-x-auto bg-black border border-white rounded-lg shadow-sm">
        <div className="max-h-[400px] overflow-y-auto"> {/* Scrollable container */}
          <table className="min-w-full divide-y divide-white">
            <thead className="bg-violet-400 sticky top-0 z-10"> {/* Sticky header */}
              <tr>
                <th className="px-4 py-2 text-left text-xs font-medium text-black uppercase tracking-wider">
                  #
                </th>
                {Object.keys(dataset[0]).map((key) => (
                  <th
                    key={key}
                    className="px-4 py-2 text-left text-xs font-medium text-black uppercase tracking-wider"
                  >
                    {key}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody className="bg-black divide-y divide-white">
              {dataset.map((row, rowIndex) => (
                <tr key={rowIndex} className={rowIndex % 2 === 0 ? 'bg-black' : 'bg-white/10'}>
                  <td className="px-4 py-2 text-sm text-white whitespace-nowrap">{rowIndex + 1}</td> {/* Row index */}
                  {Object.values(row).map((value, colIndex) => (
                    <td
                      key={colIndex}
                      className="px-4 py-2 text-sm text-white whitespace-nowrap"
                    >
                      {typeof value === 'string' && value.length > 50
                        ? `${value.substring(0, 50)}...`
                        : value}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
};
