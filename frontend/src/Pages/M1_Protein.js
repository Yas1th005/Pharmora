import React,{useState,useEffect} from 'react'
import Popup from '../Pages_Components/Model1_Protein/Popup';
import Loading from '../Pages_Components/Model1_Protein/Loading';
import MainResult from '../Pages_Components/Model1_Protein/Results/MainResult';
import drugvista from "../imgaes/DrugVista/DrugVista.webp"
import LandModel from '../components/LandModel';

export default function M1_Protein() {
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
    const [hero, setHero] = useState(true);
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
        setHero(false);
      }
    };


  return (
    <div className='flex justify-center items-center h-screen bg-black overflow-x-hidden'>
      {/* HERO SECTION */}
      {hero?
      !loading?
        <section className="p-6 bg-black text-gray-100 ">
          {/* <h1 className="text-3xl font-bold mb-8 text-white font-monument text-center tracking-[3px]">🧬 BioFold: Protein Analysis Platform</h1> */}
          <div className="grid gap-6 mx-auto text-center grid-cols-5 w-[90vw]">
            <div className="w-full px-6 py-16 rounded-[20px] sm:px-12 md:px-16 col-span-2 border border-white h-[600px]">
              <span className="block text-violet-400 font-monument tracking-[2px] text-7xl mt-2">BioFold</span>
              <p className="mb-8 mt-4 text-lg">Analyze protein sequences and structures with tools for   
                <span className="font-medium text-gray-50"> hydrophobicity, motifs, and 3D prediction</span>
              </p>
              
              <form noValidate="" action="" className="self-stretch space-y-3">
                  <select 
                  className="w-full bg-black p-2 mb-4 border border-gray-300 rounded-[10px] focus:outline-none focus:ring-2 focus:ring-blue-500"
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
                  className="w-full bg-black p-2 mb-4 border border-gray-300 rounded-[10px] h-32 focus:outline-none focus:ring-2 focus:ring-blue-500"
                  value={sequence}
                  onChange={(e) => setSequence(e.target.value)}
                  placeholder="Enter protein sequence..."
                />
                <Popup analysisOptions={analysisOptions} setAnalysisOptions={setAnalysisOptions} handleAnalyze={handleAnalyze}/>
                {error && (
                  <div className="text-red-700 px-4 py-3 rounded mb-4 text-center">
                    {error}
                  </div>
                )}
              </form>
  
            </div>
            {/* <img src="https://t3.ftcdn.net/jpg/02/65/31/80/360_F_265318078_VHJWZnzly32CcH3tUWb9ruOLf8lBRJOn.jpg" alt="" className="object-cover w-[100%] h-[600px] rounded-md xl:col-span-3 bg-gray-500" /> */}
                <div className='border border-white p-2 flex justify-center items-center w-[100%] h-[600px] col-span-3 rounded-[20px]'><LandModel/></div>
          </div>
        </section>
      :
      <Loading/>
      :
      <>
      {results && (
        <MainResult results={results} sequence={sequence}/>
      )}
        
      </>
      }
      
    </div>
  )
}
