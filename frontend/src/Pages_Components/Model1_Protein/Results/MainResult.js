import React,{useState,useEffect,useRef} from 'react';
import { TabsBtn, TabsContent, TabsProvider } from '../Tab';
import { PieChart, LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer,Pie, Cell, Legend } from 'recharts';
import { ChevronDown, ChevronUp } from 'lucide-react';

function MolViewer({ pdbData }) {
  const viewerRef = useRef(null);
  
  useEffect(() => {
    if (pdbData && viewerRef.current) {
      // Initialize 3Dmol viewer
      const viewer = window.$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: 'black',
      });
      
      // Add the PDB data
      viewer.addModel(pdbData, "pdb");
      
      // Style the protein
      viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
      
      // Zoom to fit the protein
      viewer.zoom(0.4);
      
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
      className="rounded z-0"
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

  const [showAllRows, setShowAllRows] = useState(false);


  // Brighter, more vibrant colors for better contrast on dark background
  const COLORS = [
    '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEEAD',
    '#FFD93D', '#6C5CE7', '#A8E6CF', '#FF8B94', '#A8D8EA',
    '#FFD93D', '#6C5CE7', '#A8E6CF', '#FF8B94', '#A8D8EA',
    '#FFD93D', '#6C5CE7', '#A8E6CF', '#FF8B94', '#A8D8EA'
  ];

  const visibleRows = showAllRows ? counts : counts.slice(0, 5);

  return (
    <div className="bg-black p-6 rounded-lg shadow-xl space-y-6 border border-gray-800">
      <h2 className="text-2xl font-bold text-gray-100 text-center">
        Amino Acid Distribution
      </h2>

      {/* Donut Chart */}
      <div className="w-full h-[550px] bg-black rounded-lg p-4">
        <ResponsiveContainer width="100%" height="100%">
          <PieChart>
            <Pie
              data={counts}
              cx="50%"
              cy="50%"
              innerRadius="60%"
              outerRadius="80%"
              paddingAngle={2}
              dataKey="count"
              nameKey="aa"
              label={({ name, percent }) => `${name} (${(percent * 100).toFixed(1)}%)`}
              labelLine={{ stroke: '#666' }}
            >
              {counts.map((entry, index) => (
                <Cell 
                  key={`cell-${index}`} 
                  fill={COLORS[index]} 
                  className="transition-all duration-300 hover:opacity-80 hover:scale-105"
                />
              ))}
            </Pie>
            <Tooltip 
              content={({ active, payload }) => {
                if (active && payload && payload.length) {
                  const data = payload[0].payload;
                  return (
                    <div className="bg-gray-900 p-3 rounded shadow-lg border border-gray-700">
                      <p className="text-gray-100 font-medium">{`${data.aa}`}</p>
                      <p className="text-gray-300">{`Count: ${data.count}`}</p>
                      <p className="text-gray-300">{`Percentage: ${data.percentage.toFixed(1)}%`}</p>
                    </div>
                  );
                }
                return null;
              }}
            />
            <Legend 
              layout="horizontal" 
              verticalAlign="bottom" 
              align="center"
              wrapperStyle={{ 
                color: '#fff',
                padding: '20px 0'
              }}
              formatter={(value) => <span className="text-gray-300">{value}</span>}
            />
          </PieChart>
        </ResponsiveContainer>
      </div>

      {/* Table Section */}
      <div className="mt-8 w-[80vw] mx-auto">
        <div className="overflow-x-hidden bg-black rounded-b-lg border-t border-gray-800">
          <table className="min-w-full divide-y divide-gray-800">
            <thead>
              <tr className="bg-gray-900">
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-300 uppercase tracking-wider">
                  Amino Acid
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-300 uppercase tracking-wider">
                  Count
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-300 uppercase tracking-wider">
                  Percentage
                </th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-800">
              {visibleRows.map(({ aa, count, percentage }, index) => (
                <tr 
                  key={aa}
                  className="transform transition-all duration-200 hover:bg-gray-900 hover:scale-[1.01]"
                >
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div className="flex items-center">
                      <div 
                        className="w-3 h-3 rounded-full mr-2"
                        style={{ backgroundColor: COLORS[index] }}
                      />
                      <span className="text-gray-300">{aa}</span>
                    </div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-gray-300">
                    {count}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-gray-300">
                    {percentage.toFixed(2)}%
                  </td>
                </tr>
              ))}
            </tbody>
          </table>

          {/* View More/Less Button */}
          <div className="p-4 border-t border-gray-800">
            <button
              onClick={() => setShowAllRows(!showAllRows)}
              className="w-full flex items-center justify-center px-4 py-2 bg-gray-900 text-gray-300 
                         rounded-lg hover:bg-gray-800 transition-colors duration-200 space-x-2"
            >
              <span>{showAllRows ? 'Show Less' : 'View Full Table'}</span>
              {showAllRows ? (
                <ChevronUp className="w-4 h-4" />
              ) : (
                <ChevronDown className="w-4 h-4" />
              )}
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}


function MainResult({results,sequence}) {
  console.log(results)
  return (
    <div className="p-4 bg-black relative h-screen w-[100vw]">
      <TabsProvider defaultValue={'design'}>
      <div className="flex justify-center mt-2">
          <div className="flex items-center w-fit bg-transparent p-1 text-white rounded-[10px] border">
            <TabsBtn value="design">
              <span className="relative z-[2] uppercase">Structure</span>
            </TabsBtn>
            <TabsBtn value="collaborate">
              <span className="relative z-[2] uppercase">amino acid</span>
            </TabsBtn>
            <TabsBtn value="share">
              <span className="relative z-[2] uppercase">Properties</span>
            </TabsBtn>
            <TabsBtn value="publish1">
              <span className="relative z-[2] uppercase">advanced</span>
            </TabsBtn>
          </div>
        </div>
        
        <TabsContent value="design">
          <div className="w-[80vw] mx-auto mt-[50px]">
            <section className="border border-white text-white mt-6 rounded-[20px]">
              <div className="flex flex-col-reverse mx-auto lg:flex-row">
                <div className="flex flex-col px-6 py-8 space-y-6 rounded-sm sm:p-8 lg:p-12 lg:w-1/2 xl:w-2/5 bg-transparent">
                  <div className="flex space-x-1 sm:space-x-3">
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor" className="flex-shrink-0 w-6 h-6">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M9 12l2 2 4-4M7.835 4.697a3.42 3.42 0 001.946-.806 3.42 3.42 0 014.438 0 3.42 3.42 0 001.946.806 3.42 3.42 0 013.138 3.138 3.42 3.42 0 00.806 1.946 3.42 3.42 0 010 4.438 3.42 3.42 0 00-.806 1.946 3.42 3.42 0 01-3.138 3.138 3.42 3.42 0 00-1.946.806 3.42 3.42 0 01-4.438 0 3.42 3.42 0 00-1.946-.806 3.42 3.42 0 01-3.138-3.138 3.42 3.42 0 00-.806-1.946 3.42 3.42 0 010-4.438 3.42 3.42 0 00.806-1.946 3.42 3.42 0 013.138-3.138z"></path>
                    </svg>
                    <div className="space-y-2">
                      <p className="text-xl font-bold leading-snug font-poppins tracking-widest">Structural Stability</p>
                      <p className="leading-snug text-justify">The 3D folding of a protein structure is stabilized by hydrogen bonds, hydrophobic interactions, and disulfide bridges, ensuring proper function.</p>
                    </div>
                  </div>
                  <div className="flex space-x-2 sm:space-x-3">
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor" className="flex-shrink-0 w-6 h-6">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M9 12l2 2 4-4M7.835 4.697a3.42 3.42 0 001.946-.806 3.42 3.42 0 014.438 0 3.42 3.42 0 001.946.806 3.42 3.42 0 013.138 3.138 3.42 3.42 0 00.806 1.946 3.42 3.42 0 010 4.438 3.42 3.42 0 00-.806 1.946 3.42 3.42 0 01-3.138 3.138 3.42 3.42 0 00-1.946.806 3.42 3.42 0 01-4.438 0 3.42 3.42 0 00-1.946-.806 3.42 3.42 0 01-3.138-3.138 3.42 3.42 0 00-.806-1.946 3.42 3.42 0 010-4.438 3.42 3.42 0 00.806-1.946 3.42 3.42 0 013.138-3.138z"></path>
                    </svg>
                    <div className="space-y-2 mb-5">
                      <p className="text-xl font-bold leading-snug font-poppins tracking-widest">Functional Specificity</p>
                      <p className="leading-snug text-justify">A protein’s structure dictates its interaction with molecules, enabling precise biological processes like enzyme activity and signal transduction.</p>
                    </div>
                  </div>
                  <div className="mt-5 border border-white p-3 rounded-[10px]">
                        <p className="text-lg text-center font-poppins">Confidence Score (pIDDT):</p>
                        <p className="text-2xl text-center">{results.structure.b_value}/1.00</p>
                      </div>
                      <button
                        className="mt-4 rounded-[10px] bg-violet-400 text-gray-900 py-2 px-4 font-poppins font-bold"
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
                <div className="lg:w-1/2 xl:w-3/5 dark:bg-gray-100">
                  <div className="flex items-center justify-center p-2 md:p-8 lg:p-6">
                  {results.structure && !results.structure.error && (
                    <div className="bg-black rounded-lg shadow p-6 border border-white text-white w-full rounded-[10px]">
                      <h2 className="text-xl font-bold mb-4 text-center font-poppins">Protein Structure</h2>
                      <MolViewer pdbData={results.structure.pdb_string} />
                    </div>
                  )}
                  </div>
                </div>
              </div>
            </section>
          </div>
        </TabsContent>
        <TabsContent value="collaborate">
          <div className="w-full bg-black overflow-x-hidden">
            <div className="bg-black text-white rounded-lg shadow p-6">
              <AminoAcidComposition sequence={sequence} />
            </div>
          </div>
        </TabsContent>
        <TabsContent value="share">
          <div className="w-full mt-6">
          <div className="bg-black text-white rounded-lg shadow p-6">
              <h2 className="text-3xl font-bold mb-6 font-poppins tracking-widest text-center">Molecular Properties</h2>
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4 w-[80vw] mx-auto">
                {Object.entries(results.molecular_properties).map(([key, value]) => (
                  <div key={key} className="bg-black border border-white text-white p-4 rounded-[20px]">
                    <h3 className="font-semibold capitalize text-center text-xl mb-2">{key.replace(/_/g, ' ')}</h3>
                    <Data name={key.replace(/_/g, ' ')}/>
                    <p className="text-2xl text-center">{typeof value === 'number' ? value.toFixed(2) : value} <Unit unit={key.replace(/_/g, ' ')}/></p>
                  </div>
                ))}
              </div>
            </div>
            
            {results.hydrophobicity && (
                          <div className="bg-black text-white rounded-lg shadow p-6">
                            <h2 className="text-3xl font-bold mb-6 font-poppins tracking-widest text-center">Hydrophobicity Profile</h2>
                            <div className="h-[280px]">
                              <ResponsiveContainer width="100%" height="100%">
                                <LineChart
                                  data={results.hydrophobicity.map((value, index) => ({
                                    position: index,
                                    value
                                  }))}
                                >
                                  <CartesianGrid strokeDasharray="7 9" />
                                  <XAxis dataKey="position" />
                                  <YAxis />
                                  <Tooltip />
                                  <Line type="monotone" dataKey="value" stroke="#fff" />
                                </LineChart>
                              </ResponsiveContainer>
                            </div>
                          </div>
                        )}
          </div>
        </TabsContent>
        <TabsContent value="publish1">
          <div className="w-full mt-8">
          {results.advanced_analysis && (
              <div className="bg-transparent text-white border border-white rounded-[20px] shadow px-9 py-7 w-[90vw] mx-auto">
                <h2 className="text-3xl font-bold mb-4 text-center font-poppins tracking-widest">Advanced Analysis</h2>
                
                {/* Motifs */}
                <div className="mb-6">
                  <h3 className="mb-5 font-poppins tracking-widest text-xl">Detected Motifs:</h3>
                  <div className='grid grid-cols-1 md:grid-cols-2 lg:grid-cols-2 gap-x-4 gap-y-1'>
                  {Object.entries(results.advanced_analysis.motifs).map(([motif, positions]) => (
                    <div key={motif} className="mb-4 bg-transparent border border-white p-3 rounded-[10px] h-[100px] flex flex-col justify-center items-center">
                      <p className="font-medium text-center uppercase">{motif}</p>
                      <p className="text-gray-300 text-center">
                        {positions.length > 0 ? 
                          `Found at positions: ${positions.join(', ')}` : 
                          'None found'}
                      </p>
                    </div>
                  ))}
                  </div>             
                </div>

                {/* Secondary Structure */}
                <div>
                  <h3 className="mb-2 font-poppins tracking-widest text-2xl">Secondary Structure Composition:</h3>
                  {Object.entries(results.advanced_analysis.secondary_structure).map(([structure, value]) => (
                    <div key={structure} className="mb-3">
                      <div className="flex items-center justify-between mb-1">
                        <span className="font-medium">{structure}</span>
                        <span>{(value * 100).toFixed(1)}%</span>
                      </div>
                      <div className="w-full bg-transparent border border-white rounded-full h-2">
                        <div
                          className="bg-white rounded-full h-2"
                          style={{ width: `${value * 100}%` }}
                        />
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        </TabsContent>
      </TabsProvider>
    </div>
  );
}
export default MainResult;


const Data = ({name}) => {
  if (name === "instability index") {
    return <p className="text-sm text-center mb-2">Predicts the stability of a protein in a test tube; >40 indicates instability</p>;
  } else if (name === "isoelectric point") { // Removed (pI) for proper string matching
    return <p className="text-sm text-center mb-2">The pH at which a protein has no net charge and is least soluble</p>;
  } else if(name === "molecular weight") {
    return <p className="text-sm text-center mb-2">The total mass of a protein, calculated based on its amino acid sequence</p>; // Default case
  }else{
    return <p className="text-sm text-center mb-2">The overall electrical charge of a protein at a given pH</p>;
  }
};
const Unit = ({unit}) => {
  if (unit === "instability index") {
    return <span className="text-lg "></span>;
  } else if (unit === "isoelectric point") { // Removed (pI) for proper string matching
    return <span className="text-lg ">pH</span>;
  } else if(unit === "molecular weight") {
    return <span className="text-lg ">g/mol</span>; // Default case
  }else{
    return <span className="text-lg ">C</span>;
  }
};
