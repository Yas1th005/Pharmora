import React,{useEffect,useRef} from 'react';
import { TabsBtn, TabsContent, TabsProvider } from '../Tab';
import Plot from 'react-plotly.js';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';

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
      className="border border-white rounded z-0"
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
              <th className="px-6 py-3 bg-black text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Amino Acid
              </th>
              <th className="px-6 py-3 bg-black text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Count
              </th>
              <th className="px-6 py-3 bg-black text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                Percentage
              </th>
            </tr>
          </thead>
          <tbody className="bg-black divide-y divide-gray-200">
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

function MainResult({results,sequence}) {
  return (
    <div className="p-4 bg-black relative h-screen w-[100vw]">
      <TabsProvider defaultValue={'design'}>
        <div className="flex justify-center mt-2">
          <div className="flex items-center w-fit bg-transparent p-1 text-white rounded-md border">
            <TabsBtn value="design">
              <span className="relative z-[2] uppercase">design</span>
            </TabsBtn>
            <TabsBtn value="collaborate">
              <span className="relative z-[2] uppercase">collaborate</span>
            </TabsBtn>
            <TabsBtn value="share">
              <span className="relative z-[2] uppercase">share</span>
            </TabsBtn>
            <TabsBtn value="publish">
              <span className="relative z-[2] uppercase">publish</span>
            </TabsBtn>
            <TabsBtn value="publish1">
              <span className="relative z-[2] uppercase">publish</span>
            </TabsBtn>
            <TabsBtn value="publish2">
              <span className="relative z-[2] uppercase">publish</span>
            </TabsBtn>
          </div>
        </div>

        <TabsContent value="design">
          <div className="w-[80vw] mx-auto">
          {results.structure && !results.structure.error && (
              <div className="bg-black rounded-lg shadow p-6 border border-white text-white">
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
          </div>
        </TabsContent>
        <TabsContent value="collaborate">
          <div className="w-full">
            <div className="bg-black text-white rounded-lg shadow p-6">
              <h2 className="text-xl font-bold mb-4">Amino Acid Composition</h2>
              <AminoAcidComposition sequence={sequence} />
            </div>
          </div>
        </TabsContent>
        <TabsContent value="share">
          <div className="w-full">
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
          </div>
        </TabsContent>
        <TabsContent value="publish">
          <div className="w-full">
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
          </div>
        </TabsContent>
        <TabsContent value="publish1">
          <div className="w-full">
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
          </div>
        </TabsContent>
        <TabsContent value="publish2">
          <div className="w-full">
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
        </TabsContent>
      </TabsProvider>
    </div>
  );
}
export default MainResult;
