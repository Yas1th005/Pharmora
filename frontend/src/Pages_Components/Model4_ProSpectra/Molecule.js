import React, { useEffect, useRef } from 'react';

const MoleculeViewer = ({ result }) => {
  const structure3DRef = useRef(null);

  useEffect(() => {
    // Handle 3D structure rendering
    if (result && result['3d_structure'] && result['3d_structure']['view_html'] && structure3DRef.current) {
      // Extract the script from view_html and execute it
      const htmlContent = result['3d_structure']['view_html'];
      structure3DRef.current.innerHTML = htmlContent;
      
      // Execute any scripts in the HTML content
      const scripts = structure3DRef.current.getElementsByTagName('script');
      for (let i = 0; i < scripts.length; i++) {
        const script = document.createElement('script');
        script.textContent = scripts[i].textContent;
        document.body.appendChild(script);
        document.body.removeChild(script);
      }
    }
  }, [result]);

  if (!result || !result.visualizations) return null;

  return (
    <div className="bg-black text-white rounded-lg shadow-md p-6">
      <h2 className="text-xl font-semibold mb-4">Molecule Visualization</h2>
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        {/* 2D Structure */}
        {result.visualizations['2D Structure'] && (
          <div className="border rounded-lg p-4">
            <h3 className="text-lg font-medium mb-2">2D Structure</h3>
            <div className="flex justify-center">
              <img 
                src={`data:image/png;base64,${result.visualizations['2D Structure']}`} 
                alt="2D Molecule Structure" 
                className="max-w-full h-auto"
              />
            </div>
          </div>
        )}
        
        {/* 2D Detailed Structure */}
        {result.visualizations['2D Detailed Structure'] && (
          <div className="border rounded-lg p-4">
            <h3 className="text-lg font-medium mb-2">2D Detailed Structure</h3>
            <div className="flex justify-center">
              <img 
                src={`data:image/png;base64,${result.visualizations['2D Detailed Structure']}`} 
                alt="2D Detailed Molecule Structure" 
                className="max-w-full h-auto"
              />
            </div>
          </div>
        )}
        
        {/* 3D Structure */}
        {result['3d_structure'] && !result['3d_structure']['error'] && (
          <div className="border rounded-lg p-4 md:col-span-2">
            <h3 className="text-lg font-medium mb-2">3D Structure</h3>
            <div ref={structure3DRef} className="flex justify-center h-fit"></div>
          </div>
        )}
      </div>
    </div>
  );
};

export default MoleculeViewer;