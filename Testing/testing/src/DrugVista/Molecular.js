import React, { useEffect, useRef } from 'react';

const MoleculeViewer = ({ smiles, width = "100%", height = "300px" }) => {
  const viewerRef = useRef(null);

  useEffect(() => {
    if (!smiles || !window.RDKit || !window.py3Dmol) return;

    // Clear any previous content
    if (viewerRef.current) {
      viewerRef.current.innerHTML = '';
    }

    try {
      // Convert SMILES to molecule
      const mol = window.RDKit.Molecule.fromSmiles(smiles);
      const molblock = mol.toMolBlock();
      
      // Create 3D visualization
      const viewer = window.py3Dmol.createViewer(viewerRef.current, {
        width: width,
        height: height,
        backgroundColor: 'white',
      });
      
      viewer.addModel(molblock, 'mol');
      viewer.setStyle({}, { stick: { radius: 0.2 } });
      viewer.zoomTo();
      viewer.spin(true);
      viewer.render();
      
      return () => {
        if (viewer) {
          viewer.spin(false);
        }
      };
    } catch (error) {
      console.error('Error rendering molecule:', error);
      if (viewerRef.current) {
        viewerRef.current.innerHTML = `<div class="p-4 text-red-500">Error rendering molecule: ${error.message}</div>`;
      }
    }
  }, [smiles, width, height]);

  return (
    <div className="molecule-viewer">
      {!smiles ? (
        <div className="flex items-center justify-center h-full bg-gray-100 rounded-md">
          <p className="text-gray-500">No molecule to display</p>
        </div>
      ) : (
        <div 
          ref={viewerRef} 
          className="border rounded-md shadow-sm"
          style={{ width, height }}
        />
      )}
    </div>
  );
};

export default MoleculeViewer;