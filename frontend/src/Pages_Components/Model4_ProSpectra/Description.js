import React from 'react';

const DescriptorsTable = ({ descriptors }) => {
  if (!descriptors) return null;
  
  // Group descriptors by category
  const structuralProps = [
    'Molecular Formula', 'Molecular Weight', 'Heavy Atom Count', 
    'Total Atom Count', 'Total Bond Count'
  ];
  
  const chemicalProps = [
    'LogP', 'H-Bond Donors', 'H-Bond Acceptors', 
    'Rotatable Bonds', 'Topological Polar Surface Area'
  ];
  
  const structuralDescriptors = [
    'Ring Count', 'Aromatic Ring Count', 'Saturated Ring Count',
    'Heterocyclic Ring Count'
  ];
  
  const pharmacokinetic = [
    'Molar Refractivity', 'Fraction SP3 Carbons', 'Lipinski Violations'
  ];

  const renderDescriptorSection = (title, keys) => (
    <div className="mb-4">
      <h3 className="text-lg font-medium mb-2">{title}</h3>
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-2">
        {keys.map(key => (
          descriptors[key] !== undefined && (
            <div key={key} className="bg-black border border-white p-3 rounded">
              <span className="text-sm text-white block">{key}</span>
              <span className="font-medium">{descriptors[key]}</span>
            </div>
          )
        ))}
      </div>
    </div>
  );

  return (
    <div className="bg-black text-white border-white border rounded-lg shadow-md p-6">
      <h2 className="text-xl font-semibold mb-4">Molecular Descriptors</h2>
      
      {renderDescriptorSection('Structural Properties', structuralProps)}
      {renderDescriptorSection('Chemical Properties', chemicalProps)}
      {renderDescriptorSection('Structural Descriptors', structuralDescriptors)}
      {renderDescriptorSection('Pharmacokinetic Indicators', pharmacokinetic)}
    </div>
  );
};

export default DescriptorsTable;