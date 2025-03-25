import React from 'react';

const ExampleMolecules = ({ examples, onSelect }) => {
  if (!examples || Object.keys(examples).length === 0) return null;

  return (
    <div className="bg-white rounded-lg shadow-md p-6">
      <h2 className="text-xl font-semibold mb-4">Example Molecules</h2>
      <div className="space-y-2">
        {Object.entries(examples).map(([name, smiles]) => (
          <button
            key={name}
            onClick={() => onSelect(smiles)}
            className="w-full text-left px-4 py-2 border border-gray-200 rounded-md hover:bg-blue-50 hover:border-blue-300 transition duration-200"
          >
            <div className="font-medium">{name}</div>
            <div className="text-sm text-gray-500 truncate">{smiles}</div>
          </button>
        ))}
      </div>
    </div>
  );
};

export default ExampleMolecules;