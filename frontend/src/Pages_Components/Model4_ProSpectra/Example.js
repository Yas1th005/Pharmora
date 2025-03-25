import React from 'react';

const ExampleMolecules = ({ examples, onSelect }) => {
  if (!examples || Object.keys(examples).length === 0) return null;

  return (
    <div className="bg-black text-white rounded-lg shadow-md p-2 h-[300px] col-span-1 border border-white flex flex-col justify-center items-center">
      <h2 className="text-base font-poppins text-center font-semibold mb-4">Example Molecules</h2>
      <div className="grid grid-cols-2 gap-6">
        {Object.entries(examples).map(([name, smiles]) => (
          <button
            key={name}
            onClick={() => onSelect(smiles)}
            className="w-[300px] text-left px-4 py-2 border border-gray-200 rounded-md hover:bg-blue-50 hover:border-blue-300 transition duration-200"
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