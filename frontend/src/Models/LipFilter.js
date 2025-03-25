import React, { useState } from 'react';
import MolecularPropertyFilter from './Drug_DiscComponents/Molecular';
import DescriptorDistribution from './Drug_DiscComponents/Descriptor';
import ChemicalSpace from './Drug_DiscComponents/Chemical';

const DrugDiscovery = () => {
  const [activeTab, setActiveTab] = useState('filter');

  return (
    <div className="min-h-screen bg-gray-100">
      <nav className="bg-white shadow-lg p-4">
        <div className="max-w-7xl mx-auto">
          <div className="flex flex-col sm:flex-row justify-between items-center">
            <span className="text-xl font-bold mb-4 sm:mb-0">Drug Discovery Toolkit 🔬</span>
            <div className="flex space-x-4">
              <button
                onClick={() => setActiveTab('filter')}
                className={`px-4 py-2 rounded ${
                  activeTab === 'filter'
                    ? 'bg-blue-500 text-white'
                    : 'bg-gray-200 hover:bg-gray-300'
                }`}
              >
                Molecular Property Filter
              </button>
              <button
                onClick={() => setActiveTab('distribution')}
                className={`px-4 py-2 rounded ${
                  activeTab === 'distribution'
                    ? 'bg-blue-500 text-white'
                    : 'bg-gray-200 hover:bg-gray-300'
                }`}
              >
                Descriptor Distribution
              </button>
              <button
                onClick={() => setActiveTab('space')}
                className={`px-4 py-2 rounded ${
                  activeTab === 'space'
                    ? 'bg-blue-500 text-white'
                    : 'bg-gray-200 hover:bg-gray-300'
                }`}
              >
                Chemical Space
              </button>
            </div>
          </div>
        </div>
      </nav>

      <main className="max-w-7xl mx-auto p-6">
        {activeTab === 'filter' && <MolecularPropertyFilter />}
        {activeTab === 'distribution' && <DescriptorDistribution />}
        {activeTab === 'space' && <ChemicalSpace />}
      </main>
    </div>
  );
};

export default DrugDiscovery;