import React, { useState, useEffect } from 'react';

const MolecularPropertyFilter = () => {
  const [compounds, setCompounds] = useState([]);
  const [filters, setFilters] = useState({
    mw_min: 200,
    mw_max: 500,
    logp_min: -0.4,
    logp_max: 5.5,
    tpsa_min: 0,
    tpsa_max: 140,
    hbd_min: 0,
    hbd_max: 5,
    hba_min: 0,
    hba_max: 10,
    search: '',
  });

  useEffect(() => {
    const fetchCompounds = async () => {
      const params = new URLSearchParams(filters);
      const response = await fetch(`http://localhost:5000/api/compounds?${params}`);
      const data = await response.json();
      setCompounds(data);
    };

    fetchCompounds();
  }, [filters]);

  return (
    <div className="space-y-6">
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        <div className="bg-white p-6 rounded-lg shadow">
          <h3 className="text-lg font-medium mb-4">Molecular Weight Range</h3>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="0"
              max="1000"
              value={filters.mw_min}
              onChange={(e) => setFilters(prev => ({ ...prev, mw_min: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.mw_min}</span>
          </div>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="0"
              max="1000"
              value={filters.mw_max}
              onChange={(e) => setFilters(prev => ({ ...prev, mw_max: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.mw_max}</span>
          </div>
        </div>

        <div className="bg-white p-6 rounded-lg shadow">
          <h3 className="text-lg font-medium mb-4">LogP Range</h3>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="-5"
              max="10"
              step="0.1"
              value={filters.logp_min}
              onChange={(e) => setFilters(prev => ({ ...prev, logp_min: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.logp_min}</span>
          </div>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="-5"
              max="10"
              step="0.1"
              value={filters.logp_max}
              onChange={(e) => setFilters(prev => ({ ...prev, logp_max: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.logp_max}</span>
          </div>
        </div>

        <div className="bg-white p-6 rounded-lg shadow">
          <h3 className="text-lg font-medium mb-4">TPSA Range</h3>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="0"
              max="200"
              value={filters.tpsa_min}
              onChange={(e) => setFilters(prev => ({ ...prev, tpsa_min: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.tpsa_min}</span>
          </div>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="0"
              max="200"
              value={filters.tpsa_max}
              onChange={(e) => setFilters(prev => ({ ...prev, tpsa_max: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.tpsa_max}</span>
          </div>
        </div>

        <div className="bg-white p-6 rounded-lg shadow">
          <h3 className="text-lg font-medium mb-4">H-Bond Donors</h3>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="0"
              max="15"
              value={filters.hbd_min}
              onChange={(e) => setFilters(prev => ({ ...prev, hbd_min: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.hbd_min}</span>
          </div>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="0"
              max="15"
              value={filters.hbd_max}
              onChange={(e) => setFilters(prev => ({ ...prev, hbd_max: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.hbd_max}</span>
          </div>
        </div>

        <div className="bg-white p-6 rounded-lg shadow">
          <h3 className="text-lg font-medium mb-4">H-Bond Acceptors</h3>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="0"
              max="20"
              value={filters.hba_min}
              onChange={(e) => setFilters(prev => ({ ...prev, hba_min: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.hba_min}</span>
          </div>
          <div className="flex items-center space-x-4">
            <input
              type="range"
              min="0"
              max="20"
              value={filters.hba_max}
              onChange={(e) => setFilters(prev => ({ ...prev, hba_max: Number(e.target.value) }))}
              className="w-full"
            />
            <span>{filters.hba_max}</span>
          </div>
        </div>
      </div>

      <input
        type="text"
        placeholder="Search compounds..."
        className="w-full max-w-md px-4 py-2 border rounded"
        value={filters.search}
        onChange={(e) => setFilters(prev => ({ ...prev, search: e.target.value }))}
      />

      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
        {compounds.map((compound, index) => (
          <div key={index} className="bg-white p-4 rounded-lg shadow">
            <div className="aspect-square mb-4">
              <img
                src={`data:image/png;base64,${compound.image}`}
                alt={compound.name}
                className="w-full h-full object-contain"
              />
            </div>
            <h3 className="text-lg font-medium mb-2">{compound.name}</h3>
            <p className="text-sm text-gray-500 mb-1">MW: {compound['Molecular Weight']}</p>
            <p className="text-sm text-gray-500 mb-1">LogP: {compound.LogP}</p>
            <p className="text-sm text-gray-500 mb-1">TPSA: {compound['Topological Polar Surface Area']}</p>
            <p className="text-sm text-gray-500 mb-1">
              HBD: {compound['H-Bond Donors']} | HBA: {compound['H-Bond Acceptors']}
            </p>
            <p className="text-sm text-gray-500 break-all">{compound.smiles}</p>
          </div>
        ))}
      </div>
    </div>
  );
};

export default MolecularPropertyFilter;