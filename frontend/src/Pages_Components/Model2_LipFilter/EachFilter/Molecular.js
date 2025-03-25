import React, { useState, useEffect } from 'react';
import Profile from '../Profile';
import Loading from '../../Model1_Protein/Loading';

const MolecularPropertyFilter = ({title,desc}) => {
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

  const [loading,setLoading]=useState(true);

  useEffect(() => {
    const fetchCompounds = async () => {
      const params = new URLSearchParams(filters);
      const response = await fetch(`http://localhost:5001/api/compounds?${params}`);
      const data = await response.json();
      setCompounds(data);
      setLoading(false);
    };

    fetchCompounds();
  }, [filters]);

  return (
    <div>
      {!loading?
        <div className="space-y-6">
        <section className="p-6 bg-black text-gray-100 mt-9">
          <div className="container grid gap-6 mx-auto text-center lg:grid-cols-2 xl:grid-cols-5 w-full">
            <div className="w-full px-6 py-16 rounded-md sm:px-12 md:px-16 xl:col-span-2 border border-white flex flex-col justify-center">
            <span className="block mb-2 text-violet-400 uppercase font-monument tracking-[2px] text-base">PROSPECTRA</span>
              <h1 className="text-5xl font-extrabold text-gray-50">{title}</h1>
              <p className="my-8">
              <span className="font-medium text-gray-50">Modular and versatile.</span>Fugit vero facilis dolor sit neque cupiditate minus esse accusamus cumque at.
              </p>  
            </div>
            <div className='border border-white w-full w-full h-full rounded-md xl:col-span-3 p-4'>
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-white">
              <div className="bg-transparent border border-white p-6 rounded-lg shadow">
                <h3 className="text-lg font-medium mb-4 text-white">Molecular Weight Range</h3>
                <div className="flex items-center space-x-4">
                  <input
                    type="range"
                    min="0"
                    max="1000"
                    value={filters.mw_min}
                    onChange={(e) => setFilters(prev => ({ ...prev, mw_min: Number(e.target.value) }))}
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
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
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
                  />
                  <span>{filters.mw_max}</span>
                </div>

              </div>

              <div className="bg-black border border-white p-6 rounded-lg shadow">
                <h3 className="text-lg font-medium mb-4">LogP Range</h3>
                <div className="flex items-center space-x-4">
                  <input
                    type="range"
                    min="-5"
                    max="10"
                    step="0.1"
                    value={filters.logp_min}
                    onChange={(e) => setFilters(prev => ({ ...prev, logp_min: Number(e.target.value) }))}
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
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
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
                  />
                  <span>{filters.logp_max}</span>
                </div>
              </div>

              <div className="bg-black border border-white p-6 rounded-lg shadow">
                <h3 className="text-lg font-medium mb-4">TPSA Range</h3>
                <div className="flex items-center space-x-4">
                  <input
                    type="range"
                    min="0"
                    max="200"
                    value={filters.tpsa_min}
                    onChange={(e) => setFilters(prev => ({ ...prev, tpsa_min: Number(e.target.value) }))}
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
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
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
                  />
                  <span>{filters.tpsa_max}</span>
                </div>
              </div>

              <div className="bg-black border border-white p-6 rounded-lg shadow">
                <h3 className="text-lg font-medium mb-4">H-Bond Donors</h3>
                <div className="flex items-center space-x-4">
                  <input
                    type="range"
                    min="0"
                    max="15"
                    value={filters.hbd_min}
                    onChange={(e) => setFilters(prev => ({ ...prev, hbd_min: Number(e.target.value) }))}
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
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
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
                  />
                  <span>{filters.hbd_max}</span>
                </div>
              </div>

              <div className="bg-black border border-white p-6 rounded-lg shadow">
                <h3 className="text-lg font-medium mb-4">H-Bond Acceptors</h3>
                <div className="flex items-center space-x-4">
                  <input
                    type="range"
                    min="0"
                    max="20"
                    value={filters.hba_min}
                    onChange={(e) => setFilters(prev => ({ ...prev, hba_min: Number(e.target.value) }))}
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
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
                    className="w-full appearance-none bg-transparent border border-white h-2 rounded-lg accent-white outline-none"
                  />
                  <span>{filters.hba_max}</span>
                </div>
              </div>
            </div>
            </div>
          </div>
        </section>

      <input
        type="text"
        placeholder="Search compounds..."
        className="w-[90vw] flex justify-center py-2 border rounded bg-black border border-white mx-auto px-2 text-center"
        value={filters.search}
        onChange={(e) => setFilters(prev => ({ ...prev, search: e.target.value }))}
      />

      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-2 gap-4 w-[90vw] mx-auto">
        {compounds.map((compound, index) => (
          // <div key={index} className="bg-black border border-white p-4 rounded-lg shadow">
          //   <div className="aspect-square mb-4">
          //     <img
          //       src={`data:image/png;base64,${compound.image}`}
          //       alt={compound.name}
          //       className="w-full h-full object-contain"
          //     />
          //   </div>
          //   <h3 className="text-base font-medium mb-2 text-white text-center">{compound.name}</h3>
          //   <p className="text-xs text-gray-500 mb-1">MW: {compound['Molecular Weight']}</p>
          //   <p className="text-xs text-gray-500 mb-1">LogP: {compound.LogP}</p>
          //   <p className="text-xs text-gray-500 mb-1">TPSA: {compound['Topological Polar Surface Area']}</p>
          //   <p className="text-xs text-gray-500 mb-1">
          //     HBD: {compound['H-Bond Donors']} | HBA: {compound['H-Bond Acceptors']}
          //   </p>
          //   <p className="text-xs text-gray-500 break-all">{compound.smiles}</p>
          // </div>
          <div key={index}><Profile compound={compound}/></div>
        ))}
      </div>
    </div>
      :
      <div className='absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2'><Loading/></div>}
    </div>
  );
};

export default MolecularPropertyFilter;