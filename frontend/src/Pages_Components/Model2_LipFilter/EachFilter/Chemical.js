import React, { useState, useEffect } from 'react';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ZAxis } from 'recharts';
import Loading from '../../Model1_Protein/Loading';

const ChemicalSpace = () => {
  const [compounds, setCompounds] = useState([]);
  const [xAxis, setXAxis] = useState('LogP');
  const [yAxis, setYAxis] = useState('Molecular Weight');
  const [loading,setLoading]=useState(true);

  const descriptors = [
    'LogP',
    'Molecular Weight',
    'Topological Polar Surface Area',
    'H-Bond Donors',
    'H-Bond Acceptors'
  ];

  useEffect(() => {
    const fetchCompounds = async () => {
      const response = await fetch('http://localhost:5001/api/compounds');
      const data = await response.json();
      setCompounds(data);
      setLoading(false);
    };

    fetchCompounds();
  }, []);

  const CustomTooltip = ({ active, payload }) => {
    if (active && payload && payload.length) {
      const compound = payload[0].payload;
      return (
        <div className="bg-white p-4 shadow-lg rounded-lg border">
          <p className="font-medium">{compound.name}</p>
          <p className="text-sm text-gray-600">{`${xAxis}: ${compound[xAxis]}`}</p>
          <p className="text-sm text-gray-600">{`${yAxis}: ${compound[yAxis]}`}</p>
          <p className="text-sm text-gray-600">Aromatic Rings: {compound['Aromatic Ring Count']}</p>
          <p className="text-xs text-gray-500 mt-2 break-all">{compound.smiles}</p>
        </div>
      );
    }
    return null;
  };

  return (
    <div>
      {!loading?
      <div className="space-y-6">
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        <select
          value={xAxis}
          onChange={(e) => setXAxis(e.target.value)}
          className="px-4 py-2 border rounded bg-black text-white"
        >
          {descriptors.map(d => (
            <option key={d} value={d}>{d}</option>
          ))}
        </select>

        <select
          value={yAxis}
          onChange={(e) => setYAxis(e.target.value)}
          className="px-4 py-2 border rounded bg-black text-white"
        >
          {descriptors.map(d => (
            <option key={d} value={d}>{d}</option>
          ))}
        </select>
      </div>

      <div className="bg-black border border-white p-6 rounded-lg shadow">
        <div className="h-[600px]">
          <ResponsiveContainer width="100%" height="100%">
            <ScatterChart
              margin={{
                top: 20,
                right: 20,
                bottom: 50,
                left: 50,
              }}
            >
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis 
                dataKey={xAxis} 
                type="number" 
                name={xAxis}
                label={{ value: xAxis, position: 'bottom' }}
              />
              <YAxis 
                dataKey={yAxis} 
                type="number" 
                name={yAxis}
                label={{ value: yAxis, angle: -90, position: 'insideLeft' }}
              />
              <ZAxis dataKey="Aromatic Ring Count" range={[50, 400]} />
              <Tooltip content={<CustomTooltip />} />
              <Scatter
                name="Compounds"
                data={compounds}
                fill="#3b82f6"
                fillOpacity={0.6}
              />
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      </div>

      <div className="text-sm text-gray-500">
        Note: Bubble size represents the number of aromatic rings in the molecule
      </div>
    </div>
      :
      <div className='absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2'><Loading/></div>}
    </div>
  );
};

export default ChemicalSpace;