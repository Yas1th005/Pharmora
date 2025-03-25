import React, { useState, useEffect } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';

const DescriptorDistribution = () => {
  const [descriptor, setDescriptor] = useState('Molecular Weight');
  const [data, setData] = useState(null);

  const descriptors = [
    'Molecular Weight',
    'LogP',
    'H-Bond Donors',
    'H-Bond Acceptors',
    'Topological Polar Surface Area'
  ];

  useEffect(() => {
    const fetchData = async () => {
      const response = await fetch(`http://localhost:5000/api/descriptor-stats?descriptor=${descriptor}`);
      const result = await response.json();
      
      // Create histogram data
      const values = result.values;
      const min = Math.floor(result.min);
      const max = Math.ceil(result.max);
      const binSize = (max - min) / 30;
      
      const bins = Array.from({ length: 30 }, (_, i) => ({
        x: min + (i * binSize),
        count: values.filter(v => v >= min + (i * binSize) && v < min + ((i + 1) * binSize)).length
      }));
      
      setData(bins);
    };

    fetchData();
  }, [descriptor]);

  return (
    <div className="space-y-6">
      <select
        value={descriptor}
        onChange={(e) => setDescriptor(e.target.value)}
        className="px-4 py-2 border rounded"
      >
        {descriptors.map(d => (
          <option key={d} value={d}>{d}</option>
        ))}
      </select>

      <div className="bg-white p-6 rounded-lg shadow">
        <div className="h-[400px]">
          {data && (
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis 
                  dataKey="x" 
                  label={{ value: descriptor, position: 'bottom' }}
                />
                <YAxis 
                  label={{ value: 'Count', angle: -90, position: 'insideLeft' }}
                />
                <Tooltip />
                <Line 
                  type="monotone" 
                  dataKey="count" 
                  stroke="#3b82f6" 
                  fill="#3b82f6" 
                />
              </LineChart>
            </ResponsiveContainer>
          )}
        </div>
      </div>
    </div>
  );
};

export default DescriptorDistribution;