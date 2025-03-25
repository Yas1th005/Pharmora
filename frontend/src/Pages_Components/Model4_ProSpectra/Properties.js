import React from 'react';

const PharmacologicalProperties = ({ properties }) => {
  if (!properties) return null;

  const getStatusColor = (property) => {
    const value = properties[property] || '';
    
    // Make sure value is treated as a string before using .includes()
    const valueStr = String(value);
    
    if (valueStr.includes('Excellent') || valueStr.includes('High') || valueStr.includes('Highly Suitable')) {
      return 'bg-green-100 text-green-800';
    } else if (valueStr.includes('Good') || valueStr.includes('Moderate')) {
      return 'bg-yellow-100 text-yellow-800';
    } else if (valueStr.includes('Poor') || valueStr.includes('Low') || valueStr.includes('Challenging')) {
      return 'bg-red-100 text-red-800';
    }
    return 'bg-gray-100 text-gray-800';
  };

  return (
    <div className="bg-black text-white rounded-lg shadow-md p-6">
      <h2 className="text-xl font-semibold mb-4">Pharmacological Properties</h2>
      
      <div className="space-y-4">
        {Object.entries(properties).map(([key, value]) => (
          <div key={key} className="border rounded-lg p-4">
            <h3 className="text-lg font-medium mb-2">{key.replace(/([A-Z])/g, ' $1').trim()}</h3>
            <div className={`inline-block px-3 py-1 rounded-full text-sm font-medium ${getStatusColor(key)}`}>
              {String(value)}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};

export default PharmacologicalProperties;