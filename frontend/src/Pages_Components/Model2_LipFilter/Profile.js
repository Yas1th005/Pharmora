import React from 'react'

export default function Profile({compound}) {
  return (
    <div>
        <div className="p-4 sm:flex sm:space-x-6 border border-white text-white w-full h-[250px] overflow-hidden">
            <div className="flex-shrink-0 w-full h-44 sm:h-40 sm:w-40 my-auto">
                <img src={`data:image/png;base64,${compound.image}`} alt="" className="object-cover object-center w-full h-full rounded dark:bg-gray-500" />
            </div>
            <div className="flex flex-col space-y-4 my-auto">
                <div>
                    <h2 className="text-2xl font-semibold">{compound.name}</h2>
                    <span className="text-sm dark:text-gray-600">{compound.smiles}</span>
                </div>
                <div className="space-y-1">
                    <span className="flex items-center space-x-2">
                        <span className="dark:text-gray-600">MW: {compound['Molecular Weight']}</span>
                    </span>
                    <span className="flex items-center space-x-2">
                        <span className="dark:text-gray-600">LogP: {compound.LogP}</span>
                    </span>
                    <span className="flex items-center space-x-2">
                        <span className="dark:text-gray-600">TPSA: {compound['Topological Polar Surface Area']}</span>
                    </span>
                    <span className="flex items-center space-x-2">
                        <span className="dark:text-gray-600">HBD: {compound['H-Bond Donors']} | HBA: {compound['H-Bond Acceptors']}</span>
                    </span>
                </div>
            </div>
        </div>
    </div>
  )
}
