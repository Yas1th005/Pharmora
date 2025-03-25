import React from 'react';
import M1_Protein from './Pages/M1_Protein';
import Docking from './Docker/Docking';
import Main from './ProSpectra/Main';
import DrugVista from './DrugVista/Main';
import Final from './Final/Final';

const App = () => {

  return (
    <div>
      {/* <Docking/> */}
      {/* <Main/> */}
      {/* <DrugVista/> */}
      <Final/>
    </div>
  );
};

export default App;