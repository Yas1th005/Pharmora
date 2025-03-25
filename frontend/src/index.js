import React from 'react';
import ReactDOM from 'react-dom/client';
import './index.css';
import App from './App';
import {createBrowserRouter, RouterProvider} from 'react-router-dom';
import M1_Protein from './Pages/M1_Protein';
import LinearCard from './Pages/M2_LipFilter';
import Docking from './Pages/M3_Docking';
import ProSpectra from './Pages/M4_ProSpectra';
import DrugVista from './Pages/M6_DrugVista';
import Final from './Pages/Final';

const root = ReactDOM.createRoot(document.getElementById('root'));

let allRoutes=createBrowserRouter(
  [
    {
      path:"/",
      element:<App/>
    },
    {
      path:"Model1",
      element:<M1_Protein/>
    },
    {
      path:"Model2",
      element:<LinearCard/>
    },
    {
      path:"Model3",
      element:<Docking/>
    },
    {
      path:"Model4",
      element:<ProSpectra/>
    },
    {
      path:"Model6",
      element:<DrugVista/>
    },
    {
      path:"Final",
      element:<Final/>
    }
  ]
)
root.render(
  <React.StrictMode>
    <RouterProvider router={allRoutes}/>
  </React.StrictMode>
);

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals