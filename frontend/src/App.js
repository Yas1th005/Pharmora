import logo from './logo.svg';
import './App.css';
import { Header } from './components/Header';
import Models from './components/Models';
import Footer from './components/Footer';
import LandModel from './components/LandModel';
import HeaderImg from './components/HeaderImg';
import Landing from './components/Landing';
import Intro from './components/Intro';
import Testimonials from './components/Testimonials';
import {motion} from "framer-motion"
import Faq from './components/Faqs';
import ProteinModel from './Models/Protien';
import M1_Protein from './Pages/M1_Protein';
import MainResult from './Pages_Components/Model1_Protein/Results/MainResult';
import DrugDiscovery from './Models/LipFilter';
import LinearCard from './Pages/M2_LipFilter';
// import Faq from './components/Faqs';
function App() {

  // const handleOpenStreamlit = (newTab = true) => {
  //   const url = "http://localhost:8501/Bioactivity_Predictor"; // URL for Streamlit frontend
  //   if (newTab) {
  //     window.open(url, "_blank"); // Opens in a new tab
  //   } else {
  //     window.location.href = url; // Opens in the same tab
  //   }
  // };

  // return (
  //   <div className="flex justify-center items-center h-screen bg-gray-100">
  //     <div className="text-center">
  //       <h1 className="text-2xl font-bold mb-6 text-gray-800">
  //         React and Streamlit Integration
  //       </h1>
  //       <div className="flex space-x-4">
  //         <button
  //           onClick={() => handleOpenStreamlit(true)}
  //           className="px-6 py-3 bg-blue-600 text-white font-semibold rounded-lg shadow hover:bg-blue-700 transition"
  //         >
  //           Open Streamlit in New Tab
  //         </button>
  //         <button
  //           onClick={() => handleOpenStreamlit(false)}
  //           className="px-6 py-3 bg-green-600 text-white font-semibold rounded-lg shadow hover:bg-green-700 transition"
  //         >
  //           Open Streamlit in Same Tab
  //         </button>
  //       </div>
  //     </div>
  //   </div>
  // );
  return (
    <div className="bg-black">

<Landing/>
<Intro/>
<motion.div
initial={{opacity:0, y:148}}
whileInView={{y:0,opacity:1}}
transition={{ease:"easeInOut",duration:2}}
viewport={{ once: true }}
className='mt-[150px]'
>
<Models/>
</motion.div>


<motion.div
initial={{opacity:0, y:148}}
whileInView={{y:0,opacity:1}}
transition={{ease:"easeInOut",duration:2}}
viewport={{ once: true}}
>
<Testimonials/>
</motion.div>

<motion.div
initial={{opacity:0, y:148}}
whileInView={{y:0,opacity:1}}
transition={{ease:"easeInOut",duration:2}}
viewport={{ once: true}}
>
<Faq/>
</motion.div>




<Footer/> 


{/* <ProteinModel/> */}
{/* <M1_Protein/> */}
{/* <DrugDiscovery/> */}

{/* <LinearCard/> */}
    </div>
  );
}

export default App;





