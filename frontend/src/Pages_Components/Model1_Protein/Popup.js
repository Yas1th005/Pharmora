import { AnimatePresence, motion } from "framer-motion";
import { FiAlertCircle } from "react-icons/fi";
import { useState } from "react";
import options from "../../imgaes/options.png"

const Popup = ({ analysisOptions, setAnalysisOptions, handleAnalyze }) => {
  const [isOpen, setIsOpen] = useState(false);

  return (
    <div>
      <button
        type="button"
        onClick={() => setIsOpen(true)}
        className="w-full py-2 font-semibold rounded-[10px] bg-violet-400 text-gray-900"
      >
        Proceed
      </button>
      <SpringModal
        isOpen={isOpen}
        setIsOpen={setIsOpen}
        analysisOptions={analysisOptions}
        setAnalysisOptions={setAnalysisOptions}
        handleAnalyze={handleAnalyze}
      />
    </div>
  );
};

const SpringModal = ({ isOpen, setIsOpen, analysisOptions, setAnalysisOptions, handleAnalyze }) => {
  const handleCheckboxChange = (option) => {
    setAnalysisOptions((prevOptions) =>
      prevOptions.includes(option)
        ? prevOptions.filter((o) => o !== option)
        : [...prevOptions, option]
    );
  };

  return (
    <AnimatePresence>
      {isOpen && (
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          onClick={() => setIsOpen(false)}
          className="bg-slate-900/20 backdrop-blur p-8 fixed inset-0 z-50 grid place-items-center overflow-y-scroll cursor-pointer"
        >
          <motion.div
            initial={{ scale: 0, rotate: "12.5deg" }}
            animate={{ scale: 1, rotate: "0deg" }}
            exit={{ scale: 0, rotate: "0deg" }}
            onClick={(e) => e.stopPropagation()}
            className="bg-[#1E1E1E] text-white p-6 rounded-[20px] w-full max-w-lg shadow-xl cursor-default relative overflow-hidden"
          >
            {/* <FiAlertCircle className="text-white/10 rotate-12 text-[250px] absolute z-0 -top-24 -left-24" /> */}
            <div className="relative z-10">
              {/* <div className="w-16 h-16 mb-2 rounded-full text-3xl text-indigo-600 grid place-items-center mx-auto"> */}
                {/* <FiAlertCircle /> */}
                {/* <img src={options}></img> */}
              {/* </div> */}
              <p className="font-poppins text-white tracking-widest mt-1 mb-6 text-xl font-bold">Select Analysis option</p>
              {/* <h3 className="font-semibold mb-2">Analysis Options:</h3> */}
              <div className="space-y-2">
                {["Structure Prediction", "Sequence Composition", "Hydrophobicity"].map((option) => (
                  <label key={option} className="flex items-center space-x-2">
                    <input
                      type="checkbox"
                      className="form-checkbox"
                      checked={analysisOptions.includes(option)}
                      onChange={() => handleCheckboxChange(option)}
                    />
                    <span className="font-poppins text-lg">{option}</span>
                  </label>
                ))}
              </div>
              <div className="flex gap-2 mt-4">
                <button
                  type="button"
                  onClick={() => setIsOpen(false)}
                  className="bg-transparent hover:bg-white/10 transition-colors text-white font-semibold w-full py-2 rounded-[10px]"
                >
                  Change Sequence
                </button>
                <button
                  type="button"
                  onClick={() => {
                    setIsOpen(false);
                    handleAnalyze();
                  }}
                  className="hover:opacity-90 transition-opacity bg-violet-400 text-gray-900 font-semibold w-full py-2 rounded-[10px]"
                >
                  Analyze
                </button>
              </div>
            </div>
          </motion.div>
        </motion.div>
      )}
    </AnimatePresence>
  );
};

export default Popup;


// const check=()=>{
//   return(
//     <div class="container">
//       <input style="display: none;" id="cbx" type="checkbox" />
//       <label class="check" for="cbx">
//         <svg viewBox="0 0 18 18" height="18px" width="18px">
//           <path
//             d="M1,9 L1,3.5 C1,2 2,1 3.5,1 L14.5,1 C16,1 17,2 17,3.5 L17,14.5 C17,16 16,17 14.5,17 L3.5,17 C2,17 1,16 1,14.5 L1,9 Z"
//           ></path>
//           <polyline points="1 9 7 14 15 4"></polyline>
//         </svg>
//       </label>
//     </div>
//   )
// }