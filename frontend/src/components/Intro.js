import React from 'react'
import {motion} from 'framer-motion'
import ml from '../imgaes/intro/ml.jpg'
import sol from '../imgaes/intro/solubility.jpeg'


export default function Intro() {
  return (
    <div className='flex flex-col gap-[150px]'>
        <div
        className='flex flex-cols justify-center items-center p-6 w-[85vw] mx-auto'>
            <motion.div 
            initial={{opacity:0, y:148}}
            whileInView={{y:0,opacity:1}}
            transition={{ease:"easeInOut",duration:2}}
            viewport={{ once: true }}
            className='basis-2/5 border border-white p-5 rounded-[20px] ml-[-20px] mr-[20px]'><img src={ml} className='object-fill h-[75vh] w-full rounded-[10px]'/></motion.div>
            <motion.div
            initial={{opacity:0, y:148}}
            whileInView={{y:0,opacity:1}}
            transition={{ease:"easeInOut",duration:2}}
            viewport={{ once: true }}
            className='basis-3/5 text-white text-center'>
            <p className='font-poppins text-5xl font-bold ml-[50px] text-transparent inline-block pointer-events-none bg-gradient-to-b from-neutral-400 to-white bg-clip-text'>Unlock the future of drug discovery with <span className='text-[#FF0000]'>AI-powered</span> molecular scanning, feature analysis, and compound generation — all in one platform</p>
            {/* <p className='font-roboto mt-[50px]'>We are an all-in-one bioinformatics research platform designed to streamline drug discovery and molecular research. The platform offers advanced tools such as a Bioactivity Predictor (pIC50), New Molecule Discovery for Specific Targets, Solubility Predictor, DNA Nucleotide Counter, Antimicrobial Activity Predictor for Peptides, Molecular Descriptor Calculator, and a Lipinski’s Rule of Five Filter for Drugs. It is an ideal solution for researchers and scientists aiming to explore drug-like molecules, predict biological activity, and optimize lead compounds in a seamless, user-friendly interface.</p> */}
            </motion.div>
        </div>
        {/* <div className='flex flex-cols justify-center align-center  ml-[200px] mr-[200px]'>
        <p></p>
        <p className='basis-2/3'></p>
        <div className='basis-1/3 border border-white p-[5px]'><img src={sol}/></div> */}
        <div className='flex flex-cols justify-center items-center p-6 w-[85vw] mx-auto'>
            <motion.div
            initial={{opacity:0, y:148}}
            whileInView={{y:0,opacity:1}}
            transition={{ease:"easeInOut",duration:2}}
            viewport={{ once: true}}
            className='basis-3/5 text-white text-center'>
            <p className='font-poppins text-5xl font-bold text-transparent inline-block pointer-events-none bg-gradient-to-b from-neutral-400 to-white bg-clip-text ml-[-30px] mr-[30px]'>Uncover key substructures, decode fingerprints, and reveal hidden drug discovery insights with a single click.</p>
            {/* <p className='font-roboto mt-[50px]'>We are an all-in-one bioinformatics research platform designed to streamline drug discovery and molecular research. The platform offers advanced tools such as a Bioactivity Predictor (pIC50), New Molecule Discovery for Specific Targets, Solubility Predictor, DNA Nucleotide Counter, Antimicrobial Activity Predictor for Peptides, Molecular Descriptor Calculator, and a Lipinski’s Rule of Five Filter for Drugs. It is an ideal solution for researchers and scientists aiming to explore drug-like molecules, predict biological activity, and optimize lead compounds in a seamless, user-friendly interface.</p> */}
            </motion.div>
            <motion.div
            initial={{opacity:0, y:148}}
            whileInView={{y:0,opacity:1}}
            transition={{ease:"easeInOut",duration:2}}
            viewport={{ once: true}}
            className='basis-2/5 h-fit border border-white rounded-[20px] p-5 ml-[50px]'><img src={sol} className='object-fill w-full rounded-[10px]'/></motion.div>
            
        </div>
            
    </div>
        

  )
}
