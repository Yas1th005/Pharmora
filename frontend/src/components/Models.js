'use client';
import { ReactLenis } from 'lenis/react';
import { useTransform, motion, useScroll } from 'framer-motion';
import { useRef } from 'react';
import { Routes, Link } from "react-router"

import one from "../imgaes/models/1.jpeg"
import two from "../imgaes/models/2.jpg"
import three from "../imgaes/models/3.jpg"
import four from "../imgaes/models/4.jpg"
import five from "../imgaes/models/5.jpg"
import six from "../imgaes/models/6.jpg"





const projects = [
  {
    title1: 'BioFold',
    title2:'Revolutionizing Protein Analysis',
    description:
      'BioFold is an advanced protein analysis platform designed to provide researchers with comprehensive insights into protein sequences and structures. From sequence composition and hydrophobicity analysis to motif detection and 3D structure prediction, BioFold equips you with state-of-the-art tools for in-depth exploration.',
    src: one,
    link: 'https://images.unsplash.com/photo-1605106702842-01a887a31122?q=80&w=500&auto=format&fit=crop',
    color: '#000000',
    path:`Model1`
  },
  {
    title1: 'DrugVista',
    title2:'Insightful Drug Discovery Toolkit',
    description:
      'DrugVista is an advanced computational platform designed to empower researchers and scientists in the field of drug discovery. Leveraging cutting-edge cheminformatics and machine learning tools, DrugVista provides a comprehensive suite of features to analyze, visualize, and interpret molecular properties and chemical space.',
    src: two,
    link: 'https://images.unsplash.com/photo-1605106250963-ffda6d2a4b32?w=500&auto=format&fit=crop&q=60',
    color: '#000000',
    path:'Model6'
  },
  {
    title1: 'ProSpectra',
    title2:'Comprehensive Molecular Analysis',
    description:
      'ProSpectra is a state-of-the-art molecular analysis platform designed to empower chemists, researchers, and students. By integrating advanced cheminformatics and visualization tools, ProSpectra provides a robust suite of features to explore, interpret, and compare molecular properties with ease and precision.',
    src: four,
    link: 'https://images.unsplash.com/photo-1605106901227-991bd663255c?w=500&auto=format&fit=crop',
    color: '#000000',
    path:'Model4'
  },
  {
    title1: 'MolScan',
    title2:'Drug Filtering and Analysis',
    description:
      'This toolstreamlines bioactivity data analysis for drug discovery. It extracts, processes, and visualizes chemical and molecular data from ChEMBL, applying Lipinski\'s rules and pIC50 calculations. Equipped with advanced statistical analysis and Random Forest regression, it predicts bioactivity classes and generates insights into molecule descriptors, enabling effective drug targeting and precision research.',
    src: five,
    link: 'https://images.unsplash.com/photo-1605106715994-18d3fecffb98?w=500&auto=format&fit=crop&q=60',
    color: '#000000',
    path:'/Model2'
  }
];
export default function Models() {
  const container = useRef(null);
  const { scrollYProgress } = useScroll({
    target: container,
    offset: ['start start', 'end end'],
  });
  return (
    <ReactLenis root>
      <main className="bg-black" ref={container}>
        <section className="text-white w-full bg-black">
          {projects.map((project, i) => {
            const targetScale = 1 - (projects.length - i) * 0.05;
            return (
              <Card
                key={`p_${i}`}
                i={i}
                url={project?.link}
                src={project?.src}
                title1={project?.title1}
                title2={project?.title2}
                color={project?.color}
                description={project?.description}
                progress={scrollYProgress}
                range={[i * 0.25, 1]}
                targetScale={targetScale}
                path={project?.path}
              />
            );
          })}
        </section>
      </main>
    </ReactLenis>
  );
}
export const Card = ({
  i,
  title1,
  title2,
  description,
  src,
  url,
  color,
  progress,
  range,
  targetScale,
  path
}) => {
  const container = useRef(null);
  const { scrollYProgress } = useScroll({
    target: container,
    offset: ['start end', 'start start'],
  });
  const imageScale = useTransform(scrollYProgress, [0, 1], [2, 1]);
  const scale = useTransform(progress, range, [1, targetScale]);
  return (
    <div
      ref={container}
      className="h-screen flex items-center justify-center sticky top-0 "
    >
      <motion.div
        style={{
          backgroundColor: color,
          scale,
          top: `calc(-5vh + ${i * 25}px)`,
        }}
        className={`flex flex-col relative -top-[25%] h-[550px] w-[70%] dark:bg-black bg-white rounded-3xl p-4 shadow-xl border border-neutral-200 dark:border-white/[0.1]  shadow-black/[0.1] dark:shadow-white/[0.05] flex flex-col justify-between`}
      >
        
        <div className={`flex flex-cols mt-5 gap-10`}>
          <div className={`w-[40%] relative top-[10%] p-8 basis-3/5`}>
          <h2 className="mb-[30px] font-monument text-5xl text-center font-bold tracking-widest">{title1}</h2>
          <h2 className="mb-[10px] font-poppins text-2xl text-center font-bold tracking-wider">{title2}</h2>
            <p className="font-poppins text-base ml-[20px] mt-[20px] text-justify leading-relaxed">{description}</p>
            <span className="flex items-center gap-2 pt-2 mt-[30px] ml-[20px]">
                <ButtonHover2 link={path}/>             
            </span>
          </div>
          <div
            className={`relative basis-2/5 mt-[40px] mb-[20px] mr-[20px] overflow-hidden h-[400px]`}
          >
            <motion.div
              
            >
              <img  fill src={src} alt="image" className="object-cover rounded-[10px] h-[400px] w-full" />
            </motion.div>
          </div>
        </div>
      </motion.div>
    </div>
  );
};




const ButtonHover2 = (link) => {
  console.log(link.link)
  return (
    <>

<Link to={link.link} target="_blank" class="btn-shine">Get Started</Link>


<style>{` 

/* From Uiverse.io by neerajbaniwal */ 
.btn-shine {
letter-spacing:1.5px;
  position: relative;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  padding: 12px 48px;
  color: #fff;
  background: linear-gradient(to right, #9f9f9f 0, #fff 10%, #868686 20%);
  background-position: 0;
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  animation: shine 3s infinite linear;
  animation-fill-mode: forwards;
  -webkit-text-size-adjust: none;
  font-weight: 600;
  font-size: 18px;
  text-decoration: none;
  white-space: nowrap;
  font-family: "Poppins", sans-serif;
}
@-moz-keyframes shine {
  0% {
    background-position: 0;
  }
  60% {
    background-position: 180px;
  }
  100% {
    background-position: 180px;
  }
}
@-webkit-keyframes shine {
  0% {
    background-position: 0;
  }
  60% {
    background-position: 180px;
  }
  100% {
    background-position: 180px;
  }
}
@-o-keyframes shine {
  0% {
    background-position: 0;
  }
  60% {
    background-position: 180px;
  }
  100% {
    background-position: 180px;
  }
}
@keyframes shine {
  0% {
    background-position: 0;
  }
  60% {
    background-position: 180px;
  }
  100% {
    background-position: 180px;
  }
}


`}</style>
    </>
  );
};