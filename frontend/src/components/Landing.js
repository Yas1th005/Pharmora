// import React from 'react';
// import { SparklesGblobe } from './ui/sparkles';
// import { Header } from './Header';
// import LandModel from './LandModel';
// function Landing() {
//   return (
//     <>
//       <main className=" h-[100vh] w-full  overflow-hidden bg-black text-white ">
//         <section className="container mx-auto relative h-[100vh] mt-4 border border-white/10 w-full overflow-hidden rounded-2xl">
//         <Header/>
//         <div className='flex flex-cols justfiy-center align-center'>
//             <article className="grid gap-4 text-center relative z-10 pt-10">
//                 <span className="inline-block xl:text-base text-sm border p-1 px-3 w-fit mx-auto rounded-full border-[#3273ff] bg-[#0f1c35]">
//                 Early Access
//                 </span>
//                 <h1 className="2xl:text-6xl  xl:text-5xl text-5xl font-semibold bg-gradient-to-b from-[#edeffd] to-[#7b9cda] bg-clip-text text-transparent leading-[100%] tracking-tighter">
//                 Become an Animation Expert <br /> Easily at Our Academy
//                 </h1>
//                 <span>
//                 Our expert-led courses are perfect for all skill levels. Gain{' '}
//                 <br />
//                 hands-on experience and create stunning animations <br />{' '}
//                 effortlessly. Join us today!
//                 </span>
//                 <button className="border border-blue-400 w-fit p-2 px-4 rounded-md bg-blue-900/40 hover:bg-blue-900/60  backdrop-blur-2xl mx-auto text-white">
//                 Take The Course
//                 </button>
//             </article>
//            <div className='z-[3]'><LandModel/></div> 
//         </div>
          

//           <div className="absolute bottom-0 z-[2] h-[400px] w-screen overflow-hidden [mask-image:radial-gradient(100%_50%,white,transparent)] before:absolute before:inset-0 before:bg-[radial-gradient(circle_at_bottom_center,#3273ff,transparent_90%)] before:opacity-40 after:absolute">
//             <SparklesGblobe
//               density={1800}
//               speed={1.2}
//               color="#48b6ff"
//               direction="top"
//               className="absolute inset-x-0 bottom-0 h-full w-full "
//             />
//           </div>
//         </section>
//       </main>
//     </>
//   );
// }
// export default Landing;



// import React from 'react';
// import { SparklesGblobe } from '../components/ui/sparkles';
// import Earth from '../components/ui/globe';
// import LandModel from './LandModel';
// function Landing() {
//   return (
//     <>
//       <div className="h-[100vh] bg-black text-white">
//         <article className="grid gap-4 text-center relative z-10 pt-10">
//           <span className="inline-block text-sm border p-1 px-3 w-fit mx-auto rounded-full border-[#3273ff] bg-[#0f1c35]">
//             Get Access
//           </span>
//           <h1 className="text-4xl  font-semibold bg-gradient-to-b from-[#edeffd] to-[#7b9cda] bg-clip-text text-transparent leading-[100%] tracking-tighter">
//             Design with a Global
//             <br />
//             Perspective, Innovate with Ease.
//           </h1>
//           <div className='mx-auto'><LandModel/></div>
//         </article>

//         <div className="relative -mt-32 h-80 w-screen overflow-hidden [mask-image:radial-gradient(50%_50%,white,transparent)] before:absolute before:inset-0 before:bg-[radial-gradient(circle_at_bottom_center,#3273ff,transparent_90%)] before:opacity-40 after:absolute after:-left-1/2 after:top-1/2 after:aspect-[1/0.7] after:w-[200%] after:rounded-[10%] after:border-t after:border-[#163474] after:bg-[#08132b]">
//           <SparklesGblobe
//             density={1200}
//             className="absolute inset-x-0 bottom-0 h-full w-full "
//           />
//         </div>
//       </div>
//     </>
//   );
// }
// export default Landing;




"use client";
import React from "react";
import { BackgroundBeams } from "../components/ui/background_beams";
import { Header } from "./Header";
import LandModel from "./LandModel";

export default function Landing() {
  return (
    (<div
      className="h-[100vh] w-full rounded-md bg-black relative flex flex-col items-center justify-center antialiased">
      <div className=""><Header/></div>
      <div className="flex flex-cols gap-[70px] w-[80vw] justify-center align-center">
        <div className="ml-[100px] p-4 mt-[180px]">
        <p
          className={
            ` z-50 h-32 font-custom tracking-wider inline-block whitespace-pre-line text-transparent pointer-events-none bg-gradient-to-b from-neutral-900 to-white bg-clip-text text-center text-9xl font-semibold leading-none `
          }
        >
          Pharmora
        </p>
        <p className="text-transparent inline-block pointer-events-none bg-gradient-to-b from-neutral-600 to-white bg-clip-text z-50 text-center font-semibold leading-none whitespace-pre-line text-[32px] mx-auto">
          Accelerating Innovation in Bioinformatics
        </p>
        </div>
        <div className="z-10"><LandModel/></div>
      </div>
      
      <BackgroundBeams />
    </div>)
  );
}



