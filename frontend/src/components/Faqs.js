'use client';
import React, { useState } from 'react';
import { AnimatePresence, motion } from 'framer-motion';
import { Plus } from 'lucide-react';
const tabs = [
  {
    title: 'What problem does this project solve?',
    description:
      'This project addresses the high cost, long timelines, and inefficiencies in drug discovery by integrating AI-driven solutions to make the process faster, more affordable, and scalable for biotech firms and researchers.',
    imageUrl:
      'https://images.unsplash.com/photo-1709949908058-a08659bfa922?q=80&w=1200&auto=format',
  },
  {
    title: 'How does the platform speed up drug discovery?',
    description:
      'The platform leverages Machine Learning (ML), Deep Learning, and Generative AI (such as fine-tuned Llama 3.2 and GANs) to streamline drug discovery by optimizing data processing, reducing trial-and-error cycles, and accelerating research.',
    imageUrl:
      'https://images.unsplash.com/photo-1548192746-dd526f154ed9?q=80&w=1200&auto=format',
  },
  {
    title: 'Who can benefit from this solution?',
    description:
      'Small biotech firms, academic researchers, and pharmaceutical companies looking for a cost-effective and efficient way to discover new drugs can greatly benefit from this platform.',
    imageUrl:
      'https://images.unsplash.com/photo-1693581176773-a5f2362209e6?q=80&w=1200&auto=format',
  },
  {
    title: 'What technologies power this project?',
    description:
      'The project is built using Next.js 14, React, Streamlit, Flask, Docker, PostgreSQL, Scikit-Learn, TensorFlow, and Generative AI techniques, ensuring high performance and scalability.',
    imageUrl:
      'https://images.unsplash.com/photo-1548192746-dd526f154ed9?q=80&w=1200&auto=format',
  },
  {
    title: 'Who are the developers behind this project?',
    description:
      'The project is developed by a team of second-year B.Tech CSE students from VIT Chennai, specializing in backend development, frontend development, and machine learning.',
    imageUrl:
      'https://images.unsplash.com/photo-1548192746-dd526f154ed9?q=80&w=1200&auto=format',
  },
];
function Faq() {
  const [activeIndex, setActiveIndex] = useState(0);
  const [activeItem, setActiveItem] = useState(tabs[0]);
  const handleClick = async (index) => {
    setActiveIndex(activeIndex === index ? null : index);
    const newActiveItem = tabs.find((_, i) => i === index);
    setActiveItem(newActiveItem);
  };
  return (
    <>
      <div className="container mx-auto pb-10 pt-2 w-[80vw] font-poppins">
        {/* <h1 className="uppercase text-center text-4xl text-white font-bold pt-2 pb-4">
          FAQ
        </h1> */}
        <div className="h-fit border  rounded-[20px] p-2 bg-transparent">
          {tabs.map((tab, index) => (
            <motion.div
              key={index}
              className={`overflow-hidden ${
                index !== tabs.length - 1 ? 'border-b' : ''
              }`}
              onClick={() => handleClick(index)}>
              <button
                className={`p-3 px-2 w-full cursor-pointer sm:text-base text-xs items-center transition-all font-semibold dark:text-white text-white   flex gap-2 
               `}>
                <Plus
                  className={`${
                    activeIndex === index ? 'rotate-45' : 'rotate-0 '
                  } transition-transform ease-in-out w-5 h-5  text-white`}
                />
                {tab.title}
              </button>
              <AnimatePresence mode="sync">
                {activeIndex === index && (
                  <motion.div
                    initial={{ height: 0, opacity: 0 }}
                    animate={{ height: 'auto', opacity: 1 }}
                    exit={{ height: 0, opacity: 0 }}
                    transition={{
                      duration: 0.3,
                      ease: 'easeInOut',
                      delay: 0.14,
                    }}>
                    <p
                      className={`dark:text-white text-white p-3 xl:text-base sm:text-sm text-xs pt-0 w-[90%]`}>
                      {tab.description}
                    </p>
                  </motion.div>
                )}
              </AnimatePresence>
            </motion.div>
          ))}
        </div>
      </div>
    </>
  );
}
export default Faq;
