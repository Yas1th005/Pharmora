import React from "react";
import { Header } from "./Header";

const HeaderImg = () => {
  return (
    <div
      className="relative h-screen bg-cover bg-center animate-moving-bg"
      style={{
        backgroundImage: `url('https://cdn.britannica.com/32/234732-050-0E5E77CA/DNA-strands-concept-illustration.jpg')`, // Replace with your image URL
      }}
    >
      {/* Gradient overlay to fade out the image at the bottom */}
      <div className="absolute inset-0 bg-gradient-to-b from-transparent via-transparent to-black"></div>

      {/* Content inside the background */}
      <div className="relative z-10 flex items-center justify-center h-full">
        <h1 className="text-4xl font-bold text-white">Your Content Here</h1>
      </div>
       <div className="absolute z-999 top-[30px] w-[100vw]"><Header/></div>
    </div>
  );
};

export default HeaderImg;


// https://scitechdaily.com/images/DNA-Genetics.gif