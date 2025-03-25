import React from "react";




import twitter from "../imgaes/twitter.png"
import insta from "../imgaes/instagram.png"
import github from "../imgaes/github.png"


const Footer = () => {
  return (
    <div className="mt-4 bg-background">
      <div className="flex flex-col w-full h-fit bg-background text-[#ffffff] px-14 py-14">
        <div className="w-full border-t border-gray-500 my-8"></div>
        <div className="flex flex-row px-10">
          <div className="flex flex-col gap-2 justify-center w-[35%] w-[35%]">
            <div className="flex items-center w-full gap-4">
              <div className="text-5xl font-bold font-custom tracking-wider">
                Pharmora
              </div>
            </div>
            <div className="flex mt-2 justify-start grid-cols-3 gap-[56px]">
              <div className="flex flex-cols justify-center align-center w-[35px] gap-[60px] ml-[95px]">
              <img src={insta}/>
              <img src={github}/>
              <img src={twitter}/>
              </div>
              
            </div>
          </div>
          <div className="flex flex-row w-[65%] justify-end gap-16 text-nowrap">
            <div className="grid grid-cols-2 gap-16">
              <div className="flex flex-col gap-2">
                <div className="font-bold uppercase text-[#9ca3af] pb-3">
                  Comany
                </div>{" "}
                <a href="#xxx" className="hover:underline">
                  About Us
                </a>{" "}
                <a href="#xxx" className="hover:underline">
                  Contact
                </a>{" "}
                <a href="#xxx" className="hover:underline">
                  Support
                </a>{" "}
              </div>
              <div className="flex flex-col gap-2">
                <div className="font-bold uppercase text-[#9ca3af] pb-3">
                  Legal
                </div>{" "}
                <a href="#xxx" className="hover:underline">
                  Imprint
                </a>{" "}
                <a href="#xxx" className="hover:underline">
                  Privacy Policy
                </a>{" "}
                <a href="#xxx" className="hover:underline">
                  Terms of Use
                </a>
              </div>
            </div>
          </div>
        </div>
        <div className="w-full border-t border-gray-500 my-8"></div>
        <div className="text-center">
          © 2024 Pharmora - All rights reserved.
        </div>
      </div>
    </div>
  );
};

export default Footer;
