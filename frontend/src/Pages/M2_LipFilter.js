import {
  Dialog,
  DialogTrigger,
  DialogContent,
  DialogTitle,
  DialogImage,
  DialogClose,
  DialogDescription,
  DialogContainer,
} from '../Pages_Components/Model2_LipFilter/TypeSelect';
import { Plus } from 'lucide-react';
import MolecularPropertyFilter from "../Pages_Components/Model2_LipFilter/EachFilter/Molecular";
import DescriptorDistribution from "../Pages_Components/Model2_LipFilter/EachFilter/Descriptor";
import ChemicalSpace from "../Pages_Components/Model2_LipFilter/EachFilter/Chemical";
import drug from "../imgaes/DrugVista/DrugVista.webp"

import one from "../imgaes/molscan/1.png";
import two from "../imgaes/molscan/2.png";
import three from "../imgaes/molscan/3.png";
import main from "../imgaes/models/5.jpg"

const items = [
  {
    id: 1,
    url: one,
    title: 'Molecular Property Filter',
    description:
      'Immerse yourself in our cutting-edge interactive gallery, designed to showcase a diverse array of visual content with unparalleled clarity and style. This feature allows users to effortlessly navigate through high-resolution images, from awe-inspiring landscapes to intimate portraits and abstract art. With smooth transitions, intuitive controls, and responsive design, our gallery adapts to any device, ensuring a seamless browsing experience. Dive deeper into each piece with expandable information panels, offering insights into the artist, technique, and story behind each image. ',
    tags: ['Sunrise', 'Mountains', 'Golden', 'Scenic', 'Inspiring'],
    activeTab:"filter",
  },
  {
    id: 2,
    url: two,
    title: 'Descriptor Distribution',
    description: `Embark on a virtual journey around the world with our state-of-the-art 3D globe feature. This interactive marvel allows users to explore geographical data, global trends, and worldwide connections with unprecedented ease and detail. Spin the globe with a flick of your mouse, zoom into street-level views, or soar high for a continental perspective. Our globe section integrates real-time data feeds, showcasing everything from climate patterns and population densities to economic indicators and cultural hotspots. Customizable layers let you focus on specific data sets, while intuitive tooltips provide in-depth information at every turn. `,
    tags: ['Misty', 'Path', 'Mysterious', 'Serene', 'Rugged'],
    activeTab:"distribution",
  },
  {
    id: 3,
    url: three,
    title: 'Chemical Space',
    description: `Transform your browsing experience with our mesmerizing Image Mouse Trail feature. As you move your cursor across the screen, watch in wonder as a trail of carefully curated images follows in its wake, creating a dynamic and engaging visual spectacle. This innovative feature goes beyond mere aesthetics; it's an interactive showcase of your content, products, or artwork. Each image in the trail can be clickable, leading to detailed views or related content, turning casual mouse movements into opportunities for discovery.`,
    tags: ['Pathway', 'Adventure', 'Peaks', 'Challenging', 'Breathtaking'],
    activeTab:"space",
  },
];
export default function LinearCard() {
  return (
    <div className='bg-black'>
        <div className="p-5 mx-auto sm:p-10 md:p-16 dark:bg-gray-100 dark:text-gray-800">
            <div className="flex flex-col w-[70vw] mx-auto overflow-hidden rounded items-center">
                <img src={main} alt="" className="w-full h-60 sm:h-96 dark:bg-gray-500 object-cover" />
                <div className="p-6 pb-12 m-4 mx-auto -mt-16 space-y-6 lg:max-w-2xl sm:px-10 sm:mx-12 lg:rounded-md bg-white w-[50vw]">
                    <div className="space-y-2">
                        <a rel="noopener noreferrer" href="#" className="inline-block text-2xl sm:text-3xl text-center w-full font-monument text-gray-800">Drug Discovery Toolkit 🔬</a>
                        <p className="text-lg dark:text-gray-600 text-center uppercase font-monument tracking-[2px]">PROSPECTRA</p>
                    </div>
                    <div className="dark:text-gray-800 text-center">
                        <p>ProSpectra is an advanced molecular analysis platform integrating cheminformatics and visualization for precise research.</p>
                    </div>
                </div>
            </div>
        </div>

        <div className="w-[90vw] flex gap-4 mx-auto pb-[100px]">
            {items.map((item, i) => {
                return (
                <>
                    <Dialog
                    transition={{
                        type: 'spring',
                        bounce: 0.05,
                        duration: 0.5,
                    }}>
                    <DialogTrigger
                        style={{
                        borderRadius: '12px',
                        }}
                        className="flex w-full h-64 flex-col overflow-hidden  border bg-black hover:bg-gray-950">
                        <DialogImage
                        src={item.url}
                        alt=""
                        className="h-32 w-32 object-fill mx-auto mt-8"
                        />
                        <div className="flex flex-grow flex-row items-end justify-between p-3">
                        <div>
                            <DialogTitle className="text-white text-base font-monument tracking-widest">
                            {item.title}
                            </DialogTitle>
                        </div>
                        <button className="absolute bottom-2 right-2 p-2 dark:bg-gray-900 bg-gray-400 hover:bg-gray-500 rounded-full dark:hover:bg-gray-800">
                            <Plus className="w-5 h-5" />
                        </button>
                        </div>
                    </DialogTrigger>
                    <DialogContainer className="">
                        <DialogContent
                        style={{
                            borderRadius: '24px',
                        }}
                        className="relative flex h-[100vh] mx-auto flex-col overflow-y-auto bg-black w-[100vw]">
                        {/* <DialogImage
                            src={item.url.src}
                            alt=""
                            className="h-full object-contain w-[60%] mx-auto"
                        /> */}
                        <div className="p-6">
                            {/* <DialogTitle className="text-3xl text-white text-center font-monument">
                            {item.title}
                            </DialogTitle>

                            <DialogDescription
                            disableLayoutAnimation
                            variants={{
                                initial: { opacity: 0, scale: 0.8, y: -40 },
                                animate: { opacity: 1, scale: 1, y: 0 },
                                exit: { opacity: 0, scale: 0.8, y: -50 },
                            }}>
                            <p className="mt-2 text-zinc-500 dark:text-zinc-500">
                                {item.description}
                            </p>
                            </DialogDescription> */}
                            {item.activeTab === 'filter' && <MolecularPropertyFilter title={item.title} desc={item.description}/>}
                            {item.activeTab === 'distribution' && <DescriptorDistribution />}
                            {item.activeTab === 'space' && <ChemicalSpace />}
                        </div>
                        <DialogClose className="text-zinc-50 dark:bg-gray-900 bg-transparent border border-white p-2 hover:bg-gray-500 rounded-full dark:hover:bg-gray-800" />
                        </DialogContent>
                    </DialogContainer>
                    </Dialog>
                </>
                );
            })}
            </div>
    </div>
    
  );
}
