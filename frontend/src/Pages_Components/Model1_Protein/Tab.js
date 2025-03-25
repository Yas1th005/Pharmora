'use client';
import { cn } from '../../components/lib/utils';
import { AnimatePresence, motion } from 'framer-motion';
import React, {
  createContext,
  useContext,
  useEffect,
  useState,
  isValidElement,
} from 'react';
const TabContext = createContext(undefined);
export const useTabs = () => {
  const context = useContext(TabContext);
  if (!context) {
    throw new Error('useTabs must be used within a TabsProvider');
  }
  return context;
};
export const TabsProvider = ({
  children,
  defaultValue,
  wobbly = true,
  hover = false,
}) => {
  const [activeTab, setActiveTab] = useState(defaultValue);
  const [prevIndex, setPrevIndex] = useState(0);
  const [tabsOrder, setTabsOrder] = useState([]);
  useEffect(() => {
    const order = [];
    children?.map((child) => {
      if (isValidElement(child)) {
        if (child.type === TabsContent) {
          order.push(child.props.value);
        }
      }
    });
    setTabsOrder(order);
  }, [children]);
  return (
    <TabContext.Provider
      value={{
        activeTab,
        setActiveTab,
        wobbly,
        hover,
        defaultValue,
        setPrevIndex,
        prevIndex,
        tabsOrder,
      }}
    >
      {children}
    </TabContext.Provider>
  );
};
export const TabsBtn = ({ children, className, value }) => {
  const {
    activeTab,
    setPrevIndex,
    setActiveTab,
    defaultValue,
    hover,
    wobbly,
    tabsOrder,
  } = useTabs();
  const handleClick = () => {
    setPrevIndex(tabsOrder.indexOf(activeTab));
    setActiveTab(value);
  };
  return (
    <>
      <>
        <motion.div
          className={cn(
            `cursor-pointer sm:p-2 p-1 sm:px-4 px-2 rounded-md relative `,
            className
          )}
          onFocus={() => {
            hover && handleClick();
          }}
          onMouseEnter={() => {
            hover && handleClick();
          }}
          onClick={handleClick}
        >
          {children}

          {activeTab === value && (
            <AnimatePresence mode="wait">
              <motion.div
                transition={{
                  layout: {
                    duration: 0.2,
                    ease: 'easeInOut',
                    delay: 0.2,
                  },
                }}
                layoutId={defaultValue}
                className="absolute w-full h-full left-0 top-0 bg-[#212121] rounded-md  z-[1]"
              />
            </AnimatePresence>
          )}

          {wobbly ? (
            <>
              {activeTab === value && (
                <AnimatePresence mode="wait">
                  <motion.div
                    transition={{
                      layout: {
                        duration: 0.4,
                        ease: 'easeInOut',
                        delay: 0.04,
                      },
                    }}
                    layoutId={defaultValue}
                    className="absolute w-full h-full left-0 top-0 bg-[#212121] rounded-md  z-[1] tab-shadow"
                  />
                </AnimatePresence>
              )}
              {activeTab === value && (
                <AnimatePresence mode="wait">
                  <motion.div
                    transition={{
                      layout: {
                        duration: 0.4,
                        ease: 'easeOut',
                        delay: 0.2,
                      },
                    }}
                    layoutId={`${defaultValue}b`}
                    className="absolute w-full h-full left-0 top-0 bg-[#212121] rounded-md  z-[1] tab-shadow"
                  />
                </AnimatePresence>
              )}
            </>
          ) : (
            <></>
          )}
        </motion.div>
      </>
    </>
  );
};
export const TabsContent = ({ children, className, value, yValue }) => {
  const { activeTab, tabsOrder, prevIndex } = useTabs();
  const isForward = tabsOrder.indexOf(activeTab) > prevIndex;
  return (
    <>
      <AnimatePresence mode="popLayout">
        {activeTab === value && (
          <motion.div
            initial={{ opacity: 0, y: yValue ? (isForward ? 10 : -10) : 0 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: yValue ? (isForward ? -50 : 50) : 0 }}
            transition={{
              duration: 0.3,
              ease: 'easeInOut',
              delay: 0.5,
            }}
            className={cn(' p-2 px-4 rounded-md relative', className)}
          >
            {activeTab === value ? children : null}
          </motion.div>
        )}
      </AnimatePresence>
    </>
  );
};
