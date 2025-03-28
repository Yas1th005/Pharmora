import { cn } from '../../components/lib/utils';
import {
  createContext,
  forwardRef,
  useCallback,
  useContext,
  useEffect,
  useRef,
  useState,
} from 'react';
import { useDropzone } from 'react-dropzone';
import { toast } from 'sonner';
import { Trash2 as RemoveIcon } from 'lucide-react';
const FileUploaderContext = createContext(null);
export const useFileUpload = () => {
  const context = useContext(FileUploaderContext);
  if (!context) {
    throw new Error('useFileUpload must be used within a FileUploaderProvider');
  }
  return context;
};
export const FileUploader = forwardRef(
  (
    {
      className,
      dropzoneOptions,
      value,
      onValueChange,
      reSelect,
      orientation = 'vertical',
      children,
      dir,
      ...props
    },
    ref
  ) => {
    const [isFileTooBig, setIsFileTooBig] = useState(false);
    const [isLOF, setIsLOF] = useState(false);
    const [activeIndex, setActiveIndex] = useState(-1);
    const {
      accept = {
        'image/*': ['.jpg', '.jpeg', '.png', '.gif'],
        'video/*': ['.mp4', '.MOV', '.AVI'],
      },
      maxFiles = 1,
      maxSize = 4 * 1024 * 1024,
      multiple = true,
    } = dropzoneOptions;
    const reSelectAll = maxFiles === 1 ? true : reSelect;
    const direction = dir === 'rtl' ? 'rtl' : 'ltr';
    const removeFileFromSet = useCallback(
      (i) => {
        if (!value) return;
        const newFiles = value.filter((_, index) => index !== i);
        onValueChange(newFiles);
      },
      [value, onValueChange]
    );
    const handleKeyDown = useCallback(
      (e) => {
        e.preventDefault();
        e.stopPropagation();
        if (!value) return;
        const moveNext = () => {
          const nextIndex = activeIndex + 1;
          setActiveIndex(nextIndex > value.length - 1 ? 0 : nextIndex);
        };
        const movePrev = () => {
          const nextIndex = activeIndex - 1;
          setActiveIndex(nextIndex < 0 ? value.length - 1 : nextIndex);
        };
        const prevKey =
          orientation === 'horizontal'
            ? direction === 'ltr'
              ? 'ArrowLeft'
              : 'ArrowRight'
            : 'ArrowUp';
        const nextKey =
          orientation === 'horizontal'
            ? direction === 'ltr'
              ? 'ArrowRight'
              : 'ArrowLeft'
            : 'ArrowDown';
        if (e.key === nextKey) {
          moveNext();
        } else if (e.key === prevKey) {
          movePrev();
        } else if (e.key === 'Enter' || e.key === 'Space') {
          if (activeIndex === -1) {
            dropzoneState.inputRef.current?.click();
          }
        } else if (e.key === 'Delete' || e.key === 'Backspace') {
          if (activeIndex !== -1) {
            removeFileFromSet(activeIndex);
            if (value.length - 1 === 0) {
              setActiveIndex(-1);
              return;
            }
            movePrev();
          }
        } else if (e.key === 'Escape') {
          setActiveIndex(-1);
        }
      },
      [value, activeIndex, removeFileFromSet]
    );
    const onDrop = useCallback(
      (acceptedFiles, rejectedFiles) => {
        const files = acceptedFiles;
        if (!files) {
          toast.error('file error , probably too big');
          return;
        }
        const newValues = value ? [...value] : [];
        if (reSelectAll) {
          newValues.splice(0, newValues.length);
        }
        files.forEach((file) => {
          if (newValues.length < maxFiles) {
            newValues.push(file);
          }
        });
        onValueChange(newValues);
        if (rejectedFiles.length > 0) {
          for (let i = 0; i < rejectedFiles.length; i++) {
            if (rejectedFiles[i].errors[0]?.code === 'file-too-large') {
              toast.error(
                `File is too large. Max size is ${maxSize / 1024 / 1024}MB`
              );
              break;
            }
            if (rejectedFiles[i].errors[0]?.message) {
              toast.error(rejectedFiles[i].errors[0].message);
              break;
            }
          }
        }
      },
      [reSelectAll, value]
    );
    useEffect(() => {
      if (!value) return;
      if (value.length === maxFiles) {
        setIsLOF(true);
        return;
      }
      setIsLOF(false);
    }, [value, maxFiles]);
    const opts = dropzoneOptions
      ? dropzoneOptions
      : { accept, maxFiles, maxSize, multiple };
    const dropzoneState = useDropzone({
      ...opts,
      onDrop,
      onDropRejected: () => setIsFileTooBig(true),
      onDropAccepted: () => setIsFileTooBig(false),
    });
    return (
      <FileUploaderContext.Provider
        value={{
          dropzoneState,
          isLOF,
          isFileTooBig,
          removeFileFromSet,
          activeIndex,
          setActiveIndex,
          orientation,
          direction,
        }}
      >
        <div
          ref={ref}
          tabIndex={0}
          onKeyDownCapture={handleKeyDown}
          className={cn(
            'grid w-full focus:outline-none overflow-hidden ',
            className,
            {
              'gap-2': value && value.length > 0,
            }
          )}
          dir={dir}
          {...props}
        >
          {children}
        </div>
      </FileUploaderContext.Provider>
    );
  }
);
FileUploader.displayName = 'FileUploader';
export const FileUploaderContent = forwardRef(
  ({ children, className, ...props }, ref) => {
    const { orientation } = useFileUpload();
    const containerRef = useRef(null);
    return (
      <div
        className={cn('w-full px-1')}
        ref={containerRef}
        aria-description="content file holder"
      >
        <div
          {...props}
          ref={ref}
          className={cn(
            ' rounded-xl gap-1',
            orientation === 'horizontal' ? 'grid grid-cols-2' : 'flex flex-col',
            className
          )}
        >
          {children}
        </div>
      </div>
    );
  }
);
FileUploaderContent.displayName = 'FileUploaderContent';
export const FileUploaderItem = forwardRef(
  ({ className, index, children, ...props }, ref) => {
    const { removeFileFromSet, activeIndex, direction } = useFileUpload();
    const isSelected = index === activeIndex;
    return (
      <div
        ref={ref}
        className={cn(
          'h-7 p-1 border rounded-md justify-between overflow-hidden  w-full cursor-pointer relative hover:bg-primary-foreground',
          className,
          isSelected ? 'bg-muted' : ''
        )}
        {...props}
      >
        <div className="font-medium   leading-none tracking-tight flex items-center gap-1.5 h-full w-full">
          {children}
        </div>
        <button
          type="button"
          className={cn(
            'absolute bg-primary rounded text-background p-1',
            direction === 'rtl' ? 'top-1 left-1' : 'top-[0.145em] right-1'
          )}
          onClick={() => removeFileFromSet(index)}
        >
          <span className="sr-only">remove item {index}</span>
          <RemoveIcon className="w-3 h-3 hover:stroke-destructive  duration-200 ease-in-out" />
        </button>
      </div>
    );
  }
);
FileUploaderItem.displayName = 'FileUploaderItem';
export const FileInput = forwardRef(
  ({ className, parentclass, dropmsg, children, ...props }, ref) => {
    const { dropzoneState, isFileTooBig, isLOF } = useFileUpload();
    const rootProps = isLOF ? {} : dropzoneState.getRootProps();
    return (
      <div
        ref={ref}
        {...props}
        className={cn(
          'relative w-full',
          parentclass,
          isLOF ? 'opacity-50 cursor-not-allowed' : 'cursor-pointer'
        )}
      >
        <div
          className={cn(
            'w-full rounded-lg transition-colors duration-300 ease-in-out',
            dropzoneState.isDragAccept && 'border-green-500 bg-green-50',
            dropzoneState.isDragReject && 'border-red-500 bg-red-50',
            isFileTooBig && 'border-red-500 bg-red-200',
            !dropzoneState.isDragActive &&
              'border-gray-300 hover:border-gray-400',
            className
          )}
          {...rootProps}
        >
          {children}
          {dropzoneState.isDragActive && (
            <div className="absolute inset-0 flex items-center justify-center bg-primary-foreground/60 backdrop-blur-sm rounded-lg">
              <p className="text-primary font-medium">{dropmsg}</p>
            </div>
          )}
        </div>
        <input
          ref={dropzoneState.inputRef}
          disabled={isLOF}
          {...dropzoneState.getInputProps()}
          className={cn(isLOF && 'cursor-not-allowed')}
        />
      </div>
    );
  }
);
