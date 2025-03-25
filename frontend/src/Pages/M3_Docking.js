import React, { useState, useEffect } from 'react';
import { FileUploader, FileUploaderContent, FileUploaderItem, FileInput } from '../Pages_Components/Model3_Docking/file-upload';
import { Paperclip } from 'lucide-react';

const FileSvgDraw = ({ name }) => {
  return (
    <>
      <svg
        className="w-8 h-8 mb-3 text-primary"
        aria-hidden="true"
        xmlns="http://www.w3.org/2000/svg"
        fill="none"
        viewBox="0 0 20 16">
        <path
          stroke="currentColor"
          strokeLinecap="round"
          strokeLinejoin="round"
          strokeWidth="2"
          d="M13 13h3a3 3 0 0 0 0-6h-.025A5.56 5.56 0 0 0 16 6.5 5.5 5.5 0 0 0 5.207 5.021C5.137 5.017 5.071 5 5 5a4 4 0 0 0 0 8h2.167M10 15V6m0 0L8 8m2-2 2 2"
        />
      </svg>
      <p className="mb-1 text-sm text-primary">
        <span className="font-semibold">Click to upload</span>
        &nbsp; or drag and drop
      </p>
      <p className="text-base text-primary font-poppins">{name}</p>
    </>
  );
};

const FileUploaderTest = ({ name, onFileChange }) => {
  const [files, setFiles] = useState([]);

  const dropZoneConfig = {
    maxFiles: 1,
    maxSize: 1024 * 1024 * 4,
    multiple: false,
  };

  const handleFileChange = (files) => {
    setFiles(files);
    if (files.length > 0) {
      const file = files[0];
      const reader = new FileReader();
      reader.onload = (event) => {
        onFileChange({
          name: file.name,
          data: event.target.result
        });
      };
      reader.readAsDataURL(file);
    }
  };

  return (
    <FileUploader
      value={files}
      onValueChange={handleFileChange}
      dropzoneOptions={dropZoneConfig}
      className="relative rounded-lg p-2 w-96 mx-auto">
      <FileInput className="outline-dashed outline-2 bg-background outline-primary/40">
        <div className="flex items-center justify-center flex-col pt-3 pb-4 w-full">
          <FileSvgDraw name={name} />
        </div>
      </FileInput>
      <FileUploaderContent>
        {files &&
          files.length > 0 &&
          files.map((file, i) => (
            <FileUploaderItem key={i} index={i} className="bg-background">
              <Paperclip className="h-4 w-4 flex-shrink-0 stroke-current" />
              <p className="text-ellipsis inline-block overflow-hidden text-xs w-full">
                {file.name}
              </p>
            </FileUploaderItem>
          ))}
      </FileUploaderContent>
    </FileUploader>
  );
};

function Docking() {
  const [receptorFile, setReceptorFile] = useState(null);
  const [ligandFile, setLigandFile] = useState(null);
  const [loading, setLoading] = useState(false);
  const [output, setOutput] = useState('');
  const [resultData, setResultData] = useState(null);
  const [resultFilename, setResultFilename] = useState(null);
  const [error, setError] = useState(null);
  const [serverStatus, setServerStatus] = useState('checking');

  useEffect(() => {
    checkServerStatus();
  }, []);

  const checkServerStatus = async () => {
    try {
      const response = await fetch('http://localhost:5002/api/health', {
        method: 'GET',
        headers: {
          'Accept': 'application/json'
        }
      });
      
      if (response.ok) {
        setServerStatus('online');
      } else {
        setServerStatus('error');
      }
    } catch (err) {
      console.error('Server connection error:', err);
      setServerStatus('offline');
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    
    if (!receptorFile || !ligandFile) {
      alert('Please select both receptor and ligand files');
      return;
    }
    
    setLoading(true);
    setOutput('Starting docking process...\n');
    setError(null);
    setResultData(null);
    
    try {
      const response = await fetch('http://localhost:5002/api/dock', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Accept': 'application/json'
        },
        body: JSON.stringify({
          receptorFile: receptorFile,
          ligandFile: ligandFile
        }),
      });
      
      let data;
      const contentType = response.headers.get('content-type');
      
      if (contentType && contentType.includes('application/json')) {
        data = await response.json();
      } else {
        const textResponse = await response.text();
        throw new Error(`Server returned non-JSON response: ${textResponse.substring(0, 100)}...`);
      }
      
      if (response.ok && data.success) {
        setOutput(data.output);
        if (data.resultData) {
          setResultData(data.resultData);
          setResultFilename(data.resultFile);
        }
      } else {
        setError(data.error || 'An error occurred during the docking process');
        setOutput(data.output || '');
      }
    } catch (err) {
      console.error('Request failed:', err);
      setError(`Failed to connect to the server: ${err.message}`);
    } finally {
      setLoading(false);
    }
  };

  const handleDownload = () => {
    if (!resultData || !resultFilename) return;
    
    try {
      const byteCharacters = atob(resultData);
      const byteNumbers = new Array(byteCharacters.length);
      for (let i = 0; i < byteCharacters.length; i++) {
        byteNumbers[i] = byteCharacters.charCodeAt(i);
      }
      const byteArray = new Uint8Array(byteNumbers);
      const blob = new Blob([byteArray], {type: 'application/octet-stream'});
      
      const link = document.createElement('a');
      link.href = URL.createObjectURL(blob);
      link.download = resultFilename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    } catch (err) {
      setError(`Failed to download file: ${err.message}`);
    }
  };

  return (
    <div className="bg-black p-4 text-white h-screen">
      <header className="w-[80vw] text-center font-poppins mx-auto mt-3">
        <h1 className='text-4xl font-bold tracking-wider'>DockFlow</h1>
        <p className='text-base mt-3'>Upload receptor and ligand files to perform molecular docking using AutoDock Vina</p>
      </header>
      
      <main className="App-main mt-6">
        <div className="container">
          <div className='flex justify-center items-center w-[60vw] mx-auto'>
            <FileUploaderTest name="Receptor File" onFileChange={setReceptorFile} />
            <FileUploaderTest name="Ligand File" onFileChange={setLigandFile} />
          </div>
          
          <div className="file-upload-container w-[80vw]">
            <form onSubmit={handleSubmit}>
              <button
                type="submit"
                disabled={loading || !receptorFile || !ligandFile || serverStatus !== 'online'}
                className="dock-button w-[30vw] py-2 font-semibold rounded-[10px] bg-violet-400 text-gray-900 hover:cursor-pointer ml-[35vw] mt-7"
              >
                {loading ? 'Docking in Progress...' : 'Run Docking'}
              </button>
            </form>
          </div>
          
          {error && (
            <div className="error-message">
              <h3>Error:</h3>
              <p>{error}</p>
            </div>
          )}
          
          <div className="results-container">
            {/* <div className="terminal-container">
              <h2>Terminal Output</h2>
              <div className="terminal">
                {output ? output.split('\n').map((line, i) => (
                  <div key={i} className="terminal-line">{line}</div>
                )) : (
                  <div className="terminal-placeholder">Terminal output will appear here</div>
                )}
              </div>
            </div> */}
            
            {resultData && (
              <div className='bg-black w-[100vw] mt-6'>
                <div className="w-[80vw] mx-auto h-fit bg-black">
                <h2 className='text-center text-2xl font-poppins mb-2'>Docking Results</h2>
                <div className="result-actions">
                  <button
                    onClick={handleDownload}
                    className="w-full py-2 font-semibold rounded-[10px] bg-violet-400 text-gray-900 mb-5"
                  >
                    Download Results
                  </button>
                </div>
                
                <div className="border border-white p-4">
                  <h3>Result File Preview:</h3>
                  <pre className="pdbqt-content">
                    {atob(resultData).substring(0, 5000) + 
                      (atob(resultData).length > 5000 ? '...(truncated)' : '')}
                  </pre>
                </div>
                
                <p className="result-note">
                  Note: For visualization, consider importing this PDBQT file into a molecular viewer like PyMOL or Chimera.
                </p>
              </div>
              </div>
              
            )}
          </div>
        </div>
      </main>
    </div>
  );
}

export default Docking;