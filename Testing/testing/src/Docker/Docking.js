import React, { useState, useEffect } from 'react';
import './Docking.css';

function Docking() {
  const [receptorFile, setReceptorFile] = useState(null);
  const [ligandFile, setLigandFile] = useState(null);
  const [loading, setLoading] = useState(false);
  const [output, setOutput] = useState('');
  const [resultData, setResultData] = useState(null);
  const [resultFilename, setResultFilename] = useState(null);
  const [error, setError] = useState(null);
  const [serverStatus, setServerStatus] = useState('checking');

  // Check if the server is running on component mount
  useEffect(() => {
    checkServerStatus();
  }, []);

  const checkServerStatus = async () => {
    try {
      const response = await fetch('http://localhost:5000/api/health', {
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

  // Handle file selection
  const handleFileChange = (e, fileType) => {
    const file = e.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      if (fileType === 'receptor') {
        setReceptorFile({
          name: file.name,
          data: event.target.result
        });
      } else {
        setLigandFile({
          name: file.name,
          data: event.target.result
        });
      }
    };
    reader.readAsDataURL(file);
  };

  // Handle form submission
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
      const response = await fetch('http://localhost:5000/api/dock', {
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
      
      // Try to parse response as JSON
      let data;
      const contentType = response.headers.get('content-type');
      
      if (contentType && contentType.includes('application/json')) {
        data = await response.json();
      } else {
        // If not JSON, get the text and show it for debugging
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

  // Download result file
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
    <div className="App">
      <header className="App-header">
        <h1>Molecular Docking Web Application</h1>
        <p>Upload receptor and ligand files to perform molecular docking using AutoDock Vina</p>
        {serverStatus !== 'online' && (
          <div className={`server-status ${serverStatus}`}>
            Server Status: {serverStatus === 'checking' ? 'Checking connection...' : 
                           serverStatus === 'offline' ? 'Offline - Please start the Flask server' :
                           'Error connecting to server'}
          </div>
        )}
      </header>
      
      <main className="App-main">
        <div className="container">
          <div className="file-upload-container">
            <h2>Upload Files</h2>
            <form onSubmit={handleSubmit}>
              <div className="form-group">
                <label htmlFor="receptor-file">Receptor File (PDBQT):</label>
                <input
                  type="file"
                  id="receptor-file"
                  accept=".pdbqt"
                  onChange={(e) => handleFileChange(e, 'receptor')}
                  disabled={loading || serverStatus !== 'online'}
                />
                {receptorFile && <p className="file-info">Selected: {receptorFile.name}</p>}
              </div>
              
              <div className="form-group">
                <label htmlFor="ligand-file">Ligand File (PDBQT):</label>
                <input
                  type="file"
                  id="ligand-file"
                  accept=".pdbqt"
                  onChange={(e) => handleFileChange(e, 'ligand')}
                  disabled={loading || serverStatus !== 'online'}
                />
                {ligandFile && <p className="file-info">Selected: {ligandFile.name}</p>}
              </div>
              
              <button 
                type="submit" 
                className="dock-button"
                disabled={loading || !receptorFile || !ligandFile || serverStatus !== 'online'}
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
            <div className="terminal-container">
              <h2>Terminal Output</h2>
              <div className="terminal">
                {output ? output.split('\n').map((line, i) => (
                  <div key={i} className="terminal-line">{line}</div>
                )) : (
                  <div className="terminal-placeholder">Terminal output will appear here</div>
                )}
              </div>
            </div>
            
            {resultData && (
              <div className="visualization-container">
                <h2>Docking Results</h2>
                <div className="result-actions">
                  <button onClick={handleDownload} className="download-button">
                    Download Results
                  </button>
                </div>
                
                <div className="result-preview">
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
            )}
          </div>
        </div>
      </main>
    </div>
  );
}

export default Docking;