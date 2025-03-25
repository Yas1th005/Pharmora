from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import os
import subprocess
import base64
import tempfile
from pathlib import Path, PureWindowsPath
import uuid
import traceback
import platform

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Create necessary directories with WSL-compatible paths
def ensure_directories():
    """Create necessary directories and ensure WSL compatibility"""
    directories = ['receptor', 'ligand', 'results']
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        # Ensure proper permissions in WSL
        if is_wsl():
            subprocess.run(['chmod', '777', directory], check=True)

def is_wsl():
    """Check if running in WSL"""
    try:
        with open('/proc/version', 'r') as f:
            return 'microsoft' in f.read().lower()
    except:
        return False

def convert_to_wsl_path(windows_path):
    """Convert Windows path to WSL path"""
    # Convert backslashes to forward slashes
    path = str(windows_path).replace('\\', '/')
    
    # If it's an absolute path, convert to WSL format
    if os.path.isabs(path):
        drive = path[0].lower()
        path = f"/mnt/{drive}/{path[3:]}"
    
    return path

def convert_to_windows_path(wsl_path):
    """Convert WSL path to Windows path"""
    if wsl_path.startswith('/mnt/'):
        # Extract drive letter and rest of path
        parts = wsl_path[5:].split('/', 1)
        if len(parts) > 1:
            return f"{parts[0].upper()}:\\{parts[1]}"
    return wsl_path

# Create directories on startup
ensure_directories()

@app.route('/api/health', methods=['GET'])
def health_check():
    """Simple endpoint to check if the server is running"""
    try:
        # Check if WSL is available
        wsl_status = "available" if is_wsl() else "unavailable"
        # Check if vina is installed
        vina_check = subprocess.run(
            ['which', 'vina'],
            capture_output=True,
            text=True
        )
        vina_status = "installed" if vina_check.returncode == 0 else "not installed"
        
        return jsonify({
            'status': 'ok',
            'wsl': wsl_status,
            'vina': vina_status
        })
    except Exception as e:
        return jsonify({
            'status': 'error',
            'error': str(e)
        }), 500

@app.route('/api/dock', methods=['POST'])
def dock():
    """Handle docking requests"""
    try:
        # Get JSON data from request
        data = request.get_json()
        
        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400
        
        if 'receptorFile' not in data or 'ligandFile' not in data:
            return jsonify({'error': 'Both receptor and ligand files are required'}), 400
        
        # Extract file data and names
        receptor_data = data.get('receptorFile', {}).get('data')
        receptor_name = data.get('receptorFile', {}).get('name')
        ligand_data = data.get('ligandFile', {}).get('data')
        ligand_name = data.get('ligandFile', {}).get('name')
        
        if not receptor_data or not receptor_name or not ligand_data or not ligand_name:
            return jsonify({'error': 'Invalid file data format'}), 400
        
        try:
            # Handle data URLs (with comma) or raw base64
            receptor_content = base64.b64decode(receptor_data.split(',')[1] if ',' in receptor_data else receptor_data)
            ligand_content = base64.b64decode(ligand_data.split(',')[1] if ',' in ligand_data else ligand_data)
        except Exception as e:
            return jsonify({'error': f'Base64 decoding error: {str(e)}'}), 400
        
        # Create unique filenames to prevent conflicts
        receptor_id = str(uuid.uuid4())[:8]
        ligand_id = str(uuid.uuid4())[:8]
        receptor_filename = f"{receptor_id}_{receptor_name}"
        ligand_filename = f"{ligand_id}_{ligand_name}"
        
        # Save files with WSL-compatible paths
        receptor_path = os.path.join('receptor', receptor_filename)
        ligand_path = os.path.join('ligand', ligand_filename)
        
        # Convert paths for WSL if needed
        wsl_receptor_path = convert_to_wsl_path(receptor_path)
        wsl_ligand_path = convert_to_wsl_path(ligand_path)
        
        with open(wsl_receptor_path, 'wb') as f:
            f.write(receptor_content)
        
        with open(wsl_ligand_path, 'wb') as f:
            f.write(ligand_content)
        
        # Set proper permissions for WSL
        if is_wsl():
            subprocess.run(['chmod', '666', wsl_receptor_path], check=True)
            subprocess.run(['chmod', '666', wsl_ligand_path], check=True)
        
        # Run the docking process using WSL Python
        cmd = [
            'python3' if is_wsl() else 'python',
            'dock.py',
            wsl_receptor_path,
            wsl_ligand_path
        ]
        
        # Set up environment for WSL
        env = os.environ.copy()
        if is_wsl():
            env['DISPLAY'] = ':0'
        
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        stdout, stderr = process.communicate()
        
        # Check if the process was successful
        if process.returncode == 0:
            # Extract output file path from stdout
            output_file_match = None
            for line in stdout.split('\n'):
                if "Results saved to:" in line:
                    output_file_match = line.split("Results saved to:")[1].strip()
                    break
            
            result_file = None
            result_data = None
            
            if output_file_match:
                # Convert WSL path to Windows path if needed
                wsl_output_path = convert_to_wsl_path(output_file_match)
                
                if os.path.exists(wsl_output_path):
                    # Read the result file and encode it
                    with open(wsl_output_path, 'rb') as f:
                        result_data = base64.b64encode(f.read()).decode('utf-8')
                    result_file = os.path.basename(wsl_output_path)
            
            return jsonify({
                'success': True,
                'output': stdout,
                'resultFile': result_file,
                'resultData': result_data
            })
        else:
            return jsonify({
                'success': False,
                'output': stdout,
                'error': stderr
            }), 500
            
    except Exception as e:
        app.logger.error(f"Unhandled exception: {str(e)}")
        app.logger.error(traceback.format_exc())
        return jsonify({
            'success': False,
            'error': f"Server error: {str(e)}",
            'trace': traceback.format_exc() if app.debug else None
        }), 500

@app.route('/results/<filename>', methods=['GET'])
def get_result(filename):
    """Serve result files directly"""
    try:
        # Convert path for WSL if needed
        results_dir = convert_to_wsl_path('results') if is_wsl() else 'results'
        return send_from_directory(results_dir, filename)
    except Exception as e:
        return jsonify({'error': str(e)}), 404

@app.errorhandler(404)
def not_found(e):
    return jsonify({'error': 'Not found'}), 404

@app.errorhandler(500)
def server_error(e):
    return jsonify({'error': 'Server error', 'details': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5002)