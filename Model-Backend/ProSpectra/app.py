from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, Crippen
from io import BytesIO
import base64
import json
import traceback
import os

# Import the custom modules
from analyzer import AdvancedMoleculeAnalyzer, count_heterocycles, generate_3d_structure

app = Flask(__name__, static_folder='../frontend/build')
CORS(app)

# Define API endpoints
@app.route('/api/analyze', methods=['POST'])
def analyze_molecule():
    try:
        data = request.json
        smiles = data.get('smiles', '')
        
        if not smiles:
            return jsonify({'error': 'No SMILES string provided'}), 400
            
        # Initialize the analyzer
        analyzer = AdvancedMoleculeAnalyzer(smiles)
        
        if not analyzer.mol:
            return jsonify({'error': 'Invalid SMILES string. Unable to process the molecule.'}), 400
            
        # Get all the data
        result = {
            'smiles': smiles,
            'descriptors': analyzer.descriptors,
            'pharmacological_properties': analyzer.analyze_pharmacological_properties(),
        }
        
        # Generate visualizations
        visualizations = analyzer.generate_molecule_visualization()
        result['visualizations'] = visualizations
        
        # Generate 3D structure
        try:
            pdb_block, view_html = generate_3d_structure(smiles)
            result['3d_structure'] = {
                'pdb_block': pdb_block,
                'view_html': view_html
            }
        except Exception as e:
            result['3d_structure'] = {
                'error': str(e)
            }
            
        return jsonify(result)
        
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

# Endpoint for example molecules
@app.route('/api/examples', methods=['GET'])
def get_example_molecules():
    example_molecules = {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",
        "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    }
    return jsonify(example_molecules)

# Serve React App
@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def serve(path):
    if path != "" and os.path.exists(app.static_folder + '/' + path):
        return send_from_directory(app.static_folder, path)
    else:
        return send_from_directory(app.static_folder, 'index.html')

if __name__ == '__main__':
    app.run(debug=True, port=5003)