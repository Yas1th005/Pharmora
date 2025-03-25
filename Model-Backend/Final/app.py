from flask import Flask, request, jsonify
from flask_cors import CORS
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from PIL import Image
import io
import re
import base64  # Import the base64 module
from drug_design_rag import DrugDesignRAG

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Initialize the DrugDesignRAG instance
api_key = ""
drug_designer = DrugDesignRAG(api_key=api_key)

# Sample data
sample_data = """
Importance    Bit substructure
0.0416       >= 5 unsaturated non-aromatic heteroatom-containing ring size 6
0.0303       >= 5 unsaturated non-aromatic heteroatom-containing ring size 5
0.0301       >= 2 saturated or aromatic heteroatom-containing ring size 4
0.0190       >= 1 Bi
0.0184       >= 3 saturated or aromatic nitrogen-containing ring size 6
0.0152       >= 1 Ni
"""

@app.route('/generate-drug', methods=['POST'])
def generate_drug():
    data = request.json
    substructure_input = data.get('substructure_input', sample_data)
    
    try:
        # Parse substructures
        substructures = drug_designer.parse_substructures(substructure_input)
        
        # Generate drug design
        result = drug_designer.generate_drug(substructures)
        
        if result["success"]:
            # Convert image to base64 for frontend display
            img = Draw.MolToImage(Chem.MolFromSmiles(result["smiles"]), size=(600, 400))
            buffered = io.BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()  # Encode to base64
            
            response = {
                "success": True,
                "drug_name": extract_compound_name(result["response"]),
                "smiles": result["smiles"],
                "image": img_str,
                "properties": calculate_properties(result["smiles"]),
                "response": result["response"]
            }
        else:
            response = {
                "success": False,
                "error": result.get("error", "Unknown error"),
                "response": result.get("response", "")
            }
        
        return jsonify(response)
    
    except Exception as e:
        return jsonify({"success": False, "error": str(e)})

def extract_compound_name(response_text):
    patterns = [
        r"Name:[\s]*([\w\-\s]+)",
        r"compound[\s]*name:[\s]*([\w\-\s]+)",
        r"named[\s]*([\w\-\s]+)",
        r"called[\s]*([\w\-\s]+)"
    ]
    
    for pattern in patterns:
        match = re.search(pattern, response_text, re.IGNORECASE)
        if match:
            name = match.group(1).strip()
            name = re.sub(r'[\.\s]+$', '', name)
            return name
    
    return "Novel Compound"

def calculate_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            
            lipinski_violations = 0
            if mw > 500: lipinski_violations += 1
            if logp > 5: lipinski_violations += 1
            if hbd > 5: lipinski_violations += 1
            if hba > 10: lipinski_violations += 1
            
            properties = {
                "Molecular Weight": f"{mw:.2f} g/mol",
                "LogP": f"{logp:.2f}",
                "Hydrogen Bond Donors": f"{hbd}",
                "Hydrogen Bond Acceptors": f"{hba}",
                "Rotatable Bonds": f"{rotatable_bonds}",
                "Topological Polar Surface Area": f"{tpsa:.2f} Å²",
                "Lipinski Violations": f"{lipinski_violations}/4"
            }
            
            return properties
        else:
            return None
    except Exception as e:
        return None

if __name__ == '__main__':
    app.run(debug=True,port=5005)
