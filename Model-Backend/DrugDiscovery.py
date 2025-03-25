from flask import Flask, jsonify, request
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit.Chem.Descriptors import (
    ExactMolWt, MolLogP, NumHDonors, NumHAcceptors, 
    TPSA, NumRotatableBonds, MolWt
)
from rdkit.Chem import Crippen
import pandas as pd
from flask_cors import CORS
import base64
from io import BytesIO

app = Flask(__name__)
CORS(app)

def calculate_advanced_descriptors(smiles_string):
    """Calculate comprehensive molecular descriptors"""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return None
    
    try:
        # Generate molecule image
        img = Draw.MolToImage(mol)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        return {
            'image': img_str,
            'Molecular Weight': round(Descriptors.ExactMolWt(mol), 2),
            'LogP': round(Descriptors.MolLogP(mol), 2),
            'H-Bond Donors': Descriptors.NumHDonors(mol),
            'H-Bond Acceptors': Descriptors.NumHAcceptors(mol),
            'Topological Polar Surface Area': round(Descriptors.TPSA(mol), 2),
            'Rotatable Bonds': Descriptors.NumRotatableBonds(mol),
            'Molecular Complexity': round(Descriptors.BertzCT(mol), 2),
            'Molar Refractivity': round(Crippen.MolMR(mol), 2),
            'Formal Charge': Chem.GetFormalCharge(mol),
            'Ring Count': Descriptors.RingCount(mol),
            'Aromatic Ring Count': Descriptors.NumAromaticRings(mol),
            'Fraction SP3 Carbon': round(Descriptors.FractionCSP3(mol), 2),
        }
    except Exception as e:
        return None

@app.route('/api/compounds', methods=['GET'])
def get_compounds():
    try:
        # Load and process dataset
        df = pd.read_csv(
            "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", 
            sep="\t"
        ).dropna()
        
        # Get filter parameters
        mw_min = float(request.args.get('mw_min', 0))
        mw_max = float(request.args.get('mw_max', 1000))
        logp_min = float(request.args.get('logp_min', -5))
        logp_max = float(request.args.get('logp_max', 10))
        tpsa_min = float(request.args.get('tpsa_min', 0))
        tpsa_max = float(request.args.get('tpsa_max', 200))
        hbd_min = float(request.args.get('hbd_min', 0))
        hbd_max = float(request.args.get('hbd_max', 15))
        hba_min = float(request.args.get('hba_min', 0))
        hba_max = float(request.args.get('hba_max', 20))
        search_term = request.args.get('search', '')
        
        # Calculate descriptors for all compounds
        compounds = []
        for _, row in df.iterrows():
            descriptors = calculate_advanced_descriptors(row['smiles'])
            if descriptors:
                compound = {
                    'name': row['generic_name'],
                    'smiles': row['smiles'],
                    **descriptors
                }
                
                # Apply all filters
                if (mw_min <= compound['Molecular Weight'] <= mw_max and
                    logp_min <= compound['LogP'] <= logp_max and
                    tpsa_min <= compound['Topological Polar Surface Area'] <= tpsa_max and
                    hbd_min <= compound['H-Bond Donors'] <= hbd_max and
                    hba_min <= compound['H-Bond Acceptors'] <= hba_max and
                    (not search_term or 
                     search_term.lower() in compound['name'].lower() or
                     search_term.lower() in compound['smiles'].lower())):
                    compounds.append(compound)
        
        return jsonify(compounds)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/descriptor-stats', methods=['GET'])
def get_descriptor_stats():
    try:
        df = pd.read_csv(
            "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", 
            sep="\t"
        ).dropna()
        
        descriptor = request.args.get('descriptor', 'Molecular Weight')
        
        # Calculate descriptors for all compounds
        values = []
        for _, row in df.iterrows():
            descriptors = calculate_advanced_descriptors(row['smiles'])
            if descriptors and descriptor in descriptors:
                values.append(descriptors[descriptor])
        
        return jsonify({
            'values': values,
            'mean': sum(values) / len(values),
            'min': min(values),
            'max': max(values)
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True,port=5001)