from flask import Flask, request, jsonify
from flask_cors import CORS
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re
import requests
import biotite.structure.io as bsio
from collections import Counter
import numpy as np
from io import BytesIO

app = Flask(__name__)
CORS(app)

# Keep the same constants from your Streamlit app
HYDROPHOBICITY_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

CHARGE_SCALE = {
    'D': -1, 'E': -1, 'H': 0.1, 'K': 1, 'R': 1, 'C': -0.1, 'Y': -0.1
}

DEFAULT_SEQUENCES = {
    "Default": "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ",
    "Insulin": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
    "P53 Tumor Suppressor": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPGSSRPQKK"
}

def calculate_hydrophobicity(sequence, window_size=9):
    hydrophobicity = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        avg_hydro = sum(HYDROPHOBICITY_SCALE[aa] for aa in window) / window_size
        hydrophobicity.append(avg_hydro)
    return hydrophobicity

def calculate_charge(sequence):
    return sum(CHARGE_SCALE.get(aa, 0) for aa in sequence)

def predict_protein_structure(sequence):
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    pdb_string = response.content.decode('utf-8')
    
    with open('./results/predicted.pdb', 'w') as f:
        f.write(pdb_string)
    
    struct = bsio.load_structure('./results/predicted.pdb', extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 4)
    
    return pdb_string, b_value

def advanced_sequence_analysis(sequence):
    motifs = {
        'Nuclear Localization Signal (NLS)': r'[KR][KR]X[KR]',
        'Glycosylation Sites': r'N[^P][ST]',
        'Phosphorylation Sites': r'[ST].[ST]',
        'Zinc Finger Motif': r'C.{2,4}C.{3}[LIVMC].{8}H.{3,5}H'
    }
    
    found_motifs = {}
    for name, pattern in motifs.items():
        found_motifs[name] = [m.start() for m in re.finditer(pattern, sequence)]
    
    try:
        prot_analysis = ProteinAnalysis(sequence)
        helix = prot_analysis.secondary_structure_fraction()[0]
        sheet = prot_analysis.secondary_structure_fraction()[1]
        coil = prot_analysis.secondary_structure_fraction()[2]
    except:
        helix, sheet, coil = 0, 0, 0
    
    return {
        'motifs': found_motifs,
        'secondary_structure': {
            'Helix': helix,
            'Sheet': sheet, 
            'Coil': coil
        }
    }

@app.route('/api/analyze', methods=['POST'])
def analyze_protein():
    data = request.get_json()
    sequence = data.get('sequence')
    analysis_options = data.get('analysis_options', [])
    
    if not sequence:
        return jsonify({'error': 'No sequence provided'}), 400
    
    # Validate sequence
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    if not set(sequence).issubset(valid_amino_acids):
        return jsonify({'error': 'Invalid protein sequence'}), 400
    
    result = {}
    
    # Structure prediction
    if "Structure Prediction" in analysis_options:
        try:
            pdb_string, b_value = predict_protein_structure(sequence)
            result['structure'] = {
                'pdb_string': pdb_string,
                'b_value': b_value
            }
        except Exception as e:
            result['structure'] = {'error': str(e)}
    
    # Sequence composition
    composition = Counter(sequence)
    result['composition'] = {aa: count for aa, count in composition.items()}
    
    # Hydrophobicity
    hydro_values = calculate_hydrophobicity(sequence)
    result['hydrophobicity'] = hydro_values
    
    # Advanced analysis
    advanced_results = advanced_sequence_analysis(sequence)
    result['advanced_analysis'] = advanced_results
    
    # Molecular properties
    try:
        analyzed_seq = ProteinAnalysis(sequence)
        result['molecular_properties'] = {
            'molecular_weight': analyzed_seq.molecular_weight(),
            'isoelectric_point': analyzed_seq.isoelectric_point(),
            'net_charge': calculate_charge(sequence),
            'instability_index': analyzed_seq.instability_index()
        }
    except Exception as e:
        result['molecular_properties'] = {'error': str(e)}
    
    return jsonify(result)

@app.route('/api/sequences', methods=['GET'])
def get_sequences():
    return jsonify(DEFAULT_SEQUENCES)

if __name__ == '__main__':
    app.run(debug=True,port=5000)