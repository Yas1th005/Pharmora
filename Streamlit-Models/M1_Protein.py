import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re
import requests
import biotite.structure.io as bsio
import py3Dmol
from stmol import showmol
import plotly.express as px
import plotly.graph_objects as go
from io import BytesIO

# Expanded Hydrophobicity and Charge Scales
HYDROPHOBICITY_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

CHARGE_SCALE = {
    'D': -1, 'E': -1, 'H': 0.1, 'K': 1, 'R': 1, 'C': -0.1, 'Y': -0.1
}

# Advanced Sequences for Research
DEFAULT_SEQUENCES = {
    "Default": "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ",
    "Insulin": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
    "P53 Tumor Suppressor": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPGSSRPQKK"
}

def calculate_hydrophobicity(sequence, window_size=9):
    """Calculate hydrophobicity for a given protein sequence."""
    hydrophobicity = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        avg_hydro = sum(HYDROPHOBICITY_SCALE[aa] for aa in window) / window_size
        hydrophobicity.append(avg_hydro)
    return hydrophobicity

def calculate_charge(sequence):
    """Calculate net charge of a protein sequence."""
    charge = 0
    for aa in sequence:
        charge += CHARGE_SCALE.get(aa, 0)
    return charge

def predict_protein_structure(sequence):
    """Predict protein structure using ESMAtlas API."""
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    pdb_string = response.content.decode('utf-8')
    
    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)
    
    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 4)
    
    return pdb_string, b_value

def advanced_sequence_analysis(sequence):
    """Comprehensive protein sequence analysis."""
    # Amino Acid Patterns
    motifs = {
        'Nuclear Localization Signal (NLS)': r'[KR][KR]X[KR]',
        'Glycosylation Sites': r'N[^P][ST]',
        'Phosphorylation Sites': r'[ST].[ST]',
        'Zinc Finger Motif': r'C.{2,4}C.{3}[LIVMC].{8}H.{3,5}H'
    }
    
    found_motifs = {}
    for name, pattern in motifs.items():
        found_motifs[name] = [m.start() for m in re.finditer(pattern, sequence)]
    
    # Secondary Structure Prediction
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

def render_advanced_mol(pdb):
    """Enhanced molecular rendering with more visualization options."""
    pdbview = py3Dmol.view(width=1000, height=600)
    pdbview.addModel(pdb, 'pdb')
    
    # Multiple rendering styles
    pdbview.setStyle(
        {'cartoon': {'color': 'spectrum'}},  # Rainbow cartoon
    )
    pdbview.addStyle(
        {'hetflag': True}, 
        {'stick': {'color': 'green'}}  # Highlight hetero atoms
    )
    
    pdbview.setBackgroundColor('black')
    pdbview.zoomTo()
    pdbview.spin(True)
    return pdbview

def main():
    # Advanced UI Configuration
    st.set_page_config(
        page_title="BioFold: Advanced Protein Analysis Platform",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    # Custom CSS for enhanced UI
    st.markdown("""
    <style>
    .reportview-container {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    }
    .sidebar .sidebar-content {
        background: rgba(255,255,255,0.8);
        backdrop-filter: blur(10px);
    }
    .stButton>button {
        background-color: #4CAF50;
        color: white;
        transition: all 0.3s ease;
    }
    .stButton>button:hover {
        background-color: #45a049;
        transform: scale(1.05);
    }
    </style>
    """, unsafe_allow_html=True)

    # Title and Description
    st.title("🧬 BioFold: Advanced Protein Analysis Platform")
    st.markdown("### Comprehensive Protein Sequence & Structure Insights")

    # Sidebar with Enhanced Input
    with st.sidebar:
        st.header("🔬 Protein Sequence Input")
        
        # Sequence Selection
        sequence_type = st.selectbox(
            "Select Protein Sequence", 
            list(DEFAULT_SEQUENCES.keys()) + ["Custom Sequence"]
        )
        
        # Sequence Input
        if sequence_type == "Custom Sequence":
            protein_sequence = st.text_area(
                'Enter Protein Sequence', 
                height=200, 
                placeholder="Paste your protein sequence here..."
            )
        else:
            protein_sequence = DEFAULT_SEQUENCES[sequence_type]
        
        # Analysis Options
        st.subheader("Analysis Settings")
        analysis_options = st.multiselect(
            "Select Analysis Types",
            [
                "Structure Prediction", 
                "Sequence Composition", 
                "Hydrophobicity", 
                "Motif Analysis", 
                "Secondary Structure"
            ],
            default=[
                "Structure Prediction", 
                "Sequence Composition", 
                "Hydrophobicity"
            ]
        )
        
        # Analyze Button
        analyze_button = st.button('🔍 Analyze Protein', use_container_width=True)

    # Main Analysis
    if analyze_button:
        # Validation
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        if not set(protein_sequence).issubset(valid_amino_acids):
            st.error("Invalid protein sequence. Please use standard amino acid letters.")
            return

        # Tabbed Interface for Results
        tab1, tab2, tab3 = st.tabs([
            "🏗️ Structure", 
            "📊 Composition", 
            "🔬 Advanced Analysis"
        ])

        with tab1:
            if "Structure Prediction" in analysis_options:
                st.header("Protein Structure Prediction")
                try:
                    pdb_string, b_value = predict_protein_structure(protein_sequence)
                    
                    # Advanced 3D Visualization
                    mol_view = render_advanced_mol(pdb_string)
                    showmol(mol_view, height=600, width=1000)
                    
                    st.subheader('Prediction Confidence')
                    st.metric(
                        label="plDDT Score", 
                        value=f"{b_value}/100", 
                        help="Higher score indicates more confident prediction"
                    )

                    # Download PDB Button
                    st.download_button(
                        label="Download PDB File",
                        data=pdb_string,
                        file_name='predicted_structure.pdb',
                        mime='text/plain',
                    )

                except Exception as e:
                    st.error(f"Structure prediction error: {e}")

        with tab2:
            if any(opt in analysis_options for opt in ["Sequence Composition", "Hydrophobicity"]):
                col1, col2 = st.columns(2)
                
                with col1:
                    # Amino Acid Composition
                    st.subheader("Amino Acid Composition")
                    composition = Counter(protein_sequence)
                    fig = px.pie(
                        values=list(composition.values()),
                        names=list(composition.keys()),
                        title="Amino Acid Distribution"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Hydrophobicity Plot
                    st.subheader("Hydrophobicity Analysis")
                    hydro_values = calculate_hydrophobicity(protein_sequence)
                    fig = px.line(
                        y=hydro_values, 
                        title="Hydrophobicity Landscape",
                        labels={'index':'Residue Position', 'y':'Hydrophobicity'}
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                # Additional Molecular Properties
            st.subheader("🧪 Molecular Characteristics")
            
            # Create a more detailed analysis section
            analyzed_seq = ProteinAnalysis(protein_sequence)
            
            # Grid of molecular properties
            prop_cols = st.columns(4)
            properties = [
                ("Molecular Weight", f"{analyzed_seq.molecular_weight():.2f} Da"),
                ("Isoelectric Point", f"{analyzed_seq.isoelectric_point():.2f}"),
                ("Net Charge", f"{calculate_charge(protein_sequence)}"),
                ("Instability Index", f"{analyzed_seq.instability_index():.2f}")
            ]
            
            for col, (name, value) in zip(prop_cols, properties):
                with col:
                    st.metric(label=name, value=value)

            # Amino Acid Composition Heatmap
            st.write(" ")
            st.subheader("Amino Acid Composition Heatmap")
            composition = Counter(protein_sequence)
            
            # Create a DataFrame for heatmap
            aa_order = 'ACDEFGHIKLMNPQRSTVWY'
            comp_matrix = [[composition.get(aa, 0) for aa in aa_order]]
            
            fig = px.imshow(
                comp_matrix, 
                labels=dict(x="Amino Acids", y="Frequency"),
                x=list(aa_order),
                color_continuous_scale='Viridis'
            )
            fig.update_layout(height=200)
            st.plotly_chart(fig, use_container_width=True)

        # Replace the existing tab3 content in the main() function with this improved version:

        with tab3:
            st.header("🔬 Advanced Molecular Insights")
            
            # Create columns for better layout
            col1, col2 = st.columns([2, 1])
            
            # Advanced Analysis
            advanced_analysis = advanced_sequence_analysis(protein_sequence)
            
            with col1:
                # Motif Detection with Enhanced Visualization
                st.subheader("🔍 Motif Detection")
                
                # Create an expandable section for each motif type
                motifs = advanced_analysis['motifs']
                
                # Use a card-like layout for motifs
                for motif, positions in motifs.items():
                    with st.expander(f"**{motif}**"):
                        if positions:
                            # Create a more detailed view of motif positions
                            df_motifs = pd.DataFrame({
                                'Position': positions,
                                'Sequence Snippet': [
                                    protein_sequence[max(0, pos-5):pos+5] for pos in positions
                                ]
                            })
                            
                            st.dataframe(
                                df_motifs, 
                                use_container_width=True,
                                hide_index=True,
                                column_config={
                                    "Position": st.column_config.NumberColumn(
                                        "Residue Position",
                                        help="Location of the motif in the protein sequence"
                                    ),
                                    "Sequence Snippet": st.column_config.TextColumn(
                                        "Context",
                                        help="5 residues before and after the motif"
                                    )
                                }
                            )
                        else:
                            st.info(f"No {motif} found in the sequence.")

            with col2:
                # Secondary Structure Prediction
                st.subheader("🧩 Secondary Structure")
                
                # Enhanced Pie Chart with more details
                structure_data = advanced_analysis['secondary_structure']
                
                # Create a more informative visualization
                fig = go.Figure(data=[go.Pie(
                    labels=list(structure_data.keys()),
                    values=list(structure_data.values()),
                    textinfo='label+percent',
                    insidetextorientation='radial',
                    marker_colors=['#FF6384', '#36A2EB', '#FFCE56'],  # Custom color palette
                )])
                
                fig.update_layout(
                    title_text="Secondary Structure Composition",
                    title_font_size=16,
                    margin=dict(t=50, b=0, l=0, r=0)
                )
                
                st.plotly_chart(fig, use_container_width=True)

if __name__ == "__main__":
    main()