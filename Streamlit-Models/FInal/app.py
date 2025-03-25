













import streamlit as st
import pandas as pd
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from PIL import Image
import io
import re

# Import the DrugDesignRAG class - assuming it's in a file named drug_design_rag.py
from drug_design_rag import DrugDesignRAG

# Page configuration
st.set_page_config(page_title="AI Drug Designer", page_icon="💊", layout="wide")

# App title and description
st.title("AI Drug Designer")
st.markdown("""
This app uses AI to design novel drug compounds based on important substructures.
Simply provide your Google API key and enter the substructure data to get started.
""")

# API key input (use st.secrets in production)
api_key = ""

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

# Input area for substructure data
st.subheader("Enter Substructure Data")
st.markdown("Enter your substructure data or use the sample data below.")

# Option to use sample data
use_sample = st.checkbox("Use sample data", value=True)

if use_sample:
    substructure_input = st.text_area("Substructure data:", value=sample_data, height=200)
else:
    substructure_input = st.text_area("Substructure data:", height=200,
                                     placeholder="Enter data in format: Importance    Bit substructure")

# Function to display molecule from SMILES
def display_molecule(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(600, 400))
            buf = io.BytesIO()
            img.save(buf, format='PNG')
            return Image.open(buf)
        else:
            return None
    except Exception as e:
        st.error(f"Error displaying molecule: {e}")
        return None

# Function to calculate molecular properties from SMILES
def calculate_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Calculate properties directly using RDKit
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            
            # Check Lipinski's Rule of Five
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
        st.error(f"Error calculating properties: {e}")
        return None

# Extract compound name from response text
def extract_compound_name(response_text):
    # Try several patterns to extract the compound name
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
            # Clean up the name (remove periods, etc.)
            name = re.sub(r'[\.\s]+$', '', name)
            return name
    
    # Default name if none found
    return "Novel Compound"

# Generate button
if st.button("Generate Drug Compound"):
    if not api_key:
        st.error("Please enter a Google API key to proceed.")
    elif not substructure_input.strip():
        st.error("Please enter substructure data or use the sample data.")
    else:
        try:
            with st.spinner("Initializing AI model..."):
                # Set API key and initialize model
                drug_designer = DrugDesignRAG(api_key=api_key)
            
            with st.spinner("Parsing substructure data..."):
                # Parse substructures
                substructures = drug_designer.parse_substructures(substructure_input)
                
                # Display parsed substructures
                st.subheader("Parsed Substructures")
                
                # Create DataFrame for better display
                substructure_df = pd.DataFrame(
                    [(i+1, s['importance'], s['description']) for i, s in enumerate(substructures)],
                    columns=["#", "Importance", "Description"]
                )
                st.dataframe(substructure_df, hide_index=True)
            
            with st.spinner("Generating drug compound... This may take a minute or two..."):
                # Generate drug design
                result = drug_designer.generate_drug(substructures)
                
                if result["success"]:
                    st.success("Drug compound generated successfully!")
                    
                    # Create columns for results
                    col1, col2 = st.columns([3, 2])
                    
                    # Display molecular structure
                    with col1:
                        st.subheader("Molecular Structure")
                        molecule_img = display_molecule(result["smiles"])
                        if molecule_img:
                            st.image(molecule_img, use_column_width=True)
                        else:
                            st.warning("Could not display molecule structure.")
                    
                    # Display results
                    with col2:
                        st.subheader("Drug Information")
                        
                        # Extract name from response
                        drug_name = extract_compound_name(result["response"])
                        
                        st.markdown(f"**Name:** {drug_name}")
                        st.markdown(f"**SMILES Notation:**")
                        st.code(result["smiles"])
                        
                        # Display properties calculated directly from SMILES
                        st.subheader("Molecular Properties")
                        properties = calculate_properties(result["smiles"])
                        if properties:
                            prop_df = pd.DataFrame(list(properties.items()), columns=["Property", "Value"])
                            st.dataframe(prop_df, hide_index=True)
                        else:
                            st.warning("Could not calculate molecular properties.")
                    
                    # Option to download structure as image
                    if molecule_img:
                        buf = io.BytesIO()
                        molecule_img.save(buf, format='PNG')
                        buf.seek(0)
                        
                        st.download_button(
                            label="Download Molecular Structure",
                            data=buf,
                            file_name=f"{drug_name.replace(' ', '_')}_structure.png",
                            mime="image/png"
                        )
                        
                    # Show full LLM response in an expandable section
                    with st.expander("Show AI Model Response"):
                        st.markdown(result["response"])
                        
                else:
                    st.error("Failed to generate a valid drug compound.")
                    st.write(f"Error: {result.get('error', 'Unknown error')}")
                    
                    # Display partial results if available
                    if "smiles" in result and result["smiles"]:
                        st.subheader("Partial Results - SMILES")
                        st.code(result["smiles"])
                        st.warning("This SMILES string could not be validated as a proper structure.")
                    
        except Exception as e:
            st.error(f"An error occurred: {str(e)}")

# Add an information section at the bottom
st.markdown("---")
st.markdown("""
### About this App
This application uses Generative AI to design novel drug compounds based on input substructures. 
The AI considers important medicinal chemistry principles including Lipinski's Rule of Five 
and incorporates the specified substructures to create a chemically valid molecule.

### How to Use
1. Enter your Google API key
2. Provide substructure data in the format shown in the sample
3. Click "Generate Drug Compound"
4. View the resulting molecular structure and properties

### Note
The generated compounds are conceptual and would require further optimization and testing 
before they could be considered as actual drug candidates.
""")
