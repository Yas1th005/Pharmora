import streamlit as st
import pandas as pd
import numpy as np
import mols2grid
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit.Chem.Descriptors import (
    ExactMolWt, MolLogP, NumHDonors, NumHAcceptors, 
    TPSA, NumRotatableBonds, MolWt
)
from rdkit.Chem import Crippen
import plotly.express as px
import plotly.graph_objs as go

@st.cache_data
def calculate_advanced_descriptors(smiles_string):
    """Calculate comprehensive molecular descriptors"""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return None
    
    try:
        return {
            # Basic Lipinski Descriptors
            'Molecular Weight': Descriptors.ExactMolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'H-Bond Donors': Descriptors.NumHDonors(mol),
            'H-Bond Acceptors': Descriptors.NumHAcceptors(mol),
            
            # Advanced Physicochemical Properties
            'Topological Polar Surface Area': Descriptors.TPSA(mol),
            'Rotatable Bonds': Descriptors.NumRotatableBonds(mol),
            'Molecular Complexity': Descriptors.BertzCT(mol),
            
            # Pharmacokinetic Indicators
            'Molar Refractivity': Crippen.MolMR(mol),
            'Formal Charge': Chem.GetFormalCharge(mol),
            
            # Structural Indicators
            'Ring Count': Descriptors.RingCount(mol),
            'Aromatic Ring Count': Descriptors.NumAromaticRings(mol),
            'Fraction SP3 Carbon': Descriptors.FractionCSP3(mol),
        }
    except Exception as e:
        st.warning(f"Error calculating descriptors for {smiles_string}: {e}")
        return None

def main():
    st.set_page_config(
        page_title="Advanced Drug Discovery Toolkit", 
        page_icon="💊", 
        layout="wide"
    )
    
    # Sidebar Navigation
    st.sidebar.title("🔬 Drug Discovery Toolkit")
    app_mode = st.sidebar.selectbox(
        "Choose Analysis Mode",
        [
            "Molecular Property Filtering", 
            "Descriptor Distribution", 
            "Chemical Space Visualization",
        ]
    )
    
    # Load dataset
    @st.cache_data
    def load_comprehensive_dataset():
        try:
            df = pd.read_csv(
                "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", 
                sep="\t"
            ).dropna()
            
            # Calculate advanced descriptors
            descriptors = df['smiles'].apply(calculate_advanced_descriptors).apply(pd.Series)
            df = pd.concat([df, descriptors], axis=1)
            
            return df
        except Exception as e:
            st.error(f"Dataset Loading Error: {e}")
            return pd.DataFrame()
    
    df = load_comprehensive_dataset()
    
    if df.empty:
        st.error("Could not load dataset. Please check your connection.")
        return
    
    # Analysis Modes
    if app_mode == "Molecular Property Filtering":
        st.header("🔍 Molecular Property Filter")
        
        # Multi-parameter Filtering
        cols = st.columns(5)
        with cols[0]:
            mw_range = st.slider(
                "Molecular Weight Range", 
                min_value=0, 
                max_value=1000, 
                value=(200, 500)
            )
        
        with cols[1]:
            logp_range = st.slider(
                "LogP Range", 
                min_value=-5.0, 
                max_value=10.0, 
                value=(-0.4, 5.5),
                step=0.1
            )
        
        with cols[2]:
            tpsa_range = st.slider(
                "Topological Polar Surface Area", 
                min_value=0, 
                max_value=200, 
                value=(0, 140)
            )
        
        with cols[3]:
            hdonors_range = st.slider(
                "H-Bond Donors", 
                min_value=0, 
                max_value=15, 
                value=(0, 5)
            )
        
        with cols[4]:
            hacceptors_range = st.slider(
                "H-Bond Acceptors", 
                min_value=0, 
                max_value=20, 
                value=(0, 10)
            )
        
        # Compound Search
        search_term = st.text_input("Search for a Compound (by name or SMILES)")
        
        # Advanced Filtering
        filtered_df = df[
            (df['Molecular Weight'].between(mw_range[0], mw_range[1])) &
            (df['LogP'].between(logp_range[0], logp_range[1])) &
            (df['Topological Polar Surface Area'].between(tpsa_range[0], tpsa_range[1])) &
            (df['H-Bond Donors'].between(hdonors_range[0], hdonors_range[1])) &
            (df['H-Bond Acceptors'].between(hacceptors_range[0], hacceptors_range[1]))
        ]
        
        # Apply search filter if a search term is provided
        if search_term:
            filtered_df = filtered_df[
                filtered_df['generic_name'].str.contains(search_term, case=False, na=False) |
                filtered_df['smiles'].str.contains(search_term, case=False, na=False)
            ]
        
        st.write(f"Filtered Compounds: {len(filtered_df)}")
        
        # Display filtered compounds in a grid using mols2grid
        if not filtered_df.empty:
            # Prepare DataFrame for mols2grid
            display_df = filtered_df.rename(columns={
                'smiles': 'SMILES', 
                'generic_name': 'Name'
            })
            
            # Generate mols2grid HTML
            raw_html = mols2grid.display(
                display_df, 
                subset=["img", "Name", "SMILES", "Molecular Weight", "LogP", "H-Bond Donors", "H-Bond Acceptors"],
                mapping={
                    "SMILES": "SMILES", 
                    "Name": "Name"
                }
            )._repr_html_()
            
            # Display molecular grid
            components.html(raw_html, width=1000, height=1200, scrolling=True)
        else:
            st.warning("No compounds match the selected criteria.")
    
    elif app_mode == "Descriptor Distribution":
        st.header("📊 Molecular Descriptor Distribution")
        
        # Select descriptor for distribution
        descriptor = st.selectbox(
            "Choose Descriptor for Distribution",
            [
                'Molecular Weight', 'LogP', 
                'H-Bond Donors', 'H-Bond Acceptors', 
                'Topological Polar Surface Area'
            ]
        )
        
        # Interactive Histogram
        fig = px.histogram(
            df, 
            x=descriptor, 
            title=f'Distribution of {descriptor}',
            marginal='box'
        )
        st.plotly_chart(fig)
    
    elif app_mode == "Chemical Space Visualization":
        st.header("🌈 Chemical Space Visualization")
        
        # 2D Scatter Plot of Key Descriptors
        x_axis = st.selectbox("X-Axis", ['LogP', 'Molecular Weight'], index=0)
        y_axis = st.selectbox("Y-Axis", ['Topological Polar Surface Area', 'H-Bond Donors'], index=0)
        
        fig = px.scatter(
            df, 
            x=x_axis, 
            y=y_axis, 
            color='Aromatic Ring Count',
            hover_data=['generic_name', 'smiles'],
            title=f'Chemical Space: {x_axis} vs {y_axis}'
        )
        st.plotly_chart(fig)
    
    # Add more modes as needed
    st.sidebar.markdown("---")
    st.sidebar.info(
        "🚀 Explore. Analyze. Discover.\n"
        "Advanced drug discovery toolkit for computational chemists."
    )

if __name__ == "__main__":
    main()