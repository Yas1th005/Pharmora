import streamlit as st
import pandas as pd
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, Crippen
from io import BytesIO
import base64
import py3Dmol

def count_heterocycles(mol):
    """
    Count the number of heterocyclic rings in the molecule.
    
    Parameters:
    - mol: RDKit Mol object
    
    Returns:
    - Number of heterocyclic rings
    """
    if mol is None:
        return 0

    # Get the molecule's rings
    rings = mol.GetRingInfo()
    atom_rings = rings.AtomRings()
    heterocyclic_count = 0

    for ring in atom_rings:
        # Check if the ring contains any heteroatoms
        if any(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() not in [6, 1] for atom_idx in ring):
            heterocyclic_count += 1
    
    return heterocyclic_count

class AdvancedMoleculeAnalyzer:
    def __init__(self, smiles):
        """
        Comprehensive molecule analysis toolkit
        
        Parameters:
        - smiles: Input SMILES string
        """
        self.original_smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        self.descriptors = self._generate_comprehensive_descriptors()

    def _generate_comprehensive_descriptors(self):
        """
        Generate an extensive set of molecular descriptors
        
        Returns:
        - Dictionary of molecular descriptors
        """
        if self.mol is None:
            return {}
        
        # Comprehensive descriptor calculation
        descriptors = {
            # Structural Properties
            "Molecular Formula": Chem.rdMolDescriptors.CalcMolFormula(self.mol),
            "Molecular Weight": round(Descriptors.ExactMolWt(self.mol), 2),
            "Heavy Atom Count": Descriptors.HeavyAtomCount(self.mol),
            "Total Atom Count": self.mol.GetNumAtoms(),
            "Total Bond Count": self.mol.GetNumBonds(),
            
            # Chemical Properties
            "LogP": round(Descriptors.MolLogP(self.mol), 2),
            "H-Bond Donors": Descriptors.NumHDonors(self.mol),
            "H-Bond Acceptors": Descriptors.NumHAcceptors(self.mol),
            "Rotatable Bonds": Descriptors.NumRotatableBonds(self.mol),
            "Topological Polar Surface Area": round(Descriptors.TPSA(self.mol), 2),
            
            # Advanced Structural Descriptors
            "Ring Count": Descriptors.RingCount(self.mol),
            "Aromatic Ring Count": Descriptors.NumAromaticRings(self.mol),
            "Saturated Ring Count": Descriptors.NumSaturatedRings(self.mol),
            "Heterocyclic Ring Count": count_heterocycles(self.mol),
            
            # Pharmacokinetic Indicators
            "Molar Refractivity": round(Crippen.MolMR(self.mol), 2),
            "Fraction SP3 Carbons": round(Descriptors.FractionCSP3(self.mol), 2),
        }
        
        # Lipinski's Rule of Five Evaluation
        lipinski_violations = 0
        if descriptors["Molecular Weight"] > 500: lipinski_violations += 1
        if descriptors["LogP"] > 5: lipinski_violations += 1
        if descriptors["H-Bond Donors"] > 5: lipinski_violations += 1
        if descriptors["H-Bond Acceptors"] > 10: lipinski_violations += 1
        
        descriptors["Lipinski Violations"] = lipinski_violations
        
        return descriptors
    
    def generate_molecule_visualization(self, img_size=(600, 600)):
        """
        Generate multiple molecule visualizations
        
        Returns:
        - Dictionary of molecule image representations
        """
        if self.mol is None:
            return {}
        
        visualizations = {}
        
        # 2D Structure
        img_2d = Draw.MolToImage(self.mol, size=img_size)
        buffered = BytesIO()
        img_2d.save(buffered, format="PNG")
        visualizations["2D Structure"] = base64.b64encode(buffered.getvalue()).decode()
        
        # 2D Detailed Structure
        img_2d_detailed = Draw.MolToImage(
            self.mol, 
            size=img_size, 
            kekulize=True, 
            wedgeBonds=True, 
            fitImage=True
        )
        buffered = BytesIO()
        img_2d_detailed.save(buffered, format="PNG")
        visualizations["2D Detailed Structure"] = base64.b64encode(buffered.getvalue()).decode()
        
        return visualizations
    
    def _compute_synthetic_accessibility(self):
        """
        Compute the synthetic accessibility score using SA_Score.
        
        Returns:
        - Synthetic accessibility score
        """
        try:
            from sascorer import calculateScore
        except ImportError:
            return "SA_Score module not installed. Please install sascorer."

        if self.mol is None:
            return None

        return round(calculateScore(self.mol), 2)
    
    def analyze_pharmacological_properties(self):
        """
        Assess potential pharmacological characteristics
        
        Returns:
        - Dictionary of pharmacological predictions
        """
        if self.mol is None:
            return {}
        
        # Preliminary pharmacological assessment
        pharmacological_props = {
            "Drug-likeness": self._assess_drug_likeness(),
            "Bioavailability": self._predict_bioavailability(),
            "Synthetic Accessibility": self._compute_synthetic_accessibility(),
            "Medicinal Chemistry Friendliness": self._medicinal_chemistry_assessment()
        }
        
        return pharmacological_props
    
    def _assess_drug_likeness(self):
        """
        Evaluate drug-likeness based on multiple criteria
        
        Returns:
        - Drug-likeness score and category
        """
        violations = self.descriptors.get("Lipinski Violations", 0)
        
        if violations == 0:
            return "Excellent (Follows Lipinski's Rule of Five)"
        elif violations <= 2:
            return "Good (Minor Violations)"
        else:
            return "Poor (Significant Rule Violations)"
    
    def _predict_bioavailability(self):
        """
        Predict oral bioavailability
        
        Returns:
        - Bioavailability prediction
        """
        # Simplified bioavailability prediction
        logP = self.descriptors.get("LogP", 0)
        mol_weight = self.descriptors.get("Molecular Weight", 0)
        h_donors = self.descriptors.get("H-Bond Donors", 0)
        h_acceptors = self.descriptors.get("H-Bond Acceptors", 0)
        
        # Veber's Rule and additional criteria
        if (mol_weight <= 500 and 
            logP <= 5 and 
            h_donors <= 5 and 
            h_acceptors <= 10 and 
            Descriptors.NumRotatableBonds(self.mol) <= 10):
            return "High Probability of Oral Bioavailability"
        else:
            return "Moderate to Low Oral Bioavailability"
    
    def _medicinal_chemistry_assessment(self):
        """
        Assess molecule's suitability for medicinal chemistry
        
        Returns:
        - Medicinal chemistry friendliness assessment
        """
        # Check for problematic structural features
        problematic_features = 0
        
        # Check for excessive aromatic rings
        if self.descriptors.get("Aromatic Ring Count", 0) > 4:
            problematic_features += 1
        
        # Check for extremely hydrophobic molecules
        if self.descriptors.get("LogP", 0) > 6:
            problematic_features += 1
        
        # Evaluate based on number of problematic features
        if problematic_features == 0:
            return "Highly Suitable for Medicinal Chemistry"
        elif problematic_features == 1:
            return "Moderately Suitable"
        else:
            return "Challenging for Medicinal Chemistry"
        
def generate_3d_structure(smiles):
    """
    Generate a 3D molecular structure using RDKit and py3Dmol
    
    Parameters:
    - smiles: Input SMILES string
    
    Returns:
    - PDB block of the molecule
    - py3Dmol view HTML
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid molecule structure"
    
    # Add hydrogens and generate 3D coordinates
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    # Convert to PDB
    pdb_block = Chem.MolToPDBBlock(mol)
    
    # Create py3Dmol view
    view = py3Dmol.view(width=470, height=470)
    view.addModel(pdb_block, 'pdb')
    view.setStyle({'stick':{}})
    view.zoomTo()
    
    # Generate HTML representation
    view_html = view._make_html()
    
    return pdb_block, view_html

def generate_fallback_3d_image(smiles):
    """
    Generate a fallback 2D image of the molecule as a backup
    
    Parameters:
    - smiles: Input SMILES string
    
    Returns:
    - Base64 encoded image
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Generate 2D depiction
    img = Draw.MolToImage(mol, size=(470, 470))
    
    # Convert to base64
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()

def add_3d_visualization(st, analyzer):
    """
    Add 3D molecular structure visualization
    
    Parameters:
    - st: Streamlit module
    - analyzer: AdvancedMoleculeAnalyzer instance
    """
    
    # Attempt to generate 3D structure
    try:
        pdb_block, view_html = generate_3d_structure(analyzer.original_smiles)
        
        # Columns for visualization
        col1, col2 = st.columns(2)
        
        with col1:
            # 3D Structure Visualization
            st.subheader("3D Interactive Structure")
            
            # Attempt to render 3D view
            try:
                # Use streamlit components to render the 3D view
                st.components.v1.html(view_html, height=470, width=470)
            except Exception as render_error:
                # Fallback to 2D image if 3D rendering fails
                st.warning("3D visualization failed. Showing 2D representation.")
                fallback_img = generate_fallback_3d_image(analyzer.original_smiles)
                if fallback_img:
                    st.image(f"data:image/png;base64,{fallback_img}", use_column_width=True)
        
        
        
        # Molecular Structure Details
        st.subheader("Molecular Structure Insights")
        st.markdown(f"""
        - **SMILES**: `{analyzer.original_smiles}`
        - **Molecular Weight**: {analyzer.descriptors.get('Molecular Weight', 'N/A')} g/mol
        - **Total Atoms**: {analyzer.descriptors.get('Total Atom Count', 'N/A')}
        """)
        
        # Additional Insights
        st.info("""
        🔬 **Structural Insights**:
        - 3D structure represents the molecule's optimized geometric configuration
        - Hydrogen atoms are added to complete the structural representation
        - Energy minimization performed using MMFF force field
        """)
    
    except Exception as e:
        st.error(f"Error generating molecular visualization: {str(e)}")

def main():
    # Page Configuration
    st.set_page_config(
        page_title="🧬 Advanced Molecule Analyzer", 
        page_icon="🔬", 
        layout="wide"
    )
    
    # Title and Description
    st.title("🔬 ProSpectra: Comprehensive Molecular Analysis ")
    st.markdown("""## Explore Molecular Structures with Advanced Analytics
    
    This tool provides in-depth analysis of molecular structures using SMILES notation. 
    Discover structural properties, pharmacological insights, and visualization.
    """)
    
    # Sidebar for Input and Examples
    st.sidebar.header("🧪 Molecule Input")
    
    # SMILES Input with Examples
    input_mode = st.sidebar.radio(
        "Input Method", 
        ["Manual Entry", "Predefined Examples"]
    )
    
    # Example Molecules
    example_molecules = {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",
        "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    }
    
    if input_mode == "Manual Entry":
        smiles = st.sidebar.text_input(
            "Enter SMILES String:", 
            value="CCNc1nc2ccc1n2-c1nc2c(cc1O)N2C1=C(C)N=C1c1nc2c(O)cc1c(=O)n2C1=NC2=C(N)C=C1C2"
        )
    else:
        smiles = st.sidebar.selectbox(
            "Select a Predefined Molecule", 
            list(example_molecules.keys())
        )
        smiles = example_molecules[smiles]
    
    # Analyze Button
    if st.sidebar.button("Analyze Molecule"):
        # Error handling for invalid SMILES
        try:
            # Initialize Analyzer
            analyzer = AdvancedMoleculeAnalyzer(smiles)
            
            if analyzer.mol:
                # Create Tabs for Different Views
                tab1, tab2, tab3, tab4 = st.tabs([
                    "🖼️ Molecular Visualizations",
                    "🔍 Molecular Overview", 
                    "📊 Detailed Properties", 
                    "🧪 Pharmacological Insights",
                ])
                
                with tab2:
                    # Overview Section
                    st.header("Molecular Overview")
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader("Structural Information")
                        overview_props = [
                            "Molecular Formula", 
                            "Molecular Weight", 
                            "Total Atom Count", 
                            "Total Bond Count"
                        ]
                        for prop in overview_props:
                            st.metric(prop, analyzer.descriptors.get(prop, "N/A"))
                    
                    with col2:
                        st.subheader("Key Chemical Properties")
                        chem_props = [
                            "LogP", 
                            "H-Bond Donors", 
                            "H-Bond Acceptors", 
                            "Rotatable Bonds"
                        ]
                        for prop in chem_props:
                            st.metric(prop, analyzer.descriptors.get(prop, "N/A"))
                
                with tab3:
                    # Detailed Properties
                    st.header("Comprehensive Molecular Descriptors")
                    
                    # Create DataFrame for detailed properties
                    detailed_df = pd.DataFrame.from_dict(
                        analyzer.descriptors, 
                        orient='index', 
                        columns=['Value']
                    )
                    detailed_df.index.name = 'Descriptor'
                    
                    # Interactive Table
                    st.dataframe(detailed_df, use_container_width=True)
                    
                    # Visualize selected properties
                    prop_choices = st.multiselect(
                        "Choose Properties to Visualize", 
                        detailed_df.index.tolist(),
                        default=['Molecular Weight', 'LogP', 'H-Bond Donors', 'H-Bond Acceptors']
                    )
                    
                    # Bar chart of selected properties
                    if prop_choices:
                        fig = px.bar(
                            detailed_df.loc[prop_choices], 
                            x=detailed_df.loc[prop_choices].index, 
                            y='Value',
                            title="Selected Molecular Properties"
                        )
                        st.plotly_chart(fig, use_container_width=True)
                
                with tab4:
                    # Pharmacological Insights
                    st.header("Pharmacological Assessment")
                    
                    # Pharmacological Properties
                    pharm_props = analyzer.analyze_pharmacological_properties()
                    
                    # Display as metrics
                    cols = st.columns(2)
                    for i, (prop, value) in enumerate(pharm_props.items()):
                        with cols[i % 2]:
                            st.metric(prop, value)
                    
                    # Lipinski Rule Visualization
                    st.subheader("Lipinski's Rule of Five")
                    lipinski_violations = analyzer.descriptors.get("Lipinski Violations", 0)
                    
                    violation_color = "green" if lipinski_violations == 0 else "red" if lipinski_violations > 2 else "orange"
                    st.markdown(f"""
                    <div style="background-color:{violation_color}; padding:10px; border-radius:5px;">
                    <h3 style="color:white;">Lipinski Violations: {lipinski_violations}</h3>
                    </div>
                    """, unsafe_allow_html=True)
                
                with tab1:
                    # Molecular Visualizations
                    st.header("Molecular Visualizations")
                    
                    visualizations = analyzer.generate_molecule_visualization()
                    
                    cols = st.columns(2)
                    for (name, img_base64), col in zip(visualizations.items(), cols):
                        with col:
                            st.subheader(name)
                            st.image(f"data:image/png;base64,{img_base64}", use_column_width=True)
                    if analyzer.mol:
                        add_3d_visualization(st, analyzer)
            else:
                st.error("Invalid SMILES string. Unable to process the molecule.")
        except Exception as e:
            st.error(f"An error occurred: {str(e)}")
    
    # Footer
    st.sidebar.markdown("---")
    st.sidebar.info(
        "🧬 Advanced Molecule Analyzer\n"
        "Powered by RDKit and Streamlit"
    )

if __name__ == "__main__":
    main()