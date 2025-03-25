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