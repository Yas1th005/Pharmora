import os
import re
import pandas as pd
import numpy as np
from typing import List, Dict, Any
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
import tempfile
import json

# LangChain Imports
from langchain_google_genai import GoogleGenerativeAI, GoogleGenerativeAIEmbeddings
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate
from langchain.vectorstores import FAISS
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.docstore.document import Document
from langchain.memory import ConversationBufferMemory
from langchain_google_genai import ChatGoogleGenerativeAI

class DrugDesignRAG:

        
    #     # Set up chains
    #     self.setup_chains()
    def __init__(self, api_key=None):
    # Set API Key
        if api_key:
            os.environ["GOOGLE_API_KEY"] = api_key
        elif "GOOGLE_API_KEY" not in os.environ:
            raise ValueError("Please provide a Google API key")
        
        # Initialize LangChain components
        # Fix: Use ChatGoogleGenerativeAI instead of GoogleGenerativeAI
        self.llm = ChatGoogleGenerativeAI(model="models/gemini-1.5-pro", temperature=0.7)
        

        # Fix: Ensure embedding model name is correct
        #self.embeddings = GoogleGenerativeAIEmbeddings(model="embedding-001")
        self.embeddings = GoogleGenerativeAIEmbeddings(
    model="models/text-embedding-004"  # Use full model name with prefix
)
        
        # Set up vector store - we'll load and index chemistry knowledge
        self.setup_vector_store()
        
        # Set up prompt templates
        self.setup_prompts()
        
        # Set up chains
        self.setup_chains()

    
    def setup_vector_store(self):
        """
        Set up the vector store with chemistry knowledge.
        In a real implementation, you would load more comprehensive data.
        """
        try:
            # First try to load existing vector store
            self.vector_store = FAISS.load_local("faiss_index", self.embeddings)
            print("Loaded existing vector store")
        except Exception as e:
            print(f"Creating new vector store... ({str(e)})")
            
            chemistry_knowledge = [
               
                
    {
        "name": "Pyridine",
        "description": "A six-membered aromatic heterocycle with a nitrogen atom.",
        "smiles": "C1=CC=NC=C1",
        "molecular_weight": 79.1,
        "logP": 0.65,
        "bioactivity": ["Antibacterial", "Enzyme Inhibitor"]
    },
    {
        "name": "Pyrrole",
        "description": "A five-membered aromatic heterocycle with a nitrogen atom.",
        "smiles": "C1=CC=NC1",
        "molecular_weight": 67.09,
        "logP": 0.34,
        "bioactivity": ["Antifungal", "Antiviral"]
    },
    {
        "name": "Morpholine",
        "description": "A six-membered saturated heterocycle containing an oxygen and a nitrogen atom.",
        "smiles": "C1COCCN1",
        "molecular_weight": 87.1,
        "logP": -0.42,
        "bioactivity": ["Anticancer", "Cytotoxic"]
    },
    {
        "name": "Piperazine",
        "description": "A six-membered saturated heterocycle containing two nitrogen atoms.",
        "smiles": "N1CCNCC1",
        "molecular_weight": 86.14,
        "logP": -1.23,
        "bioactivity": ["Antipsychotic", "Antidepressant"]
    },
    {
        "name": "Thiazole",
        "description": "A five-membered aromatic heterocycle containing a sulfur and a nitrogen atom.",
        "smiles": "C1=CSC=N1",
        "molecular_weight": 85.11,
        "logP": 0.37,
        "bioactivity": ["Antibacterial", "Anti-inflammatory"]
    },
    {
        "name": "Quinoline",
        "description": "A bicyclic heterocycle consisting of a benzene ring fused to a pyridine ring.",
        "smiles": "C1=CC=C2C=CC=NC2=C1",
        "molecular_weight": 129.16,
        "logP": 2.0,
        "bioactivity": ["Antimalarial", "Antitumor"]
    },
    {
        "name": "Indole",
        "description": "A bicyclic heterocycle consisting of a benzene ring fused to a pyrrole ring.",
        "smiles": "C1=CC2=CC=CN2C=C1",
        "molecular_weight": 117.15,
        "logP": 1.98,
        "bioactivity": ["Anticancer", "Neuroprotective"]
    },
    {
        "name": "Furan",
        "description": "A five-membered aromatic heterocycle containing an oxygen atom.",
        "smiles": "C1=COC=C1",
        "molecular_weight": 68.07,
        "logP": 0.38,
        "bioactivity": ["Antifungal", "Antibacterial"]
    },
    {
        "name": "Benzimidazole",
        "description": "A bicyclic heterocycle consisting of a benzene ring fused to an imidazole ring.",
        "smiles": "C1=CC2=NC=NC=C2C=C1",
        "molecular_weight": 118.13,
        "logP": 1.78,
        "bioactivity": ["Anthelmintic", "Anticancer"]
    },
    {
        "name": "Isoquinoline",
        "description": "A bicyclic heterocycle consisting of a benzene ring fused to a pyridine ring in a different orientation than quinoline.",
        "smiles": "C1=CC=C2C=NC=CC2=C1",
        "molecular_weight": 129.16,
        "logP": 2.09,
        "bioactivity": ["Antihypertensive", "Antitumor"]
    },
    {
        "name": "Indole",
        "description": "A bicyclic heterocycle consisting of a benzene ring fused to a pyrrole ring.",
        "smiles": "C1=CC=C2C=CC=CN2C1",
        "molecular_weight": 117.15,
        "logP": 1.98,
        "bioactivity": ["Anticancer", "Serotonin Receptor Modulator"]
    },
    {
        "name": "Benzofuran",
        "description": "A bicyclic heterocycle consisting of a benzene ring fused to a furan ring.",
        "smiles": "C1=CC2=C(C=CO2)C=C1",
        "molecular_weight": 118.13,
        "logP": 2.34,
        "bioactivity": ["Antifungal", "Antiviral"]
    },
    {
        "name": "Benzothiophene",
        "description": "A bicyclic heterocycle consisting of a benzene ring fused to a thiophene ring.",
        "smiles": "C1=CC2=C(C=CS2)C=C1",
        "molecular_weight": 134.19,
        "logP": 2.76,
        "bioactivity": ["Antitumor", "Antioxidant"]
    },
    
    {
        "name": "Lipinski's Rule of Five",
        "description": "Drug-like molecules generally have: molecular weight ≤ 500 Daltons, log P ≤ 5, hydrogen bond donors ≤ 5, hydrogen bond acceptors ≤ 10.",
        "relevance": "Guides early-stage drug discovery and screening."
    },
    {
        "name": "ADMET Properties",
        "description": "Absorption, Distribution, Metabolism, Excretion, and Toxicity — crucial factors determining drug viability.",
        "relevance": "Helps optimize pharmacokinetics and safety profiles."
    },
    {
        "name": "Heterocycles in Drug Discovery",
        "description": "Heterocycles are commonly found in drugs due to their ability to interact with biological targets.",
        "relevance": "Used as core scaffolds in many therapeutics."
    },
    {
        "name": "Ring Systems in Drugs",
        "description": "Provide scaffolds for positioning functional groups in three-dimensional space, influencing binding affinity.",
        "relevance": "Crucial for molecular docking and receptor-ligand interactions."
    },
    {
        "name": "Unsaturated Rings",
        "description": "Contain double bonds, introducing rigidity and planarity into a molecule.",
        "relevance": "Affects molecular geometry and electronic distribution."
    },
    {
        "name": "Heteroatoms in Drug Molecules",
        "description": "Often serve as hydrogen bond donors or acceptors, enhancing binding interactions.",
        "relevance": "Essential for protein-ligand recognition."
    },
    {
        "name": "Saturated Rings",
        "description": "Provide conformational flexibility compared to their unsaturated counterparts.",
        "relevance": "Influences the dynamic behavior of molecules in solution."
    }
]

            


            documents = []
            for entry in chemistry_knowledge:
                # Serialize the dictionary to a JSON string or format it as needed
                page_content = json.dumps(entry, indent=2)  # JSON string
                # Alternatively, you can format it as a readable string:
                # page_content = f"Name: {entry['name']}\nDescription: {entry['description']}\nSMILES: {entry['smiles']}\nMolecular Weight: {entry['molecular_weight']}\nLogP: {entry['logP']}\nBioactivity: {', '.join(entry['bioactivity'])}"
                
                # Create a Document object
                doc = Document(
                    page_content=page_content,
                    metadata={"source": "chemistry_knowledge"}
                )
                documents.append(doc)
            
            # Create text chunks to ensure they're not too large
            text_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=0)
            docs = text_splitter.split_documents(documents)
            
            # Create vector store
            self.vector_store = FAISS.from_documents(docs, self.embeddings)
            
            # Create the directory if it doesn't exist
            os.makedirs("faiss_index", exist_ok=True)
            
            # Save vector store for future use
            self.vector_store.save_local("faiss_index")
    
    def setup_prompts(self):
        """Set up prompt templates for different chains"""
        #Main drug design prompt
        self.design_prompt = PromptTemplate(
            input_variables=["substructures", "context"],
            template="""
            You are an expert medicinal chemist. Based on the following important substructures and their 
            importance scores, design a novel drug compound that incorporates these features.
            
            # Important Substructures:
            {substructures}
            
            # Relevant Chemistry Knowledge:
            {context}
            
            Please create a novel drug compound with the following specifications:
            1. Incorporate as many of the high-importance substructures as possible
            2. Ensure the compound meets Lipinski's Rule of Five
            3. Optimize for drug-likeness and synthetic feasibility
            4. Avoid known toxicophores
            
            For the output, provide:
            1. A name for the compound
            2. The SMILES notation (valid and correct chemical structure)
            
            
            
            ONLY Be precise with the SMILES notation - it must represent a valid chemical structure, without any additional texts.Dont give me any additional texts or characters.
            """
        )
        
        # Validation prompt
        self.validation_prompt = PromptTemplate(
            input_variables=["smiles", "substructures", "errors"],
            template="""
            You are an expert computational chemist. The following SMILES string was generated to represent a drug
            compound incorporating these substructures:
            
            # Substructures:
            {substructures}
            
            # SMILES:
            {smiles}
            
            # Issues detected:
            {errors}
            
            Please correct the SMILES string to create a valid chemical structure that incorporates the required substructures. Ensure the corrected structure is chemically valid and plausible.
            
            Return only the corrected SMILES string without any additional text.
            """
        )
   
    def setup_chains(self):
        """Set up LangChain chains"""
        # Fix: Ensure chains are properly initialized
        self.design_chain = LLMChain(
            llm=self.llm,
            prompt=self.design_prompt,
            verbose=True
        )
        
        self.validation_chain = LLMChain(
            llm=self.llm,
            prompt=self.validation_prompt,
            verbose=True
        )
        
    def parse_substructures(self, input_text):
        """Parse substructure information from input text"""
        print("Parsing substructure data...")
        lines = input_text.strip().split('\n')
        
        # Skip header if present
        if 'importance' in lines[0].lower() and 'substructure' in lines[0].lower():
            lines = lines[1:]
        
        substructures = []
        for line in lines:
            try:
                # Extract importance value and substructure description
                match = re.match(r'([\d\.]+)\s+(.*)', line.strip())
                if match:
                    importance = float(match.group(1))
                    description = match.group(2)
                    substructures.append({
                        'importance': importance,
                        'description': description
                    })
            except Exception as e:
                print(f"Error parsing line: {line}, {e}")
        
        # Sort by importance (descending)
        substructures = sorted(substructures, key=lambda x: x['importance'], reverse=True)
        return substructures
    
    def get_relevant_context(self, substructures):
        """Retrieve relevant knowledge from the vector store"""
        # Create query from substructures
        query = " ".join([s["description"] for s in substructures[:5]])
        
        # Search vector store
        try:
            docs = self.vector_store.similarity_search(query, k=5)
            
            # Extract and return relevant context
            context = "\n".join([doc.page_content for doc in docs])
            return context
        except Exception as e:
            print(f"Error retrieving context: {e}")
            # Return default context if search fails
            return "Heterocycles are commonly found in drugs. Lipinski's Rule of Five is important for drug-likeness."
    
    def extract_smiles(self, generated_text):
        """Extract SMILES notation from generated text"""
        # Try multiple regex patterns to extract SMILES
        patterns = [
            r'SMILES[:\s]*([^\n]+)',
            r'SMILES notation[:\s]*([^\n]+)',
            r'SMILES string[:\s]*([^\n]+)',
            r'notation[:\s]*([A-Za-z0-9\[\]\(\)\.\=\#\-\+\@\\\/@\*\%\$\{\}]+)'
        ]
        
        for pattern in patterns:
            matches = re.search(pattern, generated_text, re.IGNORECASE)
            if matches:
                return matches.group(1).strip()
        
        # If no matches, look for chemical notation patterns
        chemical_patterns = [
            r'([A-Za-z0-9\[\]\(\)\.\=\#\-\+\@\\\/@\*\%\$\{\}]{10,})'
        ]
        
        for pattern in chemical_patterns:
            matches = re.findall(pattern, generated_text)
            if matches:
                # Return the longest match as it's most likely to be a full SMILES string
                return max(matches, key=len)
        
        return None
    
    def validate_structure(self, smiles):
        """Validate chemical structure and return RDKit mol object"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None, "Invalid SMILES notation. Could not generate molecule."
            return mol, None
        except Exception as e:
            return None, f"Error validating structure: {str(e)}"
    
    def analyze_drug(self, mol, smiles):
        """Analyze drug properties"""
        try:
            # Calculate basic properties
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
            
            # Create molecule image
            img = Draw.MolToImage(mol, size=(600, 400))
            
            return {
                "mol": mol,
                "img": img,
                "smiles": smiles,
                "properties": {
                    "mw": mw,
                    "logp": logp,
                    "hbd": hbd,
                    "hba": hba,
                    "rotatable_bonds": rotatable_bonds,
                    "tpsa": tpsa,
                    "lipinski_violations": lipinski_violations
                }
            }
        except Exception as e:
            print(f"Error analyzing drug: {e}")
            return None
    
    def generate_drug(self, substructures_list):
        """Generate drug using the RAG system"""
        # Format substructures for prompt
        substructures_formatted = "\n".join([
            f"- Importance {s['importance']:.4f}: {s['description']}" 
            for s in substructures_list
        ])
        
        # Get relevant context from knowledge base
        context = self.get_relevant_context(substructures_list)
        
        # Generate drug design with error handling
        try:
            response = self.design_chain.run(
                substructures=substructures_formatted, 
                context=context
            )
        except Exception as e:
            print(f"Error generating drug design: {e}")
            return {
                "success": False,
                "response": f"Error in LLM generation: {str(e)}",
                "error": str(e)
            }
        
        # Extract SMILES from generated text
        smiles = self.extract_smiles(response)
        
        if not smiles:
            print("Could not extract SMILES notation from generated text.")
            return {
                "success": False,
                "response": response,
                "error": "Could not extract SMILES notation"
            }
        
        print(f"Extracted SMILES: {smiles}")
        
        # Validate structure
        mol, error = self.validate_structure(smiles)
        
        # If invalid, try to correct it
        if mol is None and error:
            print(f"Invalid structure: {error}. Attempting to correct...")
            try:
                corrected_smiles = self.validation_chain.run(
                    smiles=smiles,
                    substructures=substructures_formatted,
                    errors=error
                )
                
                # Try the corrected SMILES
                mol, error = self.validate_structure(corrected_smiles)
                if mol is not None:
                    smiles = corrected_smiles
                    print(f"Corrected SMILES: {smiles}")
            except Exception as e:
                print(f"Error in validation chain: {e}")
                return {
                    "success": False,
                    "response": response,
                    "error": f"Validation error: {str(e)}",
                    "smiles": smiles
                }
        
        if mol is None:
            return {
                "success": False,
                "response": response,
                "error": error,
                "smiles": smiles
            }
        
        # Analyze drug properties
        analysis = self.analyze_drug(mol, smiles)
        
        if analysis:
            return {
                "success": True,
                "response": response,
                "smiles": smiles,
                "analysis": analysis
            }
        else:
            return {
                "success": False,
                "response": response,
                "error": "Failed to analyze drug properties",
                "smiles": smiles
            }
    
    def save_results(self, result, output_dir="."):
        """Save results to files"""
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Save the full response
        with open(f"{output_dir}/drug_design_response.txt", "w") as f:
            f.write(result["response"])
        
        if result["success"]:
            # Save SMILES
            with open(f"{output_dir}/drug_smiles.txt", "w") as f:
                f.write(result["smiles"])
            
            # Save properties
            properties = result["analysis"]["properties"]
            with open(f"{output_dir}/drug_properties.txt", "w") as f:
                f.write(f"Molecular Weight: {properties['mw']:.2f} g/mol\n")
                f.write(f"LogP: {properties['logp']:.2f}\n")
                f.write(f"Hydrogen Bond Donors: {properties['hbd']}\n")
                f.write(f"Hydrogen Bond Acceptors: {properties['hba']}\n")
                f.write(f"Rotatable Bonds: {properties['rotatable_bonds']}\n")
                f.write(f"Topological Polar Surface Area: {properties['tpsa']:.2f} Å²\n")
                f.write(f"Lipinski Violations: {properties['lipinski_violations']}/4\n")
            
            # Save molecular image
            img = result["analysis"]["img"]
            img.save(f"{output_dir}/drug_structure.png")
            
            # Save a report with everything
            with open(f"{output_dir}/drug_report.md", "w") as f:
                f.write("# Drug Design Report\n\n")
                f.write("## Generated Drug\n\n")
                f.write(result["response"])
                f.write("\n\n## SMILES Notation\n\n")
                f.write(f"`{result['smiles']}`\n\n")
                f.write("## Properties\n\n")
                f.write(f"- Molecular Weight: {properties['mw']:.2f} g/mol\n")
                f.write(f"- LogP: {properties['logp']:.2f}\n")
                f.write(f"- Hydrogen Bond Donors: {properties['hbd']}\n")
                f.write(f"- Hydrogen Bond Acceptors: {properties['hba']}\n")
                f.write(f"- Rotatable Bonds: {properties['rotatable_bonds']}\n")
                f.write(f"- Topological Polar Surface Area: {properties['tpsa']:.2f} Å²\n")
                f.write(f"- Lipinski Violations: {properties['lipinski_violations']}/4\n")
        
        return f"Results saved to {output_dir}/"

def main():
    """Main function to run the drug design RAG system"""
    print("Drug Design RAG System with LangChain")
    print("======================================")
    
    # Get API key
    api_key = input("Enter your Google API key: ")
    
    # Initialize system with error handling
    try:
        system = DrugDesignRAG(api_key=api_key)
    except Exception as e:
        print(f"Error initializing DrugDesignRAG: {e}")
        print("Please check your API key and dependencies.")
        return
    
    # Sample input
    sample_input = """
Importance    Bit substructure
0.0416       >= 5 unsaturated non-aromatic heteroatom-containing ring size 6
0.0303       >= 5 unsaturated non-aromatic heteroatom-containing ring size 5
0.0301       >= 2 saturated or aromatic heteroatom-containing ring size 4
0.0190       >= 1 Bi
0.0184       >= 3 saturated or aromatic nitrogen-containing ring size 6
0.0152       >= 1 Ni


"""
    
    # Get input
    print("\nEnter your substructure data (or press Enter to use sample data):")
    user_input = input()
    
    if not user_input.strip():
        input_text = sample_input
        print("Using sample data:")
        print(sample_input)
    else:
        input_text = user_input
    
    # Parse substructures
    substructures = system.parse_substructures(input_text)
    print(f"Found {len(substructures)} substructures")
    
    # Display substructures
    for i, s in enumerate(substructures):
        print(f"{i+1}. Importance: {s['importance']:.4f} - {s['description']}")
    
    # Generate drug
    print("\nGenerating drug compound...")
    result = system.generate_drug(substructures)
    
    # Display results
    if result["success"]:
        print("\nSuccess! Drug compound generated.")
        print(f"SMILES: {result['smiles']}")
        
        properties = result["analysis"]["properties"]
        print("\nProperties:")
        print(f"- Molecular Weight: {properties['mw']:.2f} g/mol")
        print(f"- LogP: {properties['logp']:.2f}")
        print(f"- Hydrogen Bond Donors: {properties['hbd']}")
        print(f"- Hydrogen Bond Acceptors: {properties['hba']}")
        print(f"- Rotatable Bonds: {properties['rotatable_bonds']}")
        print(f"- Topological Polar Surface Area: {properties['tpsa']:.2f} Å²")
        print(f"- Lipinski Violations: {properties['lipinski_violations']}/4")
        
        # Save results
        output_dir = "drug_design_output"
        save_message = system.save_results(result, output_dir)
        print(f"\n{save_message}")
        print(f"Molecular structure image saved to {output_dir}/drug_structure.png")
    else:
        print("\nFailed to generate valid drug compound.")
        print(f"Error: {result.get('error', 'Unknown error')}")
        print("\nPartial results:")
        print(f"Generated response available in drug_design_response.txt")
        
        # Save partial results
        output_dir = "drug_design_output"
        system.save_results(result, output_dir)

if __name__ == "__main__":
    main()