from flask import Flask, request, jsonify, send_file
from rdkit.Chem import AllChem
import py3Dmol
import pandas as pd
from chembl_webresource_client.new_client import new_client
import pickle
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import numpy as np
from PIL import Image
import subprocess
import os
import base64
from io import BytesIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd
from groq import Groq
from transformers import pipeline
from flask_cors import CORS
import json
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import io
import tempfile

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Function to search for targets
def search_target(targetname):
    target = new_client.target
    target_query = target.search(targetname)
    targets = pd.DataFrame.from_dict(target_query)
    return targets

# Function to get bioactivity data
def get_bioactivity_data(target_index, targets):
    selected_target = targets.target_chembl_id[target_index]
    activity = new_client.activity
    res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
    df = pd.DataFrame.from_dict(res)
    df.to_csv('./results/bioactivity_data_1.csv', index=False)
    return df, selected_target

# Function to process bioactivity data
def process_bioactivity_data(df):
    df2 = df[df.standard_value.notna()]
    df2 = df2[df.canonical_smiles.notna()]
    df2_nr = df2.drop_duplicates(['canonical_smiles'])
    selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
    df3 = df2_nr[selection]
    df3.to_csv('./results/bioactivity_data_2.csv', index=False)
    df4 = pd.read_csv('./results/bioactivity_data_2.csv')

    bioactivity_threshold = []
    for i in df4.standard_value:
      if float(i) >= 10000:
        bioactivity_threshold.append("inactive")
      elif float(i) <= 1000:
        bioactivity_threshold.append("active")
      else:
        bioactivity_threshold.append("intermediate")

    bioactivity_class = pd.Series(bioactivity_threshold, name='bioactivity_class')
    df5 = pd.concat([df4, bioactivity_class], axis=1)
    return df5

def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1,1)
    i = 0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])

        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1

    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)

    return descriptors

# Function to calculate Lipinski descriptors
def calculate_lipinski(df5):
    df5.to_csv('./results/bioactivity_data_3.csv', index=False)
    df = pd.read_csv('./results/bioactivity_data_3.csv')
    df_lipinski = lipinski(df.canonical_smiles)
    df_combined = pd.concat([df, df_lipinski], axis=1)
    return df_combined

def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', axis=1)

    return x

def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', axis=1)

    return x

# Function to calculate pIC50
def calculate_pIC50(df_combined):
    df_norm = norm_value(df_combined)
    df_final = pIC50(df_norm)
    # Removing NaN/Infinite Values
    df_final = df_final.dropna()
    df_final = df_final[~df_final.isin([np.inf, -np.inf]).any(axis=1)]
    df_2class = df_final[df_final['bioactivity_class'] != 'intermediate']
    return df_2class

# Function to create plots
def create_plots(df_2class):
    sns.set(style='ticks')
    plt.figure(figsize=(5.5, 5.5))

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    sns.countplot(x='bioactivity_class', data=df_2class, ax=ax1, edgecolor='black')
    ax1.set_xlabel('Bioactivity class', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=14, fontweight='bold')
    ax1.set_title('Bioactivity Class Distribution')
    
    sns.scatterplot(x='MW', y='LogP', data=df_2class, hue='bioactivity_class', size='pIC50', ax=ax2, edgecolor='black', alpha=0.7)
    ax2.set_xlabel('MW', fontsize=14, fontweight='bold')
    ax2.set_ylabel('LogP', fontsize=14, fontweight='bold')
    ax2.set_title('MW vs LogP')
    
    sns.boxplot(x='bioactivity_class', y='pIC50', data=df_2class, ax=ax3)
    ax3.set_xlabel('Bioactivity class', fontsize=14, fontweight='bold')
    ax3.set_ylabel('pIC50 value', fontsize=14, fontweight='bold')
    ax3.set_title('pIC50 Distribution by Bioactivity Class')
    
    # Save plot to memory
    img_bytes = io.BytesIO()
    plt.tight_layout()
    fig.savefig(img_bytes, format='png')
    img_bytes.seek(0)
    plt.close(fig)
    
    return img_bytes

from sklearn.feature_selection import VarianceThreshold

def remove_low_variance(input_data, threshold=0.1):
    selection = VarianceThreshold(threshold)
    selection.fit(input_data)
    return input_data[input_data.columns[selection.get_support(indices=True)]]

def run_model(df_2class):
    # Save the input data
    df_2class.to_csv('./results/bioactivity_data_5.csv', index=False)
    df3 = pd.read_csv('./results/bioactivity_data_5.csv')
    
    # Prepare the SMILES file for PaDEL-Descriptor
    selection = ['canonical_smiles', 'molecule_chembl_id']
    df3_selection = df3[selection]
    df3_selection.to_csv('./results/molecule.smi', sep='\t', index=False, header=False)
    
    # Run PaDEL-Descriptor
    subprocess.run(["bash", "padel.sh"])
    
    # Load the descriptors
    df3_X = pd.read_csv('./results/descriptors_output.csv')
    df3_Y = df3['pIC50']
    
    # Combine descriptors and target values
    dataset3 = pd.concat([df3_X, df3_Y], axis=1)
    dataset3.to_csv('./results/GENAIdataset.csv', index=False)
    
    # Drop the 'Name' column from descriptors
    df3_X = df3_X.drop(columns=['Name'])
    dataset3 = pd.concat([df3_X, df3_Y], axis=1)
    dataset3.to_csv('./results/bioactivity_data_6_fingerprints.csv', index=False)
    
    # Load the final dataset
    dataset = pd.read_csv('./results/bioactivity_data_6_fingerprints.csv')
    
    # Drop rows where the target column (pIC50) has NaN values
    dataset = dataset.dropna(subset=['pIC50'])
    
    # Split into features (X) and target (Y)
    X = dataset.drop(['pIC50'], axis=1)
    Y = dataset['pIC50']
    
    # Remove low variance features
    X = remove_low_variance(X, threshold=0.1)
    X.to_csv('./results/descriptor_list.csv', index=False)
    
    # Train the model
    model = RandomForestRegressor(n_estimators=500, random_state=42)
    model.fit(X, Y)
    
    # Make predictions
    Y_pred = model.predict(X)
    mse = mean_squared_error(Y, Y_pred)
    r2 = r2_score(Y, Y_pred)
    
    # Create scatter plot for predictions
    plt.figure(figsize=(5, 5))
    plt.scatter(x=Y, y=Y_pred, c="#7CAE00", alpha=0.3)
    z = np.polyfit(Y, Y_pred, 1)
    p = np.poly1d(z)
    plt.plot(Y, p(Y), "#F8766D")
    plt.ylabel('Predicted pIC50')
    plt.xlabel('Experimental pIC50')
    
    # Save scatter plot to memory
    scatter_bytes = io.BytesIO()
    plt.savefig(scatter_bytes, format='png')
    scatter_bytes.seek(0)
    plt.close()
    
    # Extract feature importances
    importances = model.feature_importances_
    feature_importances_df = pd.DataFrame({'Feature': list(range(X.shape[1])), 'Importance': importances})
    feature_importances_df = feature_importances_df.sort_values('Importance', ascending=False)
    
    N = 15  # Number of top features to display
    top_features = feature_importances_df.head(N)
    
    # Load mapping data
    try:
        mapping_csv_path = "./results/csvCHEM.csv"
        mapping_df = pd.read_csv(mapping_csv_path)
        mapping_df.rename(columns={'Bit Position': 'Feature'}, inplace=True)
        top_features_with_mapping = top_features.merge(mapping_df, on='Feature', how='left')
    except:
        top_features_with_mapping = top_features
    
    # Store the first column (feature names) in an array
    top_feature_indices = top_features['Feature'].values
    
    # Generate SMILES using Groq and ChemGPT
    api_key = "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"  # Note: API keys should be stored securely
    client = Groq(api_key=api_key)
    
    results = {
        'model_stats': {
            'mse': mse,
            'r2': r2
        },
        'top_features': top_features_with_mapping.to_dict(orient='records'),
        'generated_molecules': []
    }
    
    try:
        # Load CSV file containing chemical substructures
        csv_path = "./results/csvCHEM.csv"
        csvChem = pd.read_csv(csv_path)
        L1 = top_feature_indices

        def get_substructure_description(substructure):
            # Prepare messages for Groq API
            messages = [
                {
                    "role": "system",
                    "content": "Given the input formula, return the description of the chemical element's presence with a specified number of atoms in one line. For example, for the input '>= 4 H', output should be 'presence of 4 hydrogen atoms'"
                },
                {
                    "role": "user", 
                    "content": substructure
                }
            ]

            # Create a chat completion using the Groq API
            try:
                chat_completion = client.chat.completions.create(
                    messages=messages,
                    model="llama3-8b-8192",
                    temperature=0.5,
                    max_tokens=1024,
                    top_p=1,
                    stop=None,
                    stream=False,
                )
                
                return chat_completion.choices[0].message.content.strip()
            except Exception as e:
                print(f"Error processing substructure {substructure}: {e}")
                return f"Could not process {substructure}"

        # Collect descriptions for the specified indices
        descriptions = []
        for i in L1:
            substructure = csvChem.loc[i, "Bit Substructure"]
            description = get_substructure_description(substructure)
            descriptions.append(f"{description}")

        # Combine all descriptions into a single string
        final_description_string = " | ".join(descriptions)

        # Generate a molecule description
        fingerprint_description = (
            f"Generate a molecule with the following features: {final_description_string}"
        )

        # Load the ChemGPT pipeline from Hugging Face Transformers
        pipe = pipeline("text-generation", model="ncfrey/ChemGPT-4.7M")

        # Generate a compound based on the fingerprint description
        generated_output = pipe(
            fingerprint_description,
            max_length=25,  # Limit the length of the generated SMILES
            num_return_sequences=1,  # Generate one compound
            do_sample=True,  # Enable randomness for diversity
            truncation=True  # Enable truncation if the text is too long
        )

        generated_smiles = [None] * len(generated_output)
        for i in range(0, len(generated_output)):
            generated_smiles[i] = generated_output[i]['generated_text']

        generated_smiles2 = [None] * len(generated_output)
        for i in range(len(generated_output)):
            smiles = generated_smiles[i]
            if "[" in smiles:
                k = smiles.index("[")  # Find the first occurrence of '['
                generated_smiles2[i] = smiles[k:]
            else:
                generated_smiles2[i] = smiles  # If '[' not found, keep the original string

        # Process generated SMILES with Groq
        for i in range(len(generated_output)):
            chat_completion = client.chat.completions.create(
                messages=[
                    {
                        "role": "system",
                        "content": "Just Give me the SMILES string(Simplified Molecular Input Line Entry System) of the following SMILES with proper bond. I want to use the response as input to another software, hence don't print ANY EXTRA message, just the SMILES string Strictly",
                    },
                    {
                        "role": "user",
                        "content": generated_smiles2[i],
                    }
                ],
                model="llama3-8b-8192",
                temperature=0.5,
                max_tokens=1024,
                top_p=1,
                stop=None,
                stream=False,
            )
            results['generated_molecules'].append({
                'description': f"Generated molecule {i+1}",
                'smiles': chat_completion.choices[0].message.content
            })
    except Exception as e:
        results['error'] = str(e)
    
    return dataset, results, scatter_bytes

# API Routes
@app.route('/api/search_target', methods=['POST'])
def api_search_target():
    data = request.json
    targetname = data.get('targetname')
    if not targetname:
        return jsonify({'error': 'Target name is required'}), 400
    
    try:
        targets = search_target(targetname)
        return jsonify({'targets': targets.to_dict(orient='records')})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/process_target', methods=['POST'])
def api_process_target():
    data = request.json
    targetname = data.get('targetname')
    target_index = data.get('target_index', 0)
    
    if not targetname:
        return jsonify({'error': 'Target name is required'}), 400
    
    try:
        # Search for targets
        targets = search_target(targetname)
        
        if target_index >= len(targets):
            return jsonify({'error': 'Invalid target index'}), 400
        
        # Get bioactivity data
        bioactivity_data, selected_target = get_bioactivity_data(target_index, targets)
        
        # Process data
        processed_data = process_bioactivity_data(bioactivity_data)
        df_combined = calculate_lipinski(processed_data)
        df_2class = calculate_pIC50(df_combined)
        
        # Create plots
        plot_bytes = create_plots(df_2class)
        plot_base64 = base64.b64encode(plot_bytes.getvalue()).decode('utf-8')
        
        # Run model
        dataset, model_results, scatter_bytes = run_model(df_2class)
        scatter_base64 = base64.b64encode(scatter_bytes.getvalue()).decode('utf-8')
        
        # Prepare response
        result = {
            'selected_target': selected_target,
            'bioactivity_data': processed_data.to_dict(orient='records'),
            'plot_data': plot_base64,
            'scatter_plot': scatter_base64,
            'dataset': dataset.to_dict(orient='records'),  # Limit to 20 rows for efficiency
            'model_results': model_results
        }
        
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/download_data', methods=['GET'])
def download_data():
    try:
        return send_file('./results/bioactivity_data_3.csv', 
                         mimetype='text/csv',
                         as_attachment=True,
                         download_name='bioactivity_data.csv')
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5004)