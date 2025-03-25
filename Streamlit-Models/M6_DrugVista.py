import streamlit as st
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
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer

# Set page configuration
st.set_page_config(layout="wide")
st.title("BioActivity Predictor")

# Sidebar
st.sidebar.title('Input Parameters')
targetname = st.sidebar.text_input("Enter Target Name:")
target_index = st.sidebar.number_input("Enter Target Index:", min_value=0, value=0, step=1)
with st.sidebar.header('Input CSV Data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
""")

# Function to search for targets
def search_target(targetname):
    target = new_client.target
    target_query = target.search(targetname)
    targets = pd.DataFrame.from_dict(target_query)
    return targets

# Function to get bioactivity data
def get_bioactivity_data(target_index, targets):
    selected_target = targets.target_chembl_id[target_index]
    st.write(f"Selected Target: {selected_target}")

    activity = new_client.activity
    res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
    df = pd.DataFrame.from_dict(res)
    df.to_csv('bioactivity_data_1.csv', index=False)
    return df

# Function to process bioactivity data
def process_bioactivity_data(df):
    df2 = df[df.standard_value.notna()]
    df2 = df2[df.canonical_smiles.notna()]
    df2_nr = df2.drop_duplicates(['canonical_smiles'])
    selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
    df3 = df2_nr[selection]
    df3.to_csv('bioactivity_data_2.csv', index=False)
    df4 = pd.read_csv('bioactivity_data_2.csv')

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

# Function to calculate Lipinski descriptors
def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
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

        if i == 0:
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i + 1

    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)
    return descriptors

def calculate_lipinski(df5):
    df5.to_csv('bioactivity_data_3.csv', index=False)
    df = pd.read_csv('bioactivity_data_3.csv')
    df_lipinski = lipinski(df.canonical_smiles)
    df_combined = pd.concat([df, df_lipinski], axis=1)
    return df_combined

def pIC50(input):
    pIC50 = []
    for i in input['standard_value_norm']:
        molar = i * (10**-9)  # Converts nM to M
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
    
    return fig

# Function to remove low variance features
def remove_low_variance(input_data, threshold=0.1):
    selection = VarianceThreshold(threshold)
    selection.fit(input_data)
    return input_data[input_data.columns[selection.get_support(indices=True)]]

# Function to train the model
def model(df_2class):
    df_2class.to_csv('bioactivity_data_5.csv')
    df3 = pd.read_csv('bioactivity_data_5.csv')
    selection = ['canonical_smiles', 'molecule_chembl_id']
    df3_selection = df3[selection]
    df3_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)
    subprocess.run(["bash", "padel.sh"])
    df3_X = pd.read_csv('descriptors_output.csv')

    # GENAI Dataset
    df3_Y = df3['pIC50']
    dataset3 = pd.concat([df3_X, df3_Y], axis=1)
    dataset3.to_csv('GENAIdataset.csv', index=False)
    df3_X = df3_X.drop(columns=['Name'])
    dataset3 = pd.concat([df3_X, df3_Y], axis=1)
    dataset3.to_csv('bioactivity_data_6_fingerprints.csv', index=False)
    dataset = pd.read_csv('bioactivity_data_6_fingerprints.csv')
    st.dataframe(dataset)

    # Drop rows where the target column (pIC50) has NaN values
    dataset = dataset.dropna(subset=['pIC50'])

    # Split into features (X) and target (Y)
    X = dataset.drop(['pIC50'], axis=1)
    Y = dataset['pIC50']

    # Remove low variance features
    X = remove_low_variance(X, threshold=0.1)
    X.to_csv('descriptor_list.csv', index=False)

    # Train the model
    model = RandomForestRegressor(n_estimators=500, random_state=42)
    model.fit(X, Y)

    # Evaluate the model
    r2 = model.score(X, Y)
    Y_pred = model.predict(X)
    r1 = mean_squared_error(Y, Y_pred)
    r2 = r2_score(Y, Y_pred)

    # Plot the results
    plt.figure(figsize=(5, 5))
    plt.scatter(x=Y, y=Y_pred, c="#7CAE00", alpha=0.3)
    z = np.polyfit(Y, Y_pred, 1)
    p = np.poly1d(z)
    plt.plot(Y, p(Y), "#F8766D")
    plt.ylabel('Predicted pIC50')
    plt.xlabel('Experimental pIC50')
    st.pyplot(plt)
    plt.savefig("scatter_plot.pdf", format="pdf")

    # Feature importance
    importances = model.feature_importances_
    feature_importances_df = pd.DataFrame({'Feature': list(range(X.shape[1])), 'Importance': importances})
    feature_importances_df = feature_importances_df.sort_values('Importance', ascending=False)

    N = 15  # Number of top features to display
    top_features = feature_importances_df.head(N)
    mapping_csv_path = "csvCHEM.csv"
    mapping_df = pd.read_csv(mapping_csv_path)
    mapping_df.rename(columns={'Bit Position': 'Feature'}, inplace=True)
    top_features_with_mapping = top_features.merge(mapping_df, on='Feature', how='left')

    st.write(f"Top {N} most important fingerprint features with substructures:")
    st.dataframe(top_features_with_mapping)

# Main app logic
if targetname:
    targets = search_target(targetname)
    st.write("Found Targets:")
    st.dataframe(targets)
    
    if target_index < len(targets):
        bioactivity_data = get_bioactivity_data(target_index, targets)
        processed_data = process_bioactivity_data(bioactivity_data)
        df_combined = calculate_lipinski(processed_data)
        df_2class = calculate_pIC50(df_combined)

        # Mean Imputation for Numerical Columns
        numerical_cols = df_2class.select_dtypes(include=['float64', 'int64']).columns
        numerical_imputer = SimpleImputer(strategy='mean')
        df_2class[numerical_cols] = numerical_imputer.fit_transform(df_2class[numerical_cols])

        fig = create_plots(df_2class)
        st.pyplot(fig)

        model_output = model(df_2class)
        st.write(model_output)

        # Download buttons for data
        csv = processed_data.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="bioactivity_data.csv">Download Bioactivity Data</a>'
        st.markdown(href, unsafe_allow_html=True)
    else:
        st.error("Invalid target index. Please select a valid index.")
else:
    st.info("Please enter a target name in the sidebar to begin.")