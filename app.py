from flask import Flask, request, jsonify, render_template

import numpy as np
import pandas as pd
import pickle
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('index.html')

# Load the trained model
model = pickle.load(open('random_forest_ddi_model.pkl', 'rb'))

# Function to fetch SMILES from PubChem
def get_smiles(drug_name):
    base_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/CanonicalSMILES/JSON"
    response = requests.get(base_url)
    if response.status_code == 200:
        try:
            smiles = response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
            return smiles
        except (KeyError, IndexError):
            return None
    else:
        return None

# Function to generate Morgan fingerprint
def smiles_to_morgan(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        arr = np.zeros((1,))
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return None

@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    drug1 = data.get('drug1')
    drug2 = data.get('drug2')
    
    # Fetch SMILES
    smiles1 = get_smiles(drug1)
    smiles2 = get_smiles(drug2)

    if not smiles1 or not smiles2:
        return jsonify({'error': 'Could not retrieve SMILES for one or both drugs.'}), 400

    # Generate fingerprints
    fp1 = smiles_to_morgan(smiles1)
    fp2 = smiles_to_morgan(smiles2)

    if fp1 is None or fp2 is None:
        return jsonify({'error': 'Invalid SMILES or unable to generate fingerprints.'}), 400

    # For now, just using one fingerprint (or later we can combine features)
    input_features = np.abs(fp1 - fp2).reshape(1, -1)

    # Predict
    prediction = model.predict(input_features)[0]
    probability = model.predict_proba(input_features)[0][prediction]

    result = {
        'interaction': 'Interaction Predicted' if prediction == 1 else 'No Major Interaction',
        'confidence': round(float(probability), 3)
    }

    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True)