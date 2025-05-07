import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
import numpy as np

# Load the SMILES data
df = pd.read_csv('drug_smiles.csv')

# Function to convert SMILES to Morgan Fingerprint
def smiles_to_morgan(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        arr = np.zeros((1,))
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return None

# Generate fingerprints
fingerprints = []
for smi in df['SMILES']:
    if pd.notna(smi):
        fingerprints.append(smiles_to_morgan(smi))
    else:
        fingerprints.append(None)

# Add fingerprints back to dataframe
df['Fingerprint'] = fingerprints

# Drop rows where fingerprint couldn't be generated
df = df.dropna(subset=['Fingerprint'])

# Save fingerprints as a separate NumPy array
X = np.array(df['Fingerprint'].tolist())

# Save X as a numpy file for training
np.save('X_features.npy', X)

# Save drug names (optional, useful for matching back later)
df[['Drug', 'SMILES']].to_csv('drug_cleaned.csv', index=False)

print("✅ Feature preparation complete. Morgan fingerprints saved as 'X_features.npy'.")
