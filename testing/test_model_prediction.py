import pickle
import numpy as np
from fetch_smiles import get_smiles
from prepare_features import smiles_to_morgan

# Load the trained model
model = pickle.load(open('random_forest_ddi_model.pkl', 'rb'))

# Example drugs
drug1 = "Atorvastatin"
drug2 = "Simvastatin"

# Get SMILES strings for the drugs
smiles1 = get_smiles(drug1)
smiles2 = get_smiles(drug2)

# Generate the Morgan fingerprints
fp1 = smiles_to_morgan(smiles1)
fp2 = smiles_to_morgan(smiles2)

# Calculate the absolute difference between the fingerprints
input_features = np.abs(fp1 - fp2).reshape(1, -1)

# Make predictions
prediction = model.predict(input_features)[0]
probability = model.predict_proba(input_features)[0][prediction]

# Print the result
print(f"Interaction: {'Interaction Predicted' if prediction == 1 else 'No Major Interaction'}")
print(f"Confidence: {round(float(probability), 3)}")