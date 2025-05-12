# Import the function from prepare_features.py
from prepare_features import smiles_to_morgan

# Sample SMILES string (representing a drug molecule)
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

# Call the function and print the Morgan fingerprint
fingerprint = smiles_to_morgan(smiles)
print(fingerprint)