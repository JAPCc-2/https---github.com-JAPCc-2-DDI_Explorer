import requests
import pandas as pd

# Example list of drugs
drug_list = [
    "Atorvastatin", "Simvastatin", "Rosuvastatin", 
    "Pravastatin", "Fluvastatin", "Ezetimibe", 
    "Lovastatin", "Pitavastatin", "Fenofibrate", 
    "Gemfibrozil", "Niacin"
]

# Function to retrieve SMILES from PubChem
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

# Fetch SMILES for all drugs
data = []
for drug in drug_list:
    smiles = get_smiles(drug)
    data.append({
        "Drug": drug,
        "SMILES": smiles
    })

# Create DataFrame
df = pd.DataFrame(data)

# Show results
print(df)

# Save to CSV
df.to_csv("drug_smiles.csv", index=False)
print("✅ Data saved to drug_smiles.csv")
