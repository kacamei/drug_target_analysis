# drug_target_analysis.py

import pandas as pd
import requests
import concurrent.futures
import csv
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

# Step 1: Fetch and organize approved drugs
print("Retrieving approved drug records from ChEMBL...")
molecule_client = new_client.molecule
approved_molecules = molecule_client.filter(max_phase=4, first_approval__isnull=False)

drug_entries = []
for entry in tqdm(approved_molecules, desc="Processing approved molecules"):
    if entry['molecule_chembl_id'] and entry['first_approval']:
        drug_entries.append((entry['first_approval'], entry.get('pref_name', ''), entry['molecule_chembl_id']))

# Sort by approval year and alphabetical name
sorted_drugs = sorted(drug_entries, key=lambda x: (x[0], x[1]))

# Save approved drug list
with open('approved_drugs_list.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(["Approval_Year", "Drug_Name", "ChEMBL_ID"])
    writer.writerows(sorted_drugs)

print(f"✔️ Stored {len(sorted_drugs)} approved drugs into approved_drugs_list.csv.\n")

# Step 2: Retrieve targets for drugs approved from 2019 onwards
print("Identifying recent drug targets (approved from 2019 onwards)...")
approved_drugs_df = pd.read_csv('approved_drugs_list.csv')
recent_approvals = approved_drugs_df[approved_drugs_df['Approval_Year'] >= 2019]

target_lookup = new_client.target
activity_lookup = new_client.activity

drug_to_protein = {}

for _, drug in tqdm(recent_approvals.iterrows(), total=recent_approvals.shape[0], desc="Fetching targets"):
    chembl_id = drug['ChEMBL_ID']
    activities = activity_lookup.filter(molecule_chembl_id=chembl_id).only(['target_chembl_id'])
    target_ids = list({a['target_chembl_id'] for a in activities if a['target_chembl_id']})

    for target_id in target_ids:
        target_details = target_lookup.get(target_id)
        if 'target_components' in target_details:
            for component in target_details['target_components']:
                accession = component.get('accession')
                if accession:
                    drug_to_protein.setdefault(chembl_id, []).append(accession)

# Save drug-to-protein mappings
protein_mappings = []
for drug_id, accessions in drug_to_protein.items():
    for accession in accessions:
        protein_mappings.append((drug_id, accession))

protein_df = pd.DataFrame(protein_mappings, columns=["ChEMBL_ID", "UniProt_Accession"])
protein_df.to_csv('drug_targets.csv', index=False)

print(f"✔️ Collected {len(protein_mappings)} protein targets into drug_targets.csv.\n")


# Step 3: Fetch UniProt keywords concurrently
def retrieve_keywords(accession):
    endpoint = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    headers = {"User-Agent": "Mozilla/5.0"}

    try:
        response = requests.get(endpoint, headers=headers, timeout=10)
        if response.status_code == 200:
            record = response.json()

            keywords = []

            if 'keywords' in record and isinstance(record['keywords'], list):
                keywords.extend([kw.get('value') for kw in record['keywords'] if 'value' in kw])

            if not keywords and 'comments' in record:
                for comment in record['comments']:
                    if comment.get('commentType') == 'FUNCTION' and 'texts' in comment:
                        for text in comment['texts']:
                            if 'value' in text:
                                keywords.append(text['value'])

            return (accession, keywords)
        else:
            print(f"⚠️ HTTP {response.status_code} for {accession}")
            return (accession, [])
    except Exception as e:
        print(f"⚠️ Error retrieving {accession}: {e}")
        return (accession, [])


print("Extracting UniProt keywords for protein targets...")

accessions = protein_df['UniProt_Accession'].unique()
keyword_annotations = {}

with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
    futures = {executor.submit(retrieve_keywords, acc): acc for acc in accessions}
    for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Fetching keywords"):
        accession, keywords = future.result()
        keyword_annotations[accession] = keywords

# Save keywords
keyword_data = []
for accession, keywords in keyword_annotations.items():
    for keyword in keywords:
        keyword_data.append((accession, keyword))

keywords_df = pd.DataFrame(keyword_data, columns=["UniProt_Accession", "Keyword"])
keywords_df.to_csv('protein_keywords.csv', index=False)

print(f"✔️ Retrieved {len(keyword_data)} keywords from {len(accessions)} proteins.")
print("✔️ Saved keyword annotations into protein_keywords.csv.\n")

# FINAL PRINT (No emoji to avoid Unicode crash)
print("🎯 Analysis complete! You can now continue with your project.")
