# Drug Target Analysis

This project retrieves and analyzes drug-target information using the ChEMBL and UniProt databases.  
It automates the extraction of approved drugs, maps them to their target proteins, and fetches functional keywords.

---

## Workflow

1. **Fetch Approved Drugs**  
   Retrieve drugs approved for clinical use from ChEMBL.

2. **Identify Drug Targets**  
   Find protein targets associated with the approved drugs.

3. **Annotate Proteins**  
   Extract UniProt keywords for each target protein.

---

## Outputs

- `approved_drugs_list.csv` — List of approved drugs from ChEMBL.
- `drug_targets.csv` — Mapping of drugs to their protein targets.
- `protein_keywords.csv` — Keywords for each target protein.

---

## Technologies Used

- Python
- ChEMBL Webresource Client
- UniProt REST API
- pandas
- tqdm
- requests
- concurrent.futures

---

## Setup

```bash
pip install chembl_webresource_client pandas tqdm requests
python drug_target_analysis.py
