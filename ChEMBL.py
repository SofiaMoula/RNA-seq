import pandas as pd
from chembl_webresource_client.new_client import new_client
import time
import aiohttp
import asyncio
import xml.etree.ElementTree as ET

print("Script started.")

# Load your DEGs
try:
    degs = pd.read_csv('path/to/file/top100genes.csv')['id'].tolist()[:100]  # Limiting to 100 genes for testing
    print("DEGs loaded.")
except Exception as e:
    print(f"Error loading DEGs: {e}")
    degs = []

# Initialize ChEMBL client
target = new_client.target
activity = new_client.activity

# Function to retrieve ChEMBL targets for a list of genes
def get_chembl_targets(genes):
    targets = []
    for i, gene in enumerate(genes):
        try:
            results = target.filter(pref_name__icontains=gene)
            for result in results:
                if result['organism'] == 'Homo sapiens':  # Filter by organism
                    targets.append({
                        'gene_symbol': gene,
                        'target_chembl_id': result['target_chembl_id'],
                        'pref_name': result['pref_name'],
                        'organism': result['organism']
                    })
        except Exception as e:
            print(f"Error retrieving targets for gene {gene}: {e}")
        if i % 10 == 0:
            print(f"Processed {i} genes out of {len(genes)}")
    return pd.DataFrame(targets)

# Asynchronous function to retrieve compounds for a target using aiohttp
async def fetch_compounds(session, target_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id={target_id}"
    compounds = []
    try:
        async with session.get(url, timeout=60) as response:  # Increased timeout
            if response.status == 200:
                content_type = response.headers.get('Content-Type')
                if 'application/json' in content_type:
                    activities = await response.json()
                    if 'activities' in activities:
                        for act in activities['activities']:
                            compounds.append({
                                'target_chembl_id': target_id,
                                'molecule_chembl_id': act.get('molecule_chembl_id', ''),
                                'activity_type': act.get('standard_type', ''),
                                'activity_value': act.get('standard_value', ''),
                                'activity_units': act.get('standard_units', ''),
                                'document_journal': act.get('journal', ''),
                                'document_year': act.get('year', '')
                            })
                    else:
                        print(f"No activities found for target {target_id}")
                elif 'application/xml' in content_type:
                    text = await response.text()
                    root = ET.fromstring(text)
                    activities = root.findall(".//activity")
                    for act in activities:
                        compounds.append({
                            'target_chembl_id': target_id,
                            'molecule_chembl_id': act.findtext('molecule_chembl_id', ''),
                            'activity_type': act.findtext('standard_type', ''),
                            'activity_value': act.findtext('standard_value', ''),
                            'activity_units': act.findtext('standard_units', ''),
                            'document_journal': act.findtext('journal', ''),
                            'document_year': act.findtext('year', '')
                        })
                else:
                    print(f"Unexpected content type for target {target_id}: {content_type}")
            else:
                print(f"Failed to fetch data for target {target_id} with status {response.status}")
    except asyncio.TimeoutError:
        print(f"Timeout error for target {target_id}")
    except Exception as e:
        print(f"Error retrieving compounds for target {target_id}: {e}")
    await asyncio.sleep(2)  # Increased delay to avoid rate limiting
    return compounds

# Main function to orchestrate the fetching process
async def main():
    async with aiohttp.ClientSession() as session:
        # Fetch targets for all genes
        chembl_targets = get_chembl_targets(degs)

        # Fetch compounds for all targets concurrently
        tasks = [fetch_compounds(session, target_id) for target_id in chembl_targets['target_chembl_id']]
        all_compounds = []
        for i, task in enumerate(asyncio.as_completed(tasks)):
            compounds = await task
            all_compounds.extend(compounds)
            if i % 10 == 0:
                print(f"Processed {i} targets out of {len(chembl_targets)}")

        compounds_df = pd.DataFrame(all_compounds)

        # Save results to CSV files
        chembl_targets.to_csv('chembl_targets.csv', index=False)
        compounds_df.to_csv('chembl_compounds.csv', index=False)

        print("Targets and compounds retrieved and saved to CSV files.")

# Run the main function
asyncio.run(main())
