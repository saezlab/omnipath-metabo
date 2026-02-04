# we extract species from pypath-mrclinksdb

import pandas as pd
import re
import os
import sys

# Add pypath to path
sys.path.insert(0, '/Users/priscillabai/pypath-tool')
from pypath.inputs.mrclinksdb import _interactions

def parse_uniprot_locations(location_string):
    """Parse UniProt location string to extract location and features"""
    # Handle NaN or non-string values
    if pd.isna(location_string) or not isinstance(location_string, str):
        return []

    # Remove outer curly braces
    cleaned = re.sub(r'^\{|\}$', '', location_string)

    # Pattern to match UniprotLocation entries
    pattern = r"UniprotLocation\(location='([^']+)', features=([^)]+)\)"
    matches = re.findall(pattern, cleaned)

    results = []
    for location, features_raw in matches:
        # Clean features
        features_clean = re.sub(r"\('|'\)|'|,\s*$|^\(|\)$", "", features_raw)
        results.append({
            'location': location,
            'features': features_clean
        })

    return results


def process_mrclinksdb(species='human', output_dir="/Users/priscillabai/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project/COSMOS_PKN"):
    """
    Process MRClinksDB interactions for a specific species.

    Args:
        species: Species name (e.g., 'human', 'mouse') or NCBI taxonomy ID
        output_dir: Base directory for output files
    """
    # Set working directory
    os.chdir(output_dir)

    # Load data from pypath mrclinksdb module
    print(f"Loading MRClinksDB data for {species}...")
    interactions = list(_interactions.mrclinksdb_interaction(organism=species))

    # Convert to DataFrame with required columns
    data = pd.DataFrame([
        {
            'Target': rec.pubchem,
            'Source': str(rec.receptor_uniprot),
            'location': rec.receptor_location
        }
        for rec in interactions
    ])

    # Load abbreviation data
    abb_data = pd.read_csv("data/TCDB/protein_location_abb.csv", sep=";")
    abb_data = abb_data.iloc[:, :2]
    abb_data.columns = ['location', 'abbreviation']
    abb_data['location'] = abb_data['location'].str.strip()

    # Parse locations for each row
    all_parsed_locations = []
    for idx, row in data.iterrows():
        parsed = parse_uniprot_locations(row['location'])
        for item in parsed:
            item['original_row'] = str(idx + 1)  # R uses 1-based indexing
            all_parsed_locations.append(item)

    final_result = pd.DataFrame(all_parsed_locations)

    # Left join with abbreviation data and filter
    final_result = final_result.merge(abb_data, on='location', how='left')
    final_result = final_result[final_result['abbreviation'].notna()]

    # Group by original_row to get unique locations
    location_summary = final_result.groupby('original_row')['abbreviation'].apply(lambda x: list(x.unique())).reset_index()
    location_summary.columns = ['original_row', 'locations']

    # Add row_id to data
    data['row_id'] = (data.index + 1).astype(str)

    # Left join with location summary
    result_table = data.merge(location_summary, left_on='row_id', right_on='original_row', how='left')

    # Filter out rows with no locations and explode locations
    result_table = result_table[result_table['locations'].notna()]
    result_table = result_table.explode('locations')

    # Reset index and create Gene/Metab identifiers
    result_table = result_table.reset_index(drop=True)
    result_table['Source_new'] = 'Gene' + (result_table.index + 1).astype(str) + '__' + result_table['Source']
    result_table['Target_new'] = 'Metab__' + result_table['Target'] + '_' + result_table['locations']

    # Select and rename
    result_table = result_table[['Target_new', 'Source_new']].copy()
    result_table.columns = ['Source', 'Target']

    # add forward reaction and reverse reaction pairs
    table_list = []
    for idx, row in result_table.iterrows():
        source = row['Source']
        target = row['Target']
        # Replace the last _[a-z]+ with _c in the source
        target_modified = re.sub(r'_[a-z]+$', '_c', source)

        # Add the two rows
        table_list.append({'Source': source, 'Target': target})
        table_list.append({'Source': target, 'Target': target_modified})

    table = pd.DataFrame(table_list)

    # Create table_extended with batch processing
    table_extended_list = []
    for i in range(0, len(table), 2):
        if i + 1 < len(table):
            source1, target1 = table.iloc[i]['Source'], table.iloc[i]['Target']
            source2, target2 = table.iloc[i+1]['Source'], table.iloc[i+1]['Target']

            # Add four rows per batch
            table_extended_list.append({'Source': source1, 'Target': target1})
            table_extended_list.append({'Source': source2, 'Target': target2})
            table_extended_list.append({'Source': target2, 'Target': target1 + '_reverse'})
            table_extended_list.append({'Source': target1 + '_reverse', 'Target': source1})

    table_extended = pd.DataFrame(table_extended_list)

    ## add gene nodes
    gene_ids = pd.concat([table_extended['Source'], table_extended['Target']]).unique()
    gene_ids = [gid for gid in gene_ids if str(gid).startswith('Gene')]

    gene_mapping_list = []
    for gene_id in gene_ids:
        source = re.sub(r'.*__', '', gene_id)
        source = re.sub(r'_reverse$', '', source)
        gene_mapping_list.append({'Source': source, 'Target': gene_id})

    gene_mapping = pd.DataFrame(gene_mapping_list)

    ## add metabolite nodes
    metab_ids = pd.concat([table_extended['Source'], table_extended['Target']]).unique()
    metab_ids = [mid for mid in metab_ids if str(mid).startswith('Metab')]

    metabolite_mapping_list = []
    for metab_id in metab_ids:
        source = re.sub(r'^Metab__|_[a-z]+$', '', metab_id)
        metabolite_mapping_list.append({'Source': source, 'Target': metab_id})

    metabolite_mapping = pd.DataFrame(metabolite_mapping_list)

    # Combine all tables
    table_extended = pd.concat([table_extended, gene_mapping, metabolite_mapping], ignore_index=True)
    table_extended['mor'] = 1

    # Write to CSV
    output_file = f"result/mrclinksdb/mrclinksdb_COMOS_python.csv"
    table_extended.to_csv(output_file, index=False)

    print(f"Output written to {output_file}")
    print(f"Total rows: {len(table_extended)}")

    return table_extended


# Main execution
if __name__ == "__main__":
    # Run with default species (human)
    result = process_mrclinksdb(species='human')
