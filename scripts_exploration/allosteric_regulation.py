import sys
sys.path.insert(0, '/Users/priscillabai/pypath-tool')

import pandas as pd
import os
from pypath.inputs.brenda._main import allosteric_regulation


def generate_brenda_allosteric(
        species='human',
        output_dir="/Users/priscillabai/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project/COSMOS_PKN/data/brenda"
    ):
    """
    Generate BRENDA allosteric regulation data for a specific species.

    Args:
        species: Species name (e.g., 'human', 'mouse') or NCBI taxonomy ID
        output_dir: Directory for output files

    Returns:
        DataFrame with allosteric regulation data
    """
    os.makedirs(output_dir, exist_ok=True)

    # Get allosteric regulation data
    records = []
    for record in allosteric_regulation(organisms=[species], limit=None):
        if record.protein:
            for protein_id in record.protein:
                records.append({
                    'compound': f"Metab__{record.compound}",
                    'protein_id': f"Gene__{protein_id}",
                    'action': record.action,
                    'id_type': record.id_type,
                    'pubmeds': record.pubmeds
                })

    # Create DataFrame with specified column order
    df = pd.DataFrame(records, columns=['compound', 'protein_id', 'action', 'id_type', 'pubmeds'])

    # Save to CSV
    output_file = f"{output_dir}/brenda_allosteric_{species}.csv"
    df.to_csv(output_file, index=False)

    return df


if __name__ == "__main__":
    result = generate_brenda_allosteric(species='human')
