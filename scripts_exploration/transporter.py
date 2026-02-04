import sys
sys.path.insert(0, '/Users/priscillabai/pypath-tool')

import pandas as pd
import re
from pathlib import Path
from pypath.inputs.tcdb import _substrates as tcdb_substrates
from pypath.inputs import slc
from pypath.utils import reflists


def parse_uniprot_locations(location_string):
    """Parse UniProt location strings and extract location and features."""
    pattern = r"UniprotLocation\(location='([^']+)', features=([^)]+)\)"
    matches = re.findall(pattern, re.sub(r'^\{|\}$', '', str(location_string)))
    if not matches:
        return pd.DataFrame(columns=['location', 'features'])
    return pd.DataFrame([
        {'location': m[0], 'features': re.sub(r"\('|'\)|'|,\s*$|^\(|\)$", "", m[1])}
        for m in matches
    ])


def load_location_mapping(filepath, sep=";", columns=None):
    """Load and clean location abbreviation mapping."""
    df = pd.read_csv(filepath, sep=sep, header=None if columns else 0)
    df = df.iloc[:, :2]
    df.columns = columns or df.columns[:2]
    if 'location' in df.columns:
        df['location'] = df['location'].str.strip()
    return df


# TCDB is a multi-species database and all entries from different species mixed together. We extract each species based on protein uniprot ID
def process_tcdb_data(abb_data, ncbi_tax_id=9606):
    """Fetch and process TCDB substrate data for a specific species.

    Args:
        abb_data: Location abbreviation mapping dataframe
        ncbi_tax_id: NCBI taxonomy ID (default: 9606 for human)
    """
    print(f"Fetching TCDB substrate data from pypath (species: {ncbi_tax_id})...")
    species_proteins = set(reflists.get_reflist('uniprot', ncbi_tax_id=ncbi_tax_id))

    # Collect species-specific protein substrate data
    data = pd.DataFrame([
        {'Source': r.transporter_uniprot, 'Target': r.substrate_id,
         'substrate_name': r.substrate_name, 'location': r.location}
        for r in tcdb_substrates.tcdb_substrate()
        if r.transporter_uniprot in species_proteins
    ])
    print(f"Total records collected: {len(data)}")

    # Parse and map locations to abbreviations
    location_results = [
        parse_uniprot_locations(loc).assign(original_row=str(idx))
        for idx, loc in enumerate(data['location'])
    ]
    location_results = [df for df in location_results if not df.empty]

    if location_results:
        final_result = (pd.concat(location_results, ignore_index=True)
                       .merge(abb_data, on='location', how='left')
                       .dropna(subset=['abbreviation']))
        location_summary = (final_result.groupby('original_row')['abbreviation']
                          .apply(lambda x: list(x.unique())).reset_index(name='locations'))
    else:
        location_summary = pd.DataFrame(columns=['original_row', 'locations'])

    # Create final TCDB dataframe with database column
    result = (data.assign(row_id=data.index.astype(str))
             .merge(location_summary, left_on='row_id', right_on='original_row', how='left')
             .dropna(subset=['locations'])
             .explode('locations')
             .assign(Source=lambda df: 'Gene__' + df['Source'],
                    Target=lambda df: 'Metab__' + df['Target'] + '_' + df['locations'],
                    database='TCDB')
             [['Target', 'Source', 'database']]
             .rename(columns={'Target': 'Source', 'Source': 'Target'}))

    print(f"TCDB processed: {len(result)} rows")
    return result

# slc is only for human
def process_slc_data(abb_data):
    """Fetch and process SLC interaction data."""
    print("Fetching SLC interaction data from pypath...")

    # Collect SLC interaction data
    data = pd.DataFrame([
        {'transporter': r.transporter.uniprot, 'substrate': r.substrate.chebi,
         'localization': r.localization}
        for r in slc.slc_interactions()
        if r.substrate.chebi and r.localization and r.localization not in ('', 'Unknown', 'None')
    ])
    print(f"Total SLC records collected: {len(data)}")

    # Create final SLCtable dataframe with database column
    result = (data.merge(abb_data, on='localization', how='left')
             .dropna(subset=['abbreviation'])
             .assign(Source=lambda df: 'Gene__' + df['transporter'],
                    Target=lambda df: 'Metab__' + df['substrate'] + '_' + df['abbreviation'],
                    database='SLC')
             [['Target', 'Source', 'database']]
             .rename(columns={'Target': 'Source', 'Source': 'Target'}))

    print(f"SLCTable processed: {len(result)} rows")
    return result


def add_reverse_reactions(table):
    """Add forward and reverse reaction pairs with _c suffix for metabolites.

    Based on TCDB_COSMOS.py lines 82-132.

    Args:
        table: DataFrame with Source, Target, and database columns

    Returns:
        Dictionary containing:
        - table_extended: Extended table with reverse reactions
        - gene_mapping: Gene node mappings
        - metabolite_mapping: Metabolite node mappings
    """
    print("Adding reverse reaction pairs...")

    # Step 1: Add forward and reverse reaction pairs
    table_list = []
    for idx, row in table.iterrows():
        source = row['Source']
        target = row['Target']
        database = row['database']

        # Replace the last _[a-z]+ with _c in the source (metabolite)
        target_modified = re.sub(r'_[a-z]+$', '_c', source)

        # Add the two rows
        table_list.append({'Source': source, 'Target': target, 'database': database})
        table_list.append({'Source': target, 'Target': target_modified, 'database': database})

    table_pairs = pd.DataFrame(table_list)
    print(f"After adding reverse pairs: {len(table_pairs)} rows")

    # Step 2: Create extended table with reverse reactions
    table_extended_list = []
    for i in range(0, len(table_pairs), 2):
        if i + 1 < len(table_pairs):
            source1, target1, db1 = table_pairs.iloc[i]['Source'], table_pairs.iloc[i]['Target'], table_pairs.iloc[i]['database']
            source2, target2, db2 = table_pairs.iloc[i+1]['Source'], table_pairs.iloc[i+1]['Target'], table_pairs.iloc[i+1]['database']

            # Add four rows per batch
            table_extended_list.append({'Source': source1, 'Target': target1, 'database': db1})
            table_extended_list.append({'Source': source2, 'Target': target2, 'database': db2})
            table_extended_list.append({'Source': target2, 'Target': target1 + '_reverse', 'database': db2})
            table_extended_list.append({'Source': target1 + '_reverse', 'Target': source1, 'database': db1})

    table_extended = pd.DataFrame(table_extended_list)
    print(f"After adding extended reverse reactions: {len(table_extended)} rows")

    # Step 3: Add gene nodes
    gene_ids = pd.concat([table_extended['Source'], table_extended['Target']]).unique()
    gene_ids = [gid for gid in gene_ids if str(gid).startswith('Gene')]

    gene_mapping_list = []
    for gene_id in gene_ids:
        source = re.sub(r'.*__', '', gene_id)
        source = re.sub(r'_reverse$', '', source)
        gene_mapping_list.append({'Source': source, 'Target': gene_id})

    gene_mapping = pd.DataFrame(gene_mapping_list)
    print(f"Gene node mappings: {len(gene_mapping)} nodes")

    # Step 4: Add metabolite nodes
    metab_ids = pd.concat([table_extended['Source'], table_extended['Target']]).unique()
    metab_ids = [mid for mid in metab_ids if str(mid).startswith('Metab')]

    metabolite_mapping_list = []
    for metab_id in metab_ids:
        source = re.sub(r'^Metab__|_[a-z]+$', '', metab_id)
        metabolite_mapping_list.append({'Source': source, 'Target': metab_id})

    metabolite_mapping = pd.DataFrame(metabolite_mapping_list)
    print(f"Metabolite node mappings: {len(metabolite_mapping)} nodes")

    return {
        'table_extended': table_extended,
        'gene_mapping': gene_mapping,
        'metabolite_mapping': metabolite_mapping
    }


def calculate_intersections(slc_table, tcdb_table):
    """Calculate intersections between SLC and TCDB datasets."""
    print("\nCalculating intersections...")

    # With location
    combo1 = slc_table['Source'] + '_' + slc_table['Target']
    combo2 = tcdb_table['Source'] + '_' + tcdb_table['Target']
    with_loc = len(set(combo1) & set(combo2))
    print(f"Intersection with location: {with_loc}")

    # Without location
    combo1_no_loc = slc_table['Source'].str.replace(r'_[^_]+$', '', regex=True) + '_' + slc_table['Target']
    combo2_no_loc = tcdb_table['Source'].str.replace(r'_[^_]+$', '', regex=True) + '_' + tcdb_table['Target']
    without_loc = len(set(combo1_no_loc) & set(combo2_no_loc))
    print(f"Intersection without location: {without_loc}")

    return with_loc, without_loc


def main(ncbi_tax_id=9606):
    """Main execution function.

    Args:
        ncbi_tax_id: NCBI taxonomy ID (default: 9606 for human)
    """
    base_path = "/Users/priscillabai/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project"

    print("Processing TCDB data...")
    tcdb_abb = load_location_mapping(
        f"{base_path}/MetaboliteDB/TCDB/protein_location_abb.csv",
        sep=";", columns=['location', 'abbreviation']
    )
    TCDB = process_tcdb_data(tcdb_abb, ncbi_tax_id=ncbi_tax_id)

    # SLC is only available for human (ncbi_tax_id=9606)
    if ncbi_tax_id == 9606:
        print("\nProcessing SLCTable data...")
        slc_abb = load_location_mapping(
            f"{base_path}/COSMOS_PKN/result/SLCtable/protein_location_abb.txt",
            sep="\t", columns=['localization', 'abbreviation']
        )
        SLCtable = process_slc_data(slc_abb)

        calculate_intersections(SLCtable, TCDB)

        # Combine TCDB and SLC, remove duplicates
        print("\nCombining TCDB and SLC data...")
        combined = pd.concat([TCDB, SLCtable], ignore_index=True).drop_duplicates()
        combined['mor'] = 1
        print(f"Combined dataset: {len(combined)} unique rows (TCDB: {len(TCDB)}, SLC: {len(SLCtable)})")
    else:
        print(f"\nSkipping SLC data (only available for human, current species: {ncbi_tax_id})")
        combined = TCDB
        combined['mor'] = 1
        SLCtable = None

    # Add reverse reactions and create mappings
    print("\n" + "="*60)
    result_dict = add_reverse_reactions(combined)
    table_extended = result_dict['table_extended']
    gene_mapping = result_dict['gene_mapping']
    metabolite_mapping = result_dict['metabolite_mapping']

    print("\nSaving results...")
    output_dir = Path(f"{base_path}/COSMOS_PKN/result/transporter")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save individual datasets
    TCDB.to_csv(output_dir / "TCDB_processed.csv", index=False)
    if SLCtable is not None:
        SLCtable.to_csv(output_dir / "SLCtable_processed.csv", index=False)

    # Save combined and extended results
    combined.to_csv(output_dir / "transporter_combined.csv", index=False)


if __name__ == "__main__":
    main()
