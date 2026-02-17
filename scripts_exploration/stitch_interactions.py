"""
STITCH Chemical-Protein Interactions for COSMOS PKN

This script processes STITCH chemical-protein interaction data for integration
into the COSMOS prior-knowledge network. It uses the new_stitch client from
pypath for improved data parsing and Cloudflare compatibility.

Key features:
- Uses merged interactions (actions + links) for complete metadata
- Filters by confidence score thresholds (final_score)
- Filters by interaction mode (binding, expression, etc.)
- Adds Gene__ and Metab__ prefixes following COSMOS conventions
- Normalizes orientation: chemical as source, protein as target
- Captures activation/inhibition effects (mor: 1/-1/0)
"""

import sys
sys.path.insert(0, '/Users/priscillabai/pypath-tool')

import pandas as pd
from pathlib import Path
from pypath.inputs.new_stitch import interactions as stitch_interactions


def process_stitch_interactions(
        ncbi_tax_id=9606,
        score_threshold=700,
        mode=('activation', 'inhibition')
    ):
    """
    Fetch and process STITCH chemical-protein interaction data.

    Args:
        ncbi_tax_id: NCBI taxonomy ID (default: 9606 for human)
        score_threshold: Minimum final_score for interactions (default: 700)
        mode: Interaction mode(s) to keep (default: ('activation', 'inhibition'))
              Pass a string or tuple/list of strings. Options per record:
              'binding', 'pred_bind', 'expression', 'activation',
              'inhibition', 'reaction', 'catalysis'. Pass None to keep all modes.

    Returns:
        DataFrame with processed STITCH interactions in COSMOS format
    """

    print(f"Fetching STITCH interactions from new_stitch client...")
    print(f"  Species (NCBI Tax ID): {ncbi_tax_id}")
    print(f"  Score threshold (final_score): {score_threshold}")
    print(f"  Mode filter: {mode if mode else 'all modes'}")

    # Collect interaction data
    records = []

    for rec in stitch_interactions(ncbi_tax_id=ncbi_tax_id):

        # Filter by score threshold using final_score
        if rec.final_score < score_threshold:
            continue

        # Filter by mode if specified (supports single string or list/tuple)
        if mode:
            allowed = (mode,) if isinstance(mode, str) else mode
            if rec.mode not in allowed:
                continue

        # Normalize orientation: chemical as source, protein as target
        if rec.source.type == 'small_molecule':
            chemical_id = rec.source.id
            protein_id = rec.target.id
        elif rec.target.type == 'small_molecule':
            chemical_id = rec.target.id
            protein_id = rec.source.id
        else:
            # Skip if no small molecule (shouldn't happen in STITCH)
            continue

        # Determine mor (mode of regulation) based on effect
        if rec.activation:
            mor = 1
        elif rec.inhibition:
            mor = -1
        else:
            mor = 0

        # Create COSMOS-formatted node IDs
        source = f'Metab__{chemical_id}'
        target = f'Gene__{protein_id}'

        records.append({
            'Source': source,
            'Target': target,
            'mor': mor,
            'database': 'STITCH'
        })

    data = pd.DataFrame(records)
    print(f"Total records collected after filtering: {len(data)}")

    if len(data) == 0:
        print("Warning: No interactions found after filtering")
        return pd.DataFrame()

    print(f"STITCH processed: {len(data)} rows")

    return data


def main(
        ncbi_tax_id=9606,
        score_threshold=700,
        mode=('activation', 'inhibition')
    ):
    """
    Main execution function for STITCH COSMOS PKN integration.

    Args:
        ncbi_tax_id: NCBI taxonomy ID (default: 9606 for human)
        score_threshold: Minimum final_score (default: 700 for high confidence)
        mode: Interaction mode(s) to keep (default: ('activation', 'inhibition'))
              Pass a string, tuple/list of strings, or None for all modes.
    """

    base_path = "/Users/priscillabai/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project"
    output_dir = Path(f"{base_path}/COSMOS_PKN/result/stitch")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process STITCH interactions
    print("\n" + "="*60)
    print("Processing STITCH interactions for COSMOS PKN")
    print("="*60)

    stitch_data = process_stitch_interactions(
        ncbi_tax_id=ncbi_tax_id,
        score_threshold=score_threshold,
        mode=mode
    )

    if len(stitch_data) == 0:
        print("No data to save")
        return

    # Save results
    print("\nSaving results...")
    if mode is None:
        mode_label = 'all'
    elif isinstance(mode, str):
        mode_label = mode
    else:
        mode_label = '_'.join(mode)
    output_file = output_dir / f"stitch_processed_score{score_threshold}_mode_{mode_label}.csv"
    stitch_data.to_csv(output_file, index=False)

    print(f"\nProcessing complete!")
    print(f"Output file: {output_file}")
    print(f"Total interactions: {len(stitch_data)}")

    # Show summary statistics
    print("\n" + "="*60)
    print("Summary Statistics:")
    print("="*60)
    print(f"Mode of Regulation (mor) distribution:")
    print(stitch_data['mor'].value_counts().sort_index())
    print(f"\nSample rows:")
    print(stitch_data.head(10).to_string())


if __name__ == "__main__":
    # Default: score threshold 700, keep activation and inhibition, human (9606)
    main(ncbi_tax_id=9606, score_threshold=700, mode=('activation', 'inhibition'))

    # Alternative usage examples:

    # Medium confidence threshold
    # main(ncbi_tax_id=9606, score_threshold=400, mode=('activation', 'inhibition'))

    # Include all interaction modes
    # main(ncbi_tax_id=9606, score_threshold=700, mode=None)

    # Single mode
    # main(ncbi_tax_id=9606, score_threshold=700, mode='activation')
