#!/usr/bin/env python

"""
One-time script to convert the SLC locations CSV from compound
location strings to a simple individual-location-to-code mapping.

The original file has rows like:
    "ER; Plasma membrane",e
    "ER; Plasma membrane",r

This script extracts the implicit mapping:
    ER,r
    Plasma membrane,e

by splitting compound locations and pairing each part with its code.
"""

import csv
import re
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / (
    'omnipath_metabo/datasets/cosmos/data'
)

src = DATA_DIR / 'slc_locations.csv'

mapping: dict[str, str] = {}

with open(src) as f:
    reader = csv.DictReader(f)

    for row in reader:
        parts = re.split(r';\s*', row['location'])
        code = row['abbreviation']

        # For single-part entries, the mapping is unambiguous
        if len(parts) == 1:
            loc = parts[0].strip()
            if loc and loc != 'Unknown':
                mapping[loc] = code

# Print the extracted mapping for review
print('Individual location -> code mapping extracted from compound entries:')
print()

for loc in sorted(mapping):
    print(f'  {loc} -> {mapping[loc]}')

# Write the new CSV
dst = DATA_DIR / 'slc_locations.csv'

with open(dst, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['location', 'abbreviation'])

    for loc in sorted(mapping):
        writer.writerow([loc, mapping[loc]])

print(f'\nWrote {len(mapping)} entries to {dst}')
