#!/usr/bin/env python3
import argparse
import pandas as pd
from antigens_type import ANTIGENS_TYPE

output_file = "classified_specificity_labeled.csv"

# ========================================
# Parse arguments
# ========================================

parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True)
args = parser.parse_args()

# ========================================
# Set type file
# ========================================

df = pd.read_csv(args.input)

for type_antigen, list_antigens in ANTIGENS_TYPE.items():
    mask = df["gene_name"].isin(list_antigens)
    df.loc[mask, "type"] = type_antigen

df.to_csv(output_file, index=None)
