import os
import re
import sys

prefix = sys.argv[1]
data_dir = sys.argv[2]
output_file = sys.argv[3]

# prefix = "p_weights.28.mixing"
# data_dir = "data_files"
# output_file = "per_sample.ldcts"

filenames = os.listdir(data_dir)

# Extract unique parts after the prefix
unique_labels = set()
pattern = re.compile(f"^^{re.escape(prefix)}\\.([0-9+_]+)\\.*")

for filename in filenames:
    match = pattern.match(filename)
    if match:
        unique_labels.add(match.group(1))
    else:
        raise ValueError(f"Filename {filename} does not match the pattern.")

# Write to output file
with open(output_file, 'w') as f:
    for label in sorted(unique_labels):
        dat = f"{prefix}.{label}"
        f.write(f"{dat}\t{data_dir}/{dat}.\n")
