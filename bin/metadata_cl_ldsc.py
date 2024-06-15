import os
import re
import sys

prefix = sys.argv[1]
data_dir = sys.argv[2]
output_file = sys.argv[3]

# Collect filenames
filenames = os.listdir(data_dir)

# Extract unique parts after the prefix
unique_labels = set()
pattern = re.compile(r'^{}\.([0-9_]+)\..*'.format(re.escape(prefix)))

for filename in filenames:
    match = pattern.match(filename)
    if match:
        unique_labels.add(match.group(1))
    else:
        raise ValueError("Filename {} does not match the pattern.".format(filename))

# Write to output file
with open(output_file, 'w') as f:
    for label in sorted(unique_labels):
        dat = "{}.{}".format(prefix, label)
        f.write("{}\t{}/{}.\n".format(dat, data_dir, dat))
