import sys
import pandas as pd

# Retrieve the file paths from command line arguments
file_paths = sys.argv[1:-1]  # Skip the script name, and take all arguments except the last one
output_file = sys.argv[-1]   # The last argument is the output file path

# Initialize an empty DataFrame for counts
counts = pd.DataFrame()

# Keep track of files with and without 'NumReads'
files_with_numreads = []
files_without_numreads = []

# Iterate through each file path and read the counts
for file_path in file_paths:
    try:
        # Attempt to read the file into a DataFrame
        current_df = pd.read_csv(file_path, sep='\t', index_col=0)

        # Check if 'NumReads' column exists
        if 'NumReads' in current_df.columns:
            # Extract the sample name from the file path
            sample_name = file_path.split("/")[-2].split("_quant")[0]
            
            # Add the counts as a new column in the DataFrame
            counts[sample_name] = current_df["NumReads"]
            
            files_with_numreads.append(file_path)
        else:
            print(f"'NumReads' column not found in {file_path}.")
            files_without_numreads.append(file_path)

    except Exception as e:
        print(f"Error processing {file_path}: {e}")

# Round the counts and convert to integer
counts = counts.round().astype(int)

# Save the combined counts to a TSV file
counts.to_csv(output_file, sep='\t')

# Print summary information
print(f"Files with 'NumReads': {len(files_with_numreads)}")
for file in files_with_numreads:
    print(f"  - {file}")

print(f"Files without 'NumReads': {len(files_without_numreads)}")
for file in files_without_numreads:
    print(f"  - {file}")
