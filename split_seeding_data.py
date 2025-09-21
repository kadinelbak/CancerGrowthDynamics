import pandas as pd
import os

# This script splits CSV files containing image data into separate files based on seeding density (20k and 30k).
# Each CSV file contains an 'Image' column with filenames that include well identifiers (e.g., C4, D5).
# The mapping of wells to seeding densities is predefined, so if you want to change please copy and change the structure
# Must be run in the same directory as the datasets folder

# Define the base path
base_path = r"C:\Users\MainFrameTower\Desktop\CancerGrowthDynamics\Processed_Datasets\Intermittent Data"

# Define the mapping of files to their well configurations
file_well_mapping = {
    "A2780cisT.csv": {
        "20k_wells": ["C4", "C5", "C6"],
        "30k_wells": ["D4", "D5", "D6"]
    },
    "A2780cisUT.csv": {
        "20k_wells": ["A4", "A5", "A6"],
        "30k_wells": ["B4", "B5", "B6"]
    },
    "A2780UT.csv": {
        "20k_wells": ["A1", "A2", "A3"],
        "30k_wells": ["B1", "B2", "B3"]
    },
    "A2780T.csv": {
        "20k_wells": ["C1", "C2", "C3"],
        "30k_wells": ["D1", "D2", "D3"]
    }
}

def extract_well_from_filename(filename):
    """Extract well identifier from the image filename"""
    # Image filenames contain well info like '_C4.tif', '_D5.tif', etc.
    # Split by underscore and get the part before .tif
    parts = filename.split('_')
    for part in parts:
        if part.endswith('.tif'):
            well = part.replace('.tif', '')
            return well
    return None

def split_csv_by_wells(csv_file, well_config):
    """Split a CSV file based on well identifiers"""
    # Read the original CSV
    df = pd.read_csv(os.path.join(base_path, csv_file))
    
    # Extract well information from image filenames
    df['Well'] = df['Image'].apply(extract_well_from_filename)
    
    # Filter for 20k seeding density wells
    df_20k = df[df['Well'].isin(well_config['20k_wells'])].copy()
    df_20k = df_20k.drop('Well', axis=1)  # Remove the helper column
    
    # Filter for 30k seeding density wells
    df_30k = df[df['Well'].isin(well_config['30k_wells'])].copy()
    df_30k = df_30k.drop('Well', axis=1)  # Remove the helper column
    
    # Save to appropriate folders
    output_20k_path = os.path.join(base_path, "20k_seeding_density", csv_file)
    output_30k_path = os.path.join(base_path, "30k_seeding_density", csv_file)
    
    df_20k.to_csv(output_20k_path, index=False)
    df_30k.to_csv(output_30k_path, index=False)
    
    print(f"Processed {csv_file}:")
    print(f"  20k density: {len(df_20k)} rows -> {output_20k_path}")
    print(f"  30k density: {len(df_30k)} rows -> {output_30k_path}")
    print()

# Process each file
for csv_file, well_config in file_well_mapping.items():
    file_path = os.path.join(base_path, csv_file)
    if os.path.exists(file_path):
        split_csv_by_wells(csv_file, well_config)
    else:
        print(f"Warning: {csv_file} not found at {file_path}")

print("Data splitting completed!")