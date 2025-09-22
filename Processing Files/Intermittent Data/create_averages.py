import pandas as pd
import os
import re
from collections import defaultdict

def parse_filename_components(filename):
    """Extract day, tile, and well information from image filename"""
    # Example: A2780cis_20and30_IC50treatedF_Day10_25_FITC_Tile-0_C4.tif
    
    # Extract day number
    day_match = re.search(r'Day(\d+)', filename)
    day = int(day_match.group(1)) if day_match else None
    
    # Extract tile number  
    tile_match = re.search(r'Tile-(\d+)', filename)
    tile = int(tile_match.group(1)) if tile_match else None
    
    # Extract well (last part before .tif)
    well_match = re.search(r'_([A-Z]\d+)\.tif', filename)
    well = well_match.group(1) if well_match else None
    
    return day, tile, well

def create_day_averages(df, output_path, file_prefix):
    """Create CSV files with averages for each day-well combination"""
    
    # Parse filename components
    df['Day'] = df['Image'].apply(lambda x: parse_filename_components(x)[0])
    df['Tile'] = df['Image'].apply(lambda x: parse_filename_components(x)[1])
    df['Well'] = df['Image'].apply(lambda x: parse_filename_components(x)[2])
    
    # Group by Day and Well, calculate mean Cells across tiles.
    # Inputs are expected to be already converted to Cells in Processed_Datasets.
    value_col = 'Cells' if 'Cells' in df.columns else 'Area µm^2'
    day_averages = df.groupby(['Day', 'Well'])[value_col].agg(['mean', 'count']).reset_index()
    day_averages.columns = ['Day', 'Well', 'Average_Cells', 'Tile_Count']
    
    # Sort by Day and Well
    day_averages = day_averages.sort_values(['Day', 'Well'])
    
    # Create only the combined file with all days
    combined_file = os.path.join(output_path, f"{file_prefix}_sampleaverages.csv")
    day_averages.to_csv(combined_file, index=False)
    print(f"Created day averages file: {combined_file} ({len(day_averages)} records)")

def create_sample_averages(df, output_path, file_prefix):
    """Create CSV files with averages for each well across all days/tiles"""
    
    # Parse filename components
    df['Day'] = df['Image'].apply(lambda x: parse_filename_components(x)[0])
    df['Tile'] = df['Image'].apply(lambda x: parse_filename_components(x)[1])
    df['Well'] = df['Image'].apply(lambda x: parse_filename_components(x)[2])
    
    # Group by Well only, calculate statistics across all tiles and days
    value_col = 'Cells' if 'Cells' in df.columns else 'Area µm^2'
    sample_averages = df.groupby('Well')[value_col].agg([
        'mean', 'std', 'count', 'min', 'max'
    ]).reset_index()
    
    sample_averages.columns = [
        'Well', 
        'Average_Cells', 
        'StdDev_Cells', 
        'Total_Measurements', 
        'Min_Cells', 
        'Max_Cells'
    ]
    
    # Round numerical values for readability
    numeric_cols = ['Average_Cells', 'StdDev_Cells', 'Min_Cells', 'Max_Cells']
    sample_averages[numeric_cols] = sample_averages[numeric_cols].round(2)
    
    # Sort by well
    sample_averages = sample_averages.sort_values('Well')
    
    # Save to CSV
    output_file = os.path.join(output_path, f"{file_prefix}_sampleaverages.csv")
    sample_averages.to_csv(output_file, index=False)
    print(f"Created sample averages file: {output_file} ({len(sample_averages)} wells)")

def process_seeding_density_folder(base_path, density_folder):
    """Process all CSV files in a seeding density folder"""
    
    density_path = os.path.join(base_path, density_folder)
    averages_path = os.path.join(density_path, "Averages")
    
    # Create the Averages directory
    os.makedirs(averages_path, exist_ok=True)
    
    print(f"\n=== Processing {density_folder} ===")
    
    # Find all CSV files in the density folder (excluding subdirectories)
    csv_files = [f for f in os.listdir(density_path) 
                 if f.endswith('.csv') and os.path.isfile(os.path.join(density_path, f))]
    
    for csv_file in csv_files:
        print(f"\nProcessing {csv_file}...")
        
        # Read the CSV file
        file_path = os.path.join(density_path, csv_file)
        df = pd.read_csv(file_path)
        
        # Create file prefix (remove .csv extension)
        file_prefix = csv_file.replace('.csv', '')
        
        # Create day averages
        create_day_averages(df.copy(), averages_path, file_prefix)

# Main execution
def main():
    base_path = r"C:\Users\MainFrameTower\Desktop\CancerGrowthDynamics\Processed_Datasets\Intermittent Data"
    
    # Process both seeding density folders
    for density_folder in ["20k_seeding_density", "30k_seeding_density"]:
        process_seeding_density_folder(base_path, density_folder)
    
    print("\n=== Processing Complete! ===")
    print("\nGenerated files:")
    print("• Day averages: Combined file with average area per well for each day")
    print("• Files are organized in 'Averages' folder within each seeding density directory")

if __name__ == "__main__":
    main()