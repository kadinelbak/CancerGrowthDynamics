#!/usr/bin/env python3
"""
Script to process cancer growth dynamics datasets.
Divides all "Area µm^2" values by 144 for normalization and renames column to "Cells µm^2".

This script processes all CSV files in the Datasets folder and its subfolders,
modifying the "Area µm^2" column by dividing each value by 144 and renaming it to "Cells µm^2".
The processed files are saved to a new "Processed_Datasets" folder,
preserving the original files.

Must be in the same directory as the "Datasets" folder.
"""

import os
import pandas as pd
import glob
from pathlib import Path

def process_csv_file(file_path, output_path):
    """
    Process a single CSV file by dividing the "Area µm^2" column by 144 and renaming to "Cells".
    
    Args:
        file_path (str): Path to the CSV file to process
        output_path (str): Path where the processed file should be saved
        
    Returns:
        bool: True if successful, False if error occurred
    """
    try:
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Check if the required column exists
        if 'Area µm^2' not in df.columns:
            print(f"Warning: 'Area µm^2' column not found in {file_path}")
            return False
        
        # Divide the Area column by 144
        df['Area µm^2'] = df['Area µm^2'] / 144
        
        # Rename the column from "Area µm^2" to "Cells µm^2"
        df = df.rename(columns={'Area µm^2': 'Cells'})
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Save the modified data to the new location
        df.to_csv(output_path, index=False)
        
        print(f"✓ Successfully processed: {os.path.basename(file_path)} → {os.path.basename(output_path)}")
        print(f"  → Normalized values (÷144) and renamed column to 'Cells'")
        return True
        
    except Exception as e:
        print(f"✗ Error processing {file_path}: {str(e)}")
        return False

def main():
    """
    Main function to process all CSV files in the Datasets folder.
    """
    # Get the script directory and construct paths
    script_dir = Path(__file__).parent
    datasets_dir = script_dir / "Datasets"
    processed_dir = script_dir / "Processed_Datasets"
    
    # Check if Datasets directory exists
    if not datasets_dir.exists():
        print(f"Error: Datasets directory not found at {datasets_dir}")
        return
    
    # Create the processed datasets directory
    processed_dir.mkdir(exist_ok=True)
    
    # Find all CSV files in the Datasets folder and subfolders
    csv_pattern = str(datasets_dir / "**" / "*.csv")
    csv_files = glob.glob(csv_pattern, recursive=True)
    
    if not csv_files:
        print("No CSV files found in the Datasets directory.")
        return
    
    print(f"Found {len(csv_files)} CSV files to process...")
    print(f"Original files will be preserved in: {datasets_dir}")
    print(f"Processed files will be saved to: {processed_dir}")
    print("-" * 60)
    
    # Process each CSV file
    successful_count = 0
    failed_count = 0
    
    for csv_file in csv_files:
        # Calculate the relative path from the datasets directory
        rel_path = Path(csv_file).relative_to(datasets_dir)
        
        # Create the corresponding output path in the processed directory
        output_path = processed_dir / rel_path
        
        if process_csv_file(csv_file, str(output_path)):
            successful_count += 1
        else:
            failed_count += 1
    
    # Print summary
    print("-" * 60)
    print(f"Processing complete!")
    print(f"✓ Successfully processed: {successful_count} files")
    if failed_count > 0:
        print(f"✗ Failed to process: {failed_count} files")
    
    print(f"\nOriginal files preserved in: {datasets_dir}")
    print(f"Processed files (Area µm^² ÷ 144 → Cells µm^²) saved to: {processed_dir}")
    print(f"\nAll 'Area µm^2' values normalized (÷144) and column renamed to 'Cells'.")

def preview_changes():
    """
    Preview what changes will be made without actually processing files.
    """
    script_dir = Path(__file__).parent
    datasets_dir = script_dir / "Datasets"
    
    # Find all CSV files
    csv_pattern = str(datasets_dir / "**" / "*.csv")
    csv_files = glob.glob(csv_pattern, recursive=True)
    
    print("Preview of files that will be processed:")
    print("-" * 60)
    
    for csv_file in csv_files:
        rel_path = Path(csv_file).relative_to(datasets_dir)
        print(f"Source:      Datasets/{rel_path}")
        print(f"Destination: Processed_Datasets/{rel_path}")
        print()
    
    print(f"Total files to process: {len(csv_files)}")

if __name__ == "__main__":
    print("Cancer Growth Dynamics - Area Data Processor")
    print("=" * 50)
    
    # Show preview
    preview_changes()
    
    # Confirm before proceeding
    response = input("\nProceed with processing? This will create 'Processed_Datasets' folder with normalized cell values and renamed columns (y/n): ").lower().strip()
    
    if response in ['y', 'yes']:
        main()
    else:
        print("Operation cancelled.")