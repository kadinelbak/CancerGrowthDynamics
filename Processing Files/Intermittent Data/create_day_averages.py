#!/usr/bin/env python3
"""
Script to create day average CSV files from sample average CSV files.
Takes sample averages for each day and creates a single data point per day
with error bar statistics (mean, standard deviation, standard error, count).
"""

import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path


def create_day_averages(sample_avg_file_path):
    """
    Create day averages from a sample averages CSV file.
    
    Parameters:
    sample_avg_file_path (str): Path to the sample averages CSV file
    
    Returns:
    pandas.DataFrame: DataFrame with day averages and error statistics
    """
    # Read the sample averages file
    df = pd.read_csv(sample_avg_file_path)
    
    # Group by Day and calculate statistics
    # Support both legacy 'Average_Area_um2' and new 'Average_Cells'
    value_col = 'Average_Cells' if 'Average_Cells' in df.columns else 'Average_Area_um2'
    day_stats = df.groupby('Day')[value_col].agg([
        'mean',      # Average of all wells for each day
        'std',       # Standard deviation
        'count',     # Number of wells/samples
        'sem'        # Standard error of the mean
    ]).reset_index()
    
    # Rename columns for clarity
    day_stats.columns = ['Day', 'Mean_Cells', 'Std_Dev_Cells', 'Sample_Count', 'Std_Error_Cells']
    
    # Calculate 95% confidence interval bounds (using normal approximation for simplicity)
    # For small samples, this is an approximation; proper t-distribution would be better
    confidence_level = 0.95
    z_critical = 1.96  # 95% confidence interval for normal distribution
    
    # Calculate confidence intervals using normal approximation
    day_stats['CI_95_Margin_Cells'] = z_critical * day_stats['Std_Error_Cells']
    day_stats['CI_95_Lower_Cells'] = day_stats['Mean_Cells'] - day_stats['CI_95_Margin_Cells']
    day_stats['CI_95_Upper_Cells'] = day_stats['Mean_Cells'] + day_stats['CI_95_Margin_Cells']
    
    return day_stats


def process_directory(directory_path):
    """
    Process all sample average CSV files in a directory.
    
    Parameters:
    directory_path (str): Path to the directory containing sample average files
    """
    # Find all sample average CSV files
    pattern = os.path.join(directory_path, "*sampleaverages.csv")
    sample_files = glob.glob(pattern)
    
    if not sample_files:
        print(f"No sample average files found in {directory_path}")
        return
    
    print(f"Found {len(sample_files)} sample average files in {directory_path}")
    
    for file_path in sample_files:
        try:
            # Create day averages
            day_averages = create_day_averages(file_path)
            
            # Generate output filename
            file_name = os.path.basename(file_path)
            output_name = file_name.replace("_sampleaverages.csv", "_dayaverages.csv")
            output_path = os.path.join(directory_path, output_name)
            
            # Save the day averages
            day_averages.to_csv(output_path, index=False)
            print(f"Created: {output_path}")
            
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")


def main():
    """
    Main function to process all directories with sample average files.
    """
    base_path = r"c:\Users\MainFrameTower\Desktop\CancerGrowthDynamics\Processed_Datasets\Intermittent Data"
    
    # Process 20k seeding density
    dir_20k = os.path.join(base_path, "20k_seeding_density", "Averages")
    if os.path.exists(dir_20k):
        print("Processing 20k seeding density files...")
        process_directory(dir_20k)
    else:
        print(f"Directory not found: {dir_20k}")
    
    # Process 30k seeding density
    dir_30k = os.path.join(base_path, "30k_seeding_density", "Averages")
    if os.path.exists(dir_30k):
        print("\nProcessing 30k seeding density files...")
        process_directory(dir_30k)
    else:
        print(f"Directory not found: {dir_30k}")


if __name__ == "__main__":
    main()