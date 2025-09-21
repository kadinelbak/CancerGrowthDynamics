# CancerGrowthDynamics

# Processing the data. Most of the files for processing data need to be moved to same directory as datasets in order to make them run. So please move them from the processing data folder to the root directory and then run each of the processing files you wish to run and make to adjust them to your individual well needs and setup. 

# Personal Notes

Please create a file that will navigate to datasets to untreated monoculture and then separate it into 30k and 20k directories based on the well label. It should be so that if the well has the label A1-A3 it gets separated into 20k -> A2780Naive. A4-A6 should be 20k->A2780cis.csv.  whereas if it has label B1-B3 it goes to 30k -> A2780Naive where as label B4-B6 goes to A2780cis. Please ensure you do not remove the original csv. and just make new ones please. 