""" This script collects the results of execution of number of different scenarios and generates insights.
"""
import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 10,           # Set font size to 10
    'font.family': 'serif',    # Use serif fonts
    'font.serif': ['Times New Roman']  # Set font to Times New Roman
})
figure_format = 'pdf'

def read_data_file(input_fp):
    """
    Reads input data file and extracts vertex locations and scores.
    Assumes the file contains 'x y Si' data where:
        x  = x coordinate
        y  = y coordinate
        Si = score

    Returns:
        location (list): List of (x, y) coordinates.
        Si (list): List of corresponding scores.
    """
    if not os.path.exists(input_fp):
        raise FileNotFoundError(f"Input file '{input_fp}' not found.")

    location = []
    Si = []

    with open(input_fp, 'r') as file:
        next(file)  # Skip the first line (header)
        for line in file:
            values = list(map(float, line.split()))
            x, y, score = values[-3], values[-2], values[-1]
            location.append([x, y])
            Si.append(score)

    return location, Si

file_dir = os.path.dirname(os.path.abspath(__file__))
result_dir = os.path.join(file_dir, f'results_for_paper_set2')
plot_dir = os.path.join(result_dir, f'plots')

science_case_list = [('Bioassessment_24_scival', 24, 'Bioassessment'),
                        ('Fracking_24_scival', 24, 'Fracking'),
                        ('Plume_24_scival', 24, 'Low flow estimation'),
                        ('RivNetworkContinuity_24_scival', 24, 'River network continuity')]

# Collect scores for all science cases
scores_all_cases = []
labels = []

for sci_idx, (sci_case_name, N, science_case_label) in enumerate(science_case_list):

    vertex_data_fp = os.path.join(result_dir, f'{sci_case_name}/vertex_data.opc')

    # read in the input data 
    c_pos, Si = read_data_file(vertex_data_fp)

    scores_all_cases.append(Si)
    labels.append(science_case_label)

    plt.figure(figsize=(4, 3))
    plt.hist(Si, bins=30, edgecolor='black', alpha=0.7)
    plt.text(0.05, 0.8, f'Sum(Scores) = {sum(Si):.2f}', transform=plt.gca().transAxes, fontsize=10, color='black')

    #plt.title(f'{science_case_label}')
    plt.xlim(0,1)
    plt.xlabel('Score')
    plt.ylabel('Frequency')
    
    plt.tight_layout() # Adjust layout to prevent cutting off labels

    plot_fp = os.path.join(plot_dir, f'{sci_case_name}_scores.{figure_format}')
    plt.savefig(plot_fp , format=figure_format)
    plt.close()

# Plot a single histogram with different colors for each science case
plt.figure(figsize=(10, 6))
plt.hist(scores_all_cases, bins=30, edgecolor='black', alpha=0.7, label=labels, stacked=True)

# Add title and labels
plt.title('Histogram of Scores for All Science Cases')
plt.xlabel('Scores')
plt.ylabel('Frequency')
plt.legend(title='Science Cases')

# Save the combined plot
plot_fp = os.path.join(plot_dir, 'combined_scores_histogram.{figure_format}')
plt.savefig(plot_fp, format=figure_format)
plt.close()