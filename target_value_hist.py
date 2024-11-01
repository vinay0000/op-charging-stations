""" This script collects the results of execution of number of different scenarios and generates insights.
"""
import os
import pickle
import pandas as pd
import numpy as np

import utils

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
    c_pos, Si = utils.read_data_file(vertex_data_fp)

    scores_all_cases.append(Si)
    labels.append(science_case_label)

    utils.plt.figure(figsize=(4, 3))
    utils.plt.hist(Si, bins=30, edgecolor='black', alpha=0.7)
    utils.plt.text(0.05, 0.8, f'Sum(Scores) = {sum(Si):.2f}', transform=utils.plt.gca().transAxes, fontsize=10, color='black')

    #utils.plt.title(f'{science_case_label}')
    utils.plt.xlim(0,1)
    utils.plt.xlabel('Score')
    utils.plt.ylabel('Frequency')
    
    utils.plt.tight_layout() # Adjust layout to prevent cutting off labels

    plot_fp = os.path.join(plot_dir, f'{sci_case_name}_scores.{utils.figure_format}')
    utils.plt.savefig(plot_fp , format=utils.figure_format)
    utils.plt.close()

# Plot a single histogram with different colors for each science case
utils.plt.figure(figsize=(10, 6))
utils.plt.hist(scores_all_cases, bins=30, edgecolor='black', alpha=0.7, label=labels, stacked=True)

# Add title and labels
utils.plt.title('Histogram of Scores for All Science Cases')
utils.plt.xlabel('Scores')
utils.plt.ylabel('Frequency')
utils.plt.legend(title='Science Cases')

# Save the combined plot
plot_fp = os.path.join(plot_dir, 'combined_scores_histogram.{utils.figure_format}')
utils.plt.savefig(plot_fp, format=utils.figure_format)
utils.plt.close()