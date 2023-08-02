#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load the gene_presence_absence.Rtab file into a pandas dataframe
rtab_path = open(sys.argv[1])
df = pd.read_csv(rtab_path, sep='\t', index_col=0)

# Count the number of genes present in each genome
counts = df.sum(axis=0)

# Sort the counts in descending order
counts = counts.sort_values(ascending=False)

# Plot the number of genes against the number of genomes
#plt.plot(np.arange(1, len(counts) + 1), counts.values)
#plt.savefig('pangenome.png')

# Define the power law function to fit to the data
def power_law(x, a, b):
    return a * np.power(x, -b)

# Fit the power law function to the data
popt, pcov = curve_fit(power_law, np.arange(1, len(counts) + 1), counts.values)

# Print the Heaps Alpha value
print('Heaps Alpha value:', popt[1])
