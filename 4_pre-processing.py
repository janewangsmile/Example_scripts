import os
os.chdir("/home/jun/Desktop/lncRNAs_ASD/")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product

#################################################
# Helper functions are defined below.
#################################################
# paths
data_path = "Dataset/feature_selected/"

# data set file names
express_file = "1_expression_with_category3_genes.csv"

label_file = "2_labels.csv"
norm_express_file = "2_expression_normalized.csv"


def process_expression():
    """
    The function reads the expression dataset, ignores string values, and normalize the floating
    point numbers based on log_2(number+1).
    """
    with open(data_path + express_file, "r") as infile:
        ncols = len(infile.readline().split(','))
    
    print(ncols)
    # Skip first 3 columns
    #selected feature:
    features = np.genfromtxt(data_path + express_file, delimiter=',', usecols=range(3, ncols), skip_header=1)
    # full feature:
    features = np.genfromtxt(data_path + express_file, delimiter=',', usecols=range(3, ncols-1), skip_header=1)

    # Normalization using logarithm
    features = np.log2(features + 1) 

    # Save the normalized data to file
    np.savetxt(data_path + norm_express_file, features, delimiter=',', fmt='%.6f')

def process_label():
    """
    The function reads the expression data and extracts the label information, i.e., ASD will be 
    treated as positive and other diseases.
    """
    tmp = np.genfromtxt(data_path + express_file, dtype=None, delimiter=',', usecols=2, skip_header=1, encoding="ascii")
    my_dict = {'"ASD"': 1, '"Disease"': -1}
    labels = [my_dict[i] for i in tmp]
    np.savetxt(data_path + label_file, labels, fmt="%i", delimiter=',')

#################################################
# Main entry starts here.
#################################################            
if __name__ == "__main__":
    process_expression()
    process_label()