# RUN EVERYTHING FROM LINES 1-868

# Libraries
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate() 

# importing modules 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Data processing
import pandas as pd
# dplyr is not directly available in Python, consider using pandas for similar functionality
import string # using the built-in 'string' library
import re # using the built-in 're' library as an alternative to stringr
# mosaic is not directly available in Python, you may need to look for equivalent packages
import itertools # using itertools as an alternative to purrr


# Graphing
from plotnine import ggplot, aes, geom_line, theme_minimal # Using plotnine for ggplot2 in Python
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots # Use Plotly's subplot functionality instead of patchworklib


# Bioinformatics
from Bio import Entrez, SeqIO
from pybiomart import Dataset # Using pybiomart as an alternative to biomaRt
# import spgs # No direct equivalent in Python, explore alternatives
import primer3


# Deployment
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

# Source functions
# import functions



# Set options (no direct equivalent in Python, handle manually)
np.set_printoptions(suppress=True)



from Bio.Seq import Seq

def generate_reverse_primers(sequence, primer_length=20, shift=1):
    # Reverse complement the sequence
    reverse_sequence = str(Seq(sequence).reverse_complement())

    # Initialize list to store reverse primers
    reverse_primers = []

    # Iterate over the sequence with the specified shift
    for i in range(0, len(reverse_sequence) - primer_length + 1, shift):
        # Extract primer of specified length
        primer = reverse_sequence[i:i + primer_length]
        # Add primer to list
        reverse_primers.append(primer)

    return reverse_primers

from filter import calculate_tm, calculate_homodimer, calculate_hairpin, calculate_dimer
import time



import re
from primer_candidates import extract_substrings, extract_substrings_far
from primer_candidates import all_text_wrangling, all_text_wrangling_reverse

def get_list(i, j, nested_tables):
    k = " ".join(nested_tables[i][j])
    k = k.split(" ")
    return k



from Practical_Script import extract_substrings


from Bio.Seq import Seq
import pandas as pd
import numpy as np
import re
from itertools import groupby

from primer_candidates import extract_substrings_far


from primer_candidates import all_text_wrangling

def all_text_wrangling_reverse(snp_wrangled, start_distance, end_distance, center, far, shift):
    grouped_sequences = (
        snp_wrangled
        .groupby(['snpID'])
        .agg(sequence_list=('sequence', list))
        .assign(substrings=lambda df: df.sequence_list.apply(lambda x: extract_substrings(x, center, start_distance, end_distance)))
        .explode('substrings')
    )

    grouped_sequences_far = (
        snp_wrangled
        .groupby(['snpID'])
        .head(1)
        .reset_index(drop=True)
        .assign(substrings=lambda df: df.sequence.apply(lambda x: extract_substrings_far(x, center, start_distance, end_distance, far, shift)))
        .explode('substrings')
    )

    grouped_sequences['faraway'] = grouped_sequences_far['substrings']
    grouped_sequences.drop(columns=['sequence_list'], inplace=True)

    grouped_sequences['direction'] = grouped_sequences.duplicated(subset=['snpID'])

    grouped_sequences['direction'] = np.where(grouped_sequences['direction'], "RIGHT", "LEFT")

    return grouped_sequences













# This will give us the primer candidates, but does not filter any out that will not work
# This function will take in the primer and shift
from pybiomart import Dataset
from primer_candidates import get_primer_candidates


import pandas as pd

from filter import get_filter

from multiplex import get_multiplex
#This is the first function that peices together all the functions above
#This function will take in the primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, and hairpin
#This function will output the final results of the dataframe, melting temperature of the primers, and top
def find_acorn(primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, hairpin):
    # get primer candidates is on line ~ 557
    df = get_primer_candidates(primer, shift)
    df = get_filter(df, desired_tm, diff, Heterodimer_tm, Homodimer, hairpin)
    result = get_multiplex(df, Heterodimer_tm, top)
    return result

def main():

    # Establish the variables below before running any functions further in the script!
    primer = "rs9462492, rs58318008, rs1421085, rs9939609, rs1121980"
    shift = 100
    desired_tm = 64
    diff = 3
    Heterodimer_tm = 50
    Homodimer = 45
    top = 2
    hairpin = 45

# RUN EVERYTHING
# this will give us the end results
# Run the functions
# output should be the dataframe?, the melting temperature of the primers, and top?
    result = find_acorn(primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, hairpin)
    print(result)
if __name__ == "__main__":
    main()

