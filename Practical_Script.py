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

def get_strong1(x, type):
    temp = ""
    if type:
        target = x[2]
    else:
        target = x[-3]

    if target == "A":
        temp = "G"
    elif target == "G":
        temp = "A"
    elif target == "C":
        temp = "T"
    elif target == "T":
        temp = "C"

    if type:
        x = x[:2] + temp + x[3:]
    else:
        x = x[:-3] + temp + x[-2:]

    return x

def get_strong2(x, type):
    temp = ""

    if type:
        target = x[2]
    else:
        target = x[-3]

    if target == "T":
        temp = "T"
        if type:
            x = x[:2] + temp + x[3:]
        else:
            x = x[:-3] + temp + x[-2:]
        return x
    else:
        return "N"

def get_medium1(x, type):
    temp = ""
    if type:
        target = x[2]
    else:
        target = x[-3]

    if target == "A":
        temp = "A"
    elif target == "G":
        temp = "G"
    elif target == "C":
        temp = "C"
    else:
        return "N"

    if type:
        x = x[:2] + temp + x[3:]
    else:
        x = x[:-3] + temp + x[-2:]

    return x

def get_weak1(x, type):
    temp = ""
    if type:
        target = x[2]
    else:
        target = x[-3]

    if target == "C":
        temp = "A"
    elif target == "A":
        temp = "C"
    elif target == "G":
        temp = "T"
    elif target == "T":
        temp = "G"

    if type:
        x = x[:2] + temp + x[3:]
    else:
        x = x[:-3] + temp + x[-2:]
    return x



import re

def list_seq(snp):
    first_position = re.search('/', snp).start()
    number_slash = snp.count('/')
    block = snp[first_position - 1: first_position - 1 + number_slash * 2 + 1]
    block = block.replace("/", "")
    block = list(block)

    k = []
    for i in block:
        k.append(snp[:first_position - 1] + i + snp[first_position - 1 + number_slash * 2 + 1:])
    
    k = [seq.replace("%", "") for seq in k]
    k = [seq for seq in k if "W" not in seq]
    return k

def get_list(i, j, nested_tables):
    k = " ".join(nested_tables[i][j])
    k = k.split(" ")
    return k

def get_endpoints(lst, current_name="", parent_names=[]):
    endpoints = []

    if isinstance(lst, dict):
        if len(lst) > 0:
            for nested_name, nested_value in lst.items():
                if isinstance(nested_value, dict):
                    nested_endpoints = get_endpoints(nested_value, f"{current_name}/{nested_name}", parent_names + [current_name])
                    endpoints.extend(nested_endpoints)
                else:
                    endpoint = {"endpoint": nested_name, "parents": parent_names + [current_name]}
                    endpoints.append(endpoint)
        else:
            endpoint = {"endpoint": current_name, "parents": parent_names}
            endpoints.append(endpoint)
    else:
        endpoint = {"endpoint": current_name, "parents": parent_names}
        endpoints.append(endpoint)
    
    return endpoints

def clean_endpoints(endpoints):
    for endpoint in endpoints:
        if len(endpoint['parents']) > 0:
            endpoint['parents'] = endpoint['parents'][1:]
        
        for j in range(len(endpoint['parents'])):
            split_string = endpoint['parents'][j].split("/")
            desired_item = split_string[-1]
            endpoint['parents'][j] = desired_item

    return endpoints

def compute_bad_nodes(endpoints, threshold):
    blacklist = []

    for endpoint in endpoints:
        result = sum(calculate_dimer(endpoint['endpoint'], parent)['temp'] > threshold for parent in endpoint['parents'])
        blacklist.append(result)
    
    bad_nodes = [endpoints[i] for i in range(len(endpoints)) if blacklist[i] == 1]
    return bad_nodes


def remove_list(lst, path):
    if len(path) == 1:
        if isinstance(lst, dict) and path[0] in lst:
            lst.pop(path[0], None)
    else:
        if isinstance(lst, dict) and path[0] in lst:
            lst[path[0]] = remove_list(lst[path[0]], path[1:])
            if isinstance(lst[path[0]], dict) and len(lst[path[0]]) == 0 and not any(lst[path[0]]):
                lst.pop(path[0], None)
    return lst

def remove_empty_lists(lst):
    if isinstance(lst, list) or isinstance(lst, dict):
        lst = {k: remove_empty_lists(v) for k, v in lst.items() if v}
    return lst

def iterate_remove(level3, bad_nodes):
    for node in bad_nodes:
        level3 = remove_list(level3, node['parents'] + [node['endpoint']])
    return level3

def incoming_list(arranged_list):
    level4 = {}
    for item in arranged_list:
        # Create a sublist with the name as the item
        level4[item] = {item: 1}
    return level4

def replace_end_nodes(lst, replace_lst):
    if isinstance(lst, dict):
        if len(lst) == 0:
            return replace_lst
        else:
            return {k: replace_end_nodes(v, replace_lst) for k, v in lst.items()}
    else:
        return replace_lst

def extract_substrings(string, center, start_distance, end_distance):
    # Empty lists to store substrings
    substrings_left = []
    substrings_right = []

    for item in string:
        # Right flanking
        for distance in range(start_distance, end_distance + 1):
            sub = item[center + 1:center + 1 + distance]
            substrings_right.extend([
                str(Seq(get_strong1(sub, True)).reverse_complement()).upper(),
                str(Seq(get_strong2(sub, True)).reverse_complement()).upper(),
                str(Seq(get_medium1(sub, True)).reverse_complement()).upper(),
                str(Seq(get_weak1(sub, True)).reverse_complement()).upper()
            ])

        # Left flanking
        for distance in range(start_distance, end_distance + 1):
            sub = item[center - distance:center + 1]
            substrings_left.extend([
                get_strong1(sub, False),
                get_strong2(sub, False),
                get_medium1(sub, False),
                get_weak1(sub, False)
            ])

        start_distance += 1
        end_distance += 1

    # Return the extracted substrings
    return {
        "left": [s for s in substrings_left if s != "N"],
        "right": [s for s in substrings_right if s != "N"]
    }


from Bio.Seq import Seq
import pandas as pd
import numpy as np
import re
from itertools import groupby

def extract_substrings_far(string, center, start_distance, end_distance, far, shift):
    # Empty lists to store substrings
    substrings_left = []
    substrings_right = []

    for i in range(far, far + shift):
        # to Right
        for distance in range(start_distance, end_distance + 1):
            sub = string[i + center:i + center + distance + 1]
            substrings_right.append(str(Seq(sub).reverse_complement()).upper())

        # Left flanking
        for distance in range(start_distance, end_distance + 1):
            substrings_left.append(string[center - distance - i:center - i + 1])

    # Return the extracted substrings
    return {'right': substrings_right, 'left': substrings_left}

def all_text_wrangling(snp_wrangled, start_distance, end_distance, center, far, shift):
    grouped_sequences = (
        snp_wrangled
        .groupby(['snpID', 'snp_character'])
        .agg(sequence_list=('sequence', list))
        .assign(substrings=lambda df: df.sequence_list.apply(lambda x: extract_substrings(x, center, start_distance, end_distance)))
        .explode('substrings')
    )

    grouped_sequences_far = (
        snp_wrangled
        .groupby(['snpID', 'snp_character'])
        .head(1)
        .reset_index(drop=True)
        .assign(substrings=lambda df: df.sequence.apply(lambda x: extract_substrings_far(x, center, start_distance, end_distance, far, shift)))
        .explode('substrings')
    )

    grouped_sequences['faraway'] = grouped_sequences_far['substrings']
    grouped_sequences.drop(columns=['sequence_list'], inplace=True)

    # what is this code below for?
    grouped_sequences['direction'] = grouped_sequences.duplicated(subset=['snpID', 'snp_character'])

    grouped_sequences['direction'] = np.where(grouped_sequences['direction'], "RIGHT", "LEFT")

    return grouped_sequences

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


def stage1_filter(df, desired_tm, diff, Homodimer, hairpin):
    for i in range(len(df.iloc[:, 2])):  # Adjusting for Python indexing

        # Homodimer
        k = [seq for seq in df.iloc[i, 2] if calculate_homodimer(seq)[1] < Homodimer]
        if len(k) > 5:
            df.iloc[i, 2] = [seq for seq in df.iloc[i, 2] if calculate_homodimer(seq)[2] < Homodimer]
        else:
            print(f"Homodimer - Bottle neck {df.iloc[i, 0]}")
            calculated_values = [calculate_homodimer(seq) for seq in df.iloc[i, 2]]
            differences = [abs(val[1] - Homodimer) for val in calculated_values]
            min_diff_indices = sorted(range(len(differences)), key=lambda x: differences[x])[:5]
            df.iloc[i, 2] = [df.iloc[i, 2][index] for index in min_diff_indices]

        # Hairpin
        k = [seq for seq in df.iloc[i, 2] if calculate_hairpin(seq)[1] < hairpin]
        if len(k) > 5:
            df.iloc[i, 2] = [seq for seq in df.iloc[i, 2] if calculate_hairpin(seq)[1] < hairpin]
        else:
            print(f"Hairpin - Bottle neck {df.iloc[i, 0]}")
            calculated_values = [calculate_hairpin(seq) for seq in df.iloc[i, 2]]
            differences = [abs(val[1] - hairpin) for val in calculated_values]
            min_diff_indices = sorted(range(len(differences)), key=lambda x: differences[x])[:5]
            df.iloc[i, 2] = [df.iloc[i, 2][index] for index in min_diff_indices]

        # Filter Tm above target
        k = [seq for seq in df.iloc[i, 2] if calculate_tm(seq) < desired_tm + diff]
        if len(k) > 5:
            df.iloc[i, 2] = k
        else:
            print(f"Tm_above - Bottle neck {df.iloc[i, 0]}")
            calculated_values = [calculate_tm(seq) for seq in df.iloc[i, 2]]
            differences = [abs(val - (desired_tm + diff)) for val in calculated_values]
            min_diff_indices = sorted(range(len(differences)), key=lambda x: differences[x])[:5]
            df.iloc[i, 2] = [df.iloc[i, 2][index] for index in min_diff_indices]

        # Filter Tm below target
        k = [seq for seq in df.iloc[i, 2] if calculate_tm(seq) > desired_tm - diff]
        if len(k) > 5:
            df.iloc[i, 2] = k
        else:
            print(f"TM below - Bottle neck {df.iloc[i, 0]}")
            calculated_values = [calculate_tm(seq) for seq in df.iloc[i, 2]]
            differences = [abs(val - (desired_tm - diff)) for val in calculated_values]
            min_diff_indices = sorted(range(len(differences)), key=lambda x: differences[x])[:5]
            df.iloc[i, 2] = [df.iloc[i, 2][index] for index in min_diff_indices]

    for i in range(len(df.iloc[:, 3])):
        if len(df.iloc[i, 3]) != 0:
            df.iloc[i, 3] = [seq for seq in df.iloc[i, 3] if calculate_tm(seq) > desired_tm - diff]
        if len(df.iloc[i, 3]) != 0:
            df.iloc[i, 3] = [seq for seq in df.iloc[i, 3] if calculate_hairpin(seq)[1] < hairpin]

    for i in range(len(df.iloc[:, 3])):
        if len(df.iloc[i, 3]) != 0:
            df.iloc[i, 3] = [seq for seq in df.iloc[i, 3] if calculate_homodimer(seq)[1] < Homodimer]
            df.iloc[i, 3] = [seq for seq in df.iloc[i, 3] if calculate_tm(seq) < desired_tm + diff]

    for i in range(len(df) - 1, -1, -1):
        if len(df.iloc[i, 3]) == 0:
            df = df.drop(i)

    return df


def extract_top_n(nested_list, n):
    modified_list = [inner_list[:n] if len(inner_list) >= n else inner_list for inner_list in nested_list]
    return modified_list

def get_display_tree(level3, keep):
    endpoints = get_endpoints(level3)

    # Endpoints come back a little messy
    endpoints = clean_endpoints(endpoints)

    display_tree = []
    for i in range(keep):
        display_tree.append(endpoints[i])

    display_tree = pd.DataFrame(display_tree)

    first_row = display_tree.iloc[0]
    display_tree = display_tree.iloc[1:]

    display_tree = pd.concat([display_tree, first_row.to_frame().T])

    display_tree.columns = [f"Option {i + 1}" for i in range(keep)]

    return display_tree

def soulofmultiplex(df, Heterodimer_tm):
    list_3 = []

    for i in range(len(df.iloc[:, 0])):
        list_3.append(df.iloc[i, 4])
        list_3.append(df.iloc[i, 5])

    # Arrange the list from small to big
    arranged_list = list_3

    # Prepare the initial list for multiplexing
    level2 = incoming_list(arranged_list[0])
    level3 = replace_end_nodes(incoming_list(arranged_list[0]), incoming_list(arranged_list[1]))

    if len(arranged_list) != 2:
        level3 = replace_end_nodes(level3, incoming_list(arranged_list[2]))
        print(len(arranged_list))

        for i in range(3, len(arranged_list)):
            start_time = time.time()

            # Get all the end points from the tree
            endpoints = get_endpoints(level3)

            # Endpoints come back a little messy
            endpoints = clean_endpoints(endpoints)
            print(f"Start with {len(endpoints)}")

            # Evaluate all the end points to their parents
            bad_nodes = compute_bad_nodes(endpoints, Heterodimer_tm)
            print(f"We are removing: {len(bad_nodes)}")

            # Remove bad nodes if there are any
            if len(bad_nodes) != 0:
                level3 = Iterate_remove(level3, bad_nodes)
                level3 = remove_empty_lists(level3)

            # If all nodes are bad, return None
            if len(endpoints) == len(bad_nodes):
                print("All nodes are removed during the process")
                return None

            print(f"After trimming: {len(get_endpoints(level3))}")

            # Stop adding list if we are at the last level
            level4 = incoming_list(arranged_list[i])
            print(f"New list: {len(level4)}")

            level3 = replace_end_nodes(level3, level4)
            print(f"level3 + level4: {len(get_endpoints(level3))}")

            # Summarize results for this level
            print(f"How far are we: {i}")
            print(f"Time {round(time.time() - start_time, 1)}")
            print("--------------------------")

    level5 = get_display_tree(level3, 3)
    repeated_list = [item for sublist in [[val] * 2 for val in df.iloc[:, 0]] for item in sublist]
    suffix = ["_forward", "_reverse"]
    modified_list = [f"{repeated_list[i]}{suffix[i % len(suffix)]}" for i in range(len(repeated_list))]
    level5.index = modified_list

    level5.loc['Tm'] = [round(np.mean([calculate_tm(x) for x in col]), 2) for col in level5]

    return level5


def get_tm_for_all_primers(level5):
    import pandas as pd
    import numpy as np

    level5_with_tm_result = pd.DataFrame(index=level5.index)

    # Apply the 'calculate_tm' function to each column of the dataframe
    for col in level5.columns:
        # Calculate TM for the column and round the result
        tm_results = [round(calculate_tm(seq), 2) for seq in level5[col]]

        # Combine the original column with the TM results
        combined = pd.DataFrame({col: level5[col], f"{col}_tm": tm_results})

        # Bind the new combined columns to the result dataframe
        level5_with_tm_result = pd.concat([level5_with_tm_result, combined], axis=1)

    # Remove the first column if it contains only NA values from the placeholder creation
    level5_with_tm_result = level5_with_tm_result.dropna(axis=1, how='all')
    level5_with_tm_result.index = level5.index
    level5_with_tm_result = level5_with_tm_result.to_numpy()
    return level5_with_tm_result

# This will give us the primer candidates, but does not filter any out that will not work
# This function will take in the primer and shift
# This function will output the dataframe of the primer candidates
def get_primer_candidates(primer, shift):
    from pybiomart import Dataset

    # Explore options 800 bp away from the SNP location upstream and downstream
    center = 800
    hairpin = 45
    far = 200
    start_distance = 15
    end_distance = 28

    # Accessing database
    print("Execute MART API")
    snp_list = primer.split(" ")
    upStream = center
    downStream = center

    snpmart = Dataset(name="hsapiens_snp", host="http://www.ensembl.org")
    snp_sequence = snpmart.query(attributes=['refsnp_id', 'snp'],
                                 filters={'snp_filter': snp_list, 'upstream_flank': upStream, 'downstream_flank': downStream})

    # Create a new data frame (df) to store the snp variation sequences
    snp_wrangled = pd.DataFrame(columns=['snpID', 'sequence'])

    # New dataframe contains all snp variation sequences plus a snp id column
    for j in snp_sequence['refsnp_id']:
        for i in list_seq(snp_sequence[snp_sequence['refsnp_id'] == j]['snp'].values[0]):
            snp_wrangled = snp_wrangled.append({'snpID': j, 'sequence': i}, ignore_index=True)

    # Rename columns and data frame
    snp_wrangled['snp_character'] = snp_wrangled['sequence'].str[800]

    # I have a long long string. I want to get the left 18~25 characters and
    # between 300 ~ 800 units away, I want another 18 ~ 25
    df = all_text_wrangling(snp_wrangled, start_distance, end_distance, center, far, shift)

    print("Primer generated")
    return df


import pandas as pd

def get_filter(df, desired_tm, diff, Heterodimer_tm, Homodimer, hairpin):
    print("Python get filter activated")
    # Applied filters before multiplexing
    df = stage1_filter(df, desired_tm, diff, Homodimer, hairpin)
    print(df)

    print("Filtered")

    # Count how many candidates there are for each primer group
    df['substrings_count'] = df['substrings'].apply(len)
    df['faraway_count'] = df['faraway'].apply(len)
    df = df[['snpID', 'substrings_count', 'faraway_count'] + [col for col in df.columns if col not in ['snpID', 'substrings_count', 'faraway_count']]]

    # Display the updated nested dataframe
    return df

def get_multiplex(df, Heterodimer_tm, top):
    print("Tree search")
    df
    # Keep only certain amount of candidates
    df.iloc[:, 4] = extract_top_n(df.iloc[:, 4], 2)
    df.iloc[:, 5] = extract_top_n(df.iloc[:, 5], 2)

    # Technical debt
    df_rm = df.drop_duplicates(subset=['snpID'])

    df_rm = df.groupby('snpID').apply(lambda x: x[x['substrings_count'] == x['substrings_count'].max()]).reset_index(drop=True)
    print(df_rm)

    level5 = soulofmultiplex(df_rm, Heterodimer_tm)
    print(level5)

    level5_with_tm_result = get_tm_for_all_primers(level5)
    return level5_with_tm_result

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

