import pandas as pd
import numpy as np
from Bio.Seq import Seq
import re
import biopython
from pybiomart import server
import time
import requests
# This function will output the dataframe of the primer candidates
# The function will take the primer and shift as input
# The function will return the dataframe of the primer candidates
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

get_primer_candidates("rs9462492, rs58318008, rs1421085, rs9939609, rs1121980", 100)
# This will give us the primer candidates, but does not filter any out that will not work