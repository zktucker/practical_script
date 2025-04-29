from filter import calculate_tm, calculate_dimer
import pandas as pd
import numpy as np
import time
import re

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


def extract_top_n(nested_list, n):
    modified_list = [inner_list[:n] if len(inner_list) >= n else inner_list for inner_list in nested_list]
    return modified_list

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
                level3 = iterate_remove(level3, bad_nodes)
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
