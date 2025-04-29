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


def calculate_tm(sequence):
    # Placeholder for actual Tm calculation logic
    return 75.0  # Example Tm value, replace with actual calculation
def calculate_homodimer(sequence):
    # Placeholder for actual homodimer calculation logic
    return 75.0, 0.0, sequence  # Example values, replace with actual calculation
def calculate_hairpin(sequence):
    # Placeholder for actual hairpin calculation logic
    return 75.0, 0.0  # Example values, replace with actual calculation
