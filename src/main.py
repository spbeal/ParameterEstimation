############################
# Samuel Beal
# University of Idaho CS415
############################

from rich import print
from collections import Counter
from itertools import combinations
import math

############################

# Creates all the pair keys (i, a), (a, a) etc...
def count_pairs(sequences):
    pair_counts = Counter()
    length = len(sequences[0])  # sequence length
    for col in range(length):
        column = [seq[col] for seq in sequences]
        for a, b in combinations(column, 2):
            key = tuple(sorted((a, b)))
            pair_counts[key] += 1
    return pair_counts

# sums values and divides by number of observed
def compute_frequencies(pair_counts):
    total_pairs = sum(pair_counts.values())
    return {k: v / total_pairs for k, v in pair_counts.items()}

# Sums the number of characters
def compute_background_freq(sequences):
    counts = Counter("".join(sequences))
    total = sum(counts.values())
    return {item: c / total for item, c in counts.items()}

# Log odds matrix and scores
def build_substitution_matrix(freqs, background_freqs):
    matrix = {}
    for (a, b), observed in freqs.items(): # pairs
        expected = background_freqs[a] * background_freqs[b] # number of x * number of y
        score = math.log2(observed / expected) if expected > 0 else float('-inf') # log odds
        matrix[(a, b)] = round(score, 2)
    return matrix

# reads data file 1
def read_sequences(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

############################

# reads data file 2
def parse_labeled_sequences(file_path):
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if line.strip()]
        seqs = lines[::2]
        labels = lines[1::2]
    return seqs, labels

# analyzes every single row of the sequence of characters
def count_emissions(seqs, labels):
    emissions = {0: Counter(), 1: Counter(), 2: Counter()}
    for seq, label in zip(seqs, labels):
        for item, state in zip(seq, label):
            emissions[int(state)][item] += 1
    return emissions

# normalize
def normalize_emissions(emissions):
    emission_probs = {}
    for state, counts in emissions.items():
        total = sum(counts.values())
        emission_probs[state] = {item: round(c / total, 4) for item, c in counts.items()}
    return emission_probs

# Compute transitions based on all the labels 012 in datafile 2
def compute_transitions(labels):
    transitions = Counter() # counts transitions
    state_totals = Counter() # counts states we have been in
    for label in labels:
        for i in range(len(label) - 1):
            a, b = int(label[i]), int(label[i+1]) # 0 to 1
            transitions[(a, b)] += 1
            state_totals[a] += 1
    # Normalize
    T = {k: round(v / state_totals[k[0]], 4) for k, v in transitions.items()}
    return T

############################

# Print substitution matrix
# def print_substitution_matrix(subst_matrix):
#     print("\n=== Substitution Matrix ===")
#     print(sub_matrix)
def print_substitution_matrix(subst_matrix):
    print("\n=== Substitution Matrix ===")

    # Extract unique characters from substitution matrix keys
    chars = sorted(set(a for pair in subst_matrix for a in pair))

    # Print header row
    header = "    " + "  ".join(f"{c:>5}" for c in chars)
    print(header)
    print("    " + "-" * (len(header) - 4))

    # Print each row
    for a in chars:
        row = f"{a:>2} |"
        for b in chars:
            # Sorted tuple ensures symmetric access
            key = tuple(sorted((a, b)))
            value = subst_matrix.get(key, 0.0)
            row += f"{value:6.1f} "
        print(row)


# Print emission probabilities nicely
# def print_emission(emission_probs):
#     print("\n=== Emission Probabilities ===")
#     for state in sorted(emission_probs.keys()):
#         print(f"State {state}:")
#         probs = emission_probs[state]
#         line = ", ".join([f"{k}:{v}" for k, v in sorted(probs.items())])
#         print(f"  {line}")
def print_emission(emission_probs):
    print("\n=== Emission Probability Table ===")

    # Get all unique emission characters across states
    all_symbols = sorted(set(k for state_probs in emission_probs.values() for k in state_probs))

    # Header
    header = "State | " + "  ".join(f"{s:>5}" for s in all_symbols)
    print(header)
    print("-" * len(header))

    # Rows for each state
    for state in sorted(emission_probs):
        row = f"  {state}   | "
        for symbol in all_symbols:
            prob = emission_probs[state].get(symbol, 0.0)
            row += f"{prob:6.3f} "
        print(row)


# Print all possible transition probabilities, including zero probabilities
def print_transition(transition_probs):
    print("\n=== Transition Probabilities ===")
    states = [0, 1, 2]
    for from_state in states:
        for to_state in states:
            prob = transition_probs.get((from_state, to_state), 0.0)
            print(f"{from_state} -> {to_state}: {prob}")


############################

# Part A
sequence1 = read_sequences("../input/DataFile1-1.txt") # List, each index contains a row of the file
pairs = count_pairs(sequence1)
background = compute_background_freq(sequence1)
freqs = compute_frequencies(pairs)
sub_matrix = build_substitution_matrix(freqs, background)

# Part B
sequence2, labels = parse_labeled_sequences("../input/DataFile2.txt") # List, each index contains a row of the file
emissions = count_emissions(sequence2, labels)
emission_probs = normalize_emissions(emissions)
transition_probs = compute_transitions(labels)

# Output
print_substitution_matrix(sub_matrix)
print_emission(emission_probs)
print_transition(transition_probs)

############################