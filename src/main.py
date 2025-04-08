############################

from collections import Counter
from itertools import combinations
import math

def count_pairs(sequences):
    pair_counts = Counter()
    length = len(sequences[0])  # sequence length
    for col in range(length):
        column = [seq[col] for seq in sequences]
        for a, b in combinations(column, 2):
            key = tuple(sorted((a, b)))
            pair_counts[key] += 1
    return pair_counts

def compute_frequencies(pair_counts):
    total_pairs = sum(pair_counts.values())
    return {k: v / total_pairs for k, v in pair_counts.items()}

def compute_background_freq(sequences):
    counts = Counter("".join(sequences))
    total = sum(counts.values())
    return {aa: c / total for aa, c in counts.items()}

def build_substitution_matrix(freqs, background_freqs):
    matrix = {}
    for (a, b), observed in freqs.items():
        expected = (
            background_freqs[a] * background_freqs[b]
            if a != b else
            background_freqs[a] * background_freqs[b]
        )
        score = math.log2(observed / expected) if expected > 0 else float('-inf')
        matrix[(a, b)] = round(score, 2)
    return matrix

def read_sequences(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

############################

def parse_labeled_sequences(file_path):
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if line.strip()]
        seqs = lines[::2]
        labels = lines[1::2]
    return seqs, labels

def count_emissions(seqs, labels):
    emissions = {0: Counter(), 1: Counter(), 2: Counter()}
    for seq, label in zip(seqs, labels):
        for param, state in zip(seq, label):
            emissions[int(state)][param] += 1
    return emissions

def normalize_emissions(emissions):
    emission_probs = {}
    for state, counts in emissions.items():
        total = sum(counts.values())
        emission_probs[state] = {aa: round(c / total, 4) for aa, c in counts.items()}
    return emission_probs

def compute_transitions(labels):
    transitions = Counter()
    state_totals = Counter()
    for label in labels:
        for i in range(len(label) - 1):
            a, b = int(label[i]), int(label[i+1])
            transitions[(a, b)] += 1
            state_totals[a] += 1
    # Normalize
    T = {k: round(v / state_totals[k[0]], 4) for k, v in transitions.items()}
    return T

############################

# Print substitution matrix
def print_substitution_matrix(subst_matrix):
    print("Substitution Matrix:")
    print(sub_matrix)

# Print emission probabilities
def print_emission(emission_probs):
   print("Emission Probabilities:")
   print(emission_probs)

# Print transition probabilities 
def print_transition(transition_probs):
    print("Transition Probabilities:")
    for (from_state, to_state), prob in transition_probs.items():
        print(f"  {from_state} -> {to_state}: {prob}")

############################

# Part A
sequence1 = read_sequences("../input/DataFile1-1.txt") # List, each index contains a row of the file
pairs = count_pairs(sequence1)
freqs = compute_frequencies(pairs)
background = compute_background_freq(sequence1)
sub_matrix = build_substitution_matrix(freqs, background)

# Part B
sequence2, labels = parse_labeled_sequences("../input/DataFile2.txt") # List, each index contains a row of the file
emissions = count_emissions(sequence2, labels)
emission_probs = normalize_emissions(emissions)
transition_probs = compute_transitions(labels)

# Output
from rich import print
print_substitution_matrix(sub_matrix)
print_emission(emission_probs)
print_transition(transition_probs)

############################