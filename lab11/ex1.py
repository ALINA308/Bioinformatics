import numpy as np
import pandas as pd

# 1. Input Sequences
s1 = "ATCGATTCGATATCATACACGTAT"  # CpG Island (+)
s2 = "CTCGACTAGTATGAAGTCCACGCTTG" # Non-Island (-)
test_seq = "CAGGTTGGAAACGTAA"

bases = ['A', 'C', 'G', 'T']
base_to_idx = {base: i for i, base in enumerate(bases)}

def get_transition_matrix(seq, epsilon=1e-5):
    
    # Initialize 4x4 matrix with epsilon (pseudocount) to avoid log(0)
    matrix = np.full((4, 4), epsilon)

    # Count transitions
    for i in range(len(seq) - 1):
        row = base_to_idx[seq[i]]
        col = base_to_idx[seq[i+1]]
        matrix[row][col] += 1

    # Normalize rows to get probabilities
    row_sums = matrix.sum(axis=1, keepdims=True)
    prob_matrix = matrix / row_sums
    return prob_matrix

# 2. Generate Probabilities
prob_pos = get_transition_matrix(s1)
prob_neg = get_transition_matrix(s2)

# 3. Create Log-Likelihood Matrix (Base 2)
# Formula: log2(Prob_pos / Prob_neg)
ll_matrix = np.log2(prob_pos / prob_neg)

# Convert to DataFrame for readability
df_ll = pd.DataFrame(ll_matrix, index=bases, columns=bases)

print("Log-Likelihood Matrix:")
print("=" * 60)
print(df_ll.to_string())
print("\n")

# 4. Score the Test Sequence
total_score = 0
print(f"{'Transition':<12} | {'Log-Likelihood Score':<20}")
print("-" * 35)

for i in range(len(test_seq) - 1):
    b1, b2 = test_seq[i], test_seq[i+1]
    score = ll_matrix[base_to_idx[b1]][base_to_idx[b2]]
    total_score += score
    print(f"{b1} -> {b2:<8} | {score:>18.4f}")

print("-" * 35)
print(f"Final Total Score: {total_score:.4f}")

# 5. Result Decision
if total_score > 0:
    print("\nRESULT: The sequence is likely a CpG Island (+).")
else:
    print("\nRESULT: The sequence is likely NOT a CpG Island (-).")
