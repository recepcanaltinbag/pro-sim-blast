import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Amino asit sırası (HMMER sabit sırası)
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
               'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
               'T', 'V', 'W', 'Y']

def parse_hmm(file_path):
    with open(file_path) as f:
        lines = f.readlines()

    hmm_started = False
    data = []
    for line in lines:
        if line.startswith("HMM"):
            hmm_started = True
            continue
        if not hmm_started:
            continue
        if line.startswith("//"):
            break
        parts = line.strip().split()
        if len(parts) == 20:
            scores = [float(x) for x in parts]
            data.append(scores)

    df = pd.DataFrame(data, columns=amino_acids)
    return df

# Parse HMM file
df = parse_hmm("RieskeDB.hmm")

# Normalize scores (optional for visualization)
df_norm = df.apply(lambda x: x - x.max(), axis=1)  # convert log-odds to relative scale

# Heatmap
plt.figure(figsize=(14, 6))
sns.heatmap(df_norm.T, cmap="coolwarm", xticklabels=False)
plt.title("HMM Profile Heatmap (Log-Odds Scores)")
plt.ylabel("Amino Acids")
plt.xlabel("Position in HMM")
plt.tight_layout()
plt.show()

