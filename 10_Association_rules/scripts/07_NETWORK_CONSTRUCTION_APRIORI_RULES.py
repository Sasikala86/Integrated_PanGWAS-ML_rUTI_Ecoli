# ==========================================
# PPI-LIKE NETWORK FROM APRIORI RULES
# ==========================================

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import ast
import re
from matplotlib import cm
from matplotlib.colors import Normalize

# FILE PATHS
rules_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/rUTI_specific/FINAL_biologically_significant_rUTI_rules.csv"
output_file = "/home/sastra/Sasikala/Pangenome/ML_Accessory/results/rUTI_specific/rUTI_gene_network.png"

# LOAD DATA
df = pd.read_csv(rules_file)
df["antecedents"] = df["antecedents"].apply(ast.literal_eval)

# CLEAN GENES (REMOVE Panaroo clusters + _1,_2)
def clean_gene(g):
    if g.startswith("group_"):
        return g
    return re.sub(r'_\d+$', '', g)

def expand_and_clean(lst):
    clean_genes = []
    for g in lst:
        parts = g.split("~~~")
        for p in parts:
            clean_genes.append(clean_gene(p.strip()))
    return list(set(clean_genes))

df["genes"] = df["antecedents"].apply(expand_and_clean)

# BUILD GRAPH
G = nx.Graph()

for _, row in df.iterrows():
    genes = row["genes"]
    lift = row["lift"]
    conf = row["confidence"]

    for i in range(len(genes)):
        for j in range(i+1, len(genes)):
            g1, g2 = genes[i], genes[j]

            if G.has_edge(g1, g2):
                G[g1][g2]["weight"] = max(G[g1][g2]["weight"], lift)
                G[g1][g2]["confidence"] = max(G[g1][g2]["confidence"], conf)
            else:
                G.add_edge(g1, g2, weight=lift, confidence=conf)

# DRAW NETWORK
fig, ax = plt.subplots(figsize=(12, 10))

pos = nx.spring_layout(G, k=0.5, seed=42)

edges = G.edges(data=True)
weights = [d["weight"] for (_, _, d) in edges]
confidences = [d["confidence"] for (_, _, d) in edges]

norm = Normalize(vmin=min(confidences), vmax=max(confidences))
cmap = cm.viridis

edge_colors = [cmap(norm(c)) for c in confidences]

nx.draw_networkx_edges(
    G, pos,
    width=[w * 2 for w in weights],
    edge_color=edge_colors,
    ax=ax
)

nx.draw_networkx_nodes(
    G, pos,
    node_size=800,
    node_color="lightblue",
    ax=ax
)

nx.draw_networkx_labels(
    G, pos,
    font_size=8,
    ax=ax
)

# ✅ ONLY ONE COLORBAR (CORRECT)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax)
cbar.set_label("Confidence")

# FINAL OUTPUT
ax.set_title("rUTI Gene Association Network\nEdge thickness = Lift | Edge color = Confidence")
ax.axis('off')

plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.show()

print("Network saved to:", output_file)
