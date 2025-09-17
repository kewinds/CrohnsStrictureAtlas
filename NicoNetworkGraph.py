import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Load the data
file_path = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Annotated/Stricture/Stricture_IF_Analysis/niche_prediction_linear/Niche_interactions_edge_weights_R30.csv"
df = pd.read_csv(file_path)

# Filter interactions involving CTHRC1+ mFib
filtered_df = df[(df['Source'] == 'CTHRC1+ mFib') | (df['Target'] == 'CTHRC1+ mFib')]

# Create directed graph
G = nx.DiGraph()

# Collect edges with weights from dataframe, ensuring bidirectional edges are captured
edges_with_weights = {}

for _, row in filtered_df.iterrows():
    source, target, weight = row['Source'], row['Target'], row['Weight']
    if (source, target) not in edges_with_weights:
        edges_with_weights[(source, target)] = {}
    edges_with_weights[(source, target)]['forward_weight'] = weight
    if (target, source) in edges_with_weights:
        edges_with_weights[(target, source)]['backward_weight'] = weight
    else:
        edges_with_weights[(target, source)] = {'backward_weight': weight}

# Handle cases where only one direction exists and ensure keys are correctly set
for edge in filtered_df.itertuples(index=False):
    u, v, w = edge.Source, edge.Target, edge.Weight
    if (v, u) not in edges_with_weights:
        edges_with_weights[(v, u)] = {'forward_weight': 0}

# Add bidirectional edges with weights properly annotated
for edge, weights in edges_with_weights.items():
    source, target = edge
    if "forward_weight" in weights:
        G.add_edge(source, target, weight=weights["forward_weight"])
    if "backward_weight" in weights:
        G.add_edge(target, source, weight=weights["backward_weight"])

# Set up plot
plt.figure(figsize=(14, 10))
pos = nx.spring_layout(G, seed=42, k=0.9)

# Create categorical color palette
nodes = list(G.nodes())
palette = plt.cm.tab20  # Using matplotlib's tab20 categorical palette
node_colors = [palette(i % palette.N) for i in range(len(nodes))]

# Draw nodes
nx.draw_networkx_nodes(
    G, pos, 
    node_size=2500, 
    node_color=node_colors,
    edgecolors='black',
    linewidths=2
)

# Draw edges with color differentiation and scaled thickness
edge_colors = []
edge_widths = []
arrow_sizes = []
for u, v in G.edges():
    weight = G[u][v]['weight']
    width = weight * 10  # Scale weight for visibility
    arrow_size = weight * 40  # Scale arrow size for visibility
    edge_widths.append(width)
    arrow_sizes.append(arrow_size)
    if u == 'CTHRC1+ mFib':  # Outgoing
        edge_colors.append('#0000FF')
    elif v == 'CTHRC1+ mFib':  # Incoming
        edge_colors.append('#FF0000')
    else:
        edge_colors.append('#333333')  # Default color for other edges

nx.draw_networkx_edges(
    G, pos,
    connectionstyle='arc3,rad=0.1',
    width=edge_widths,
    edge_color=edge_colors,
    arrowstyle='-|>',
    arrowsize=arrow_sizes,
    node_size=2500
)

# Node labels with contrasting text colors
nx.draw_networkx_labels(
    G, pos, 
    font_size=10,  # Decreased font size
    font_weight='bold',
    font_family='sans-serif',
    bbox=dict(
        facecolor='white',
        edgecolor='none',
        alpha=0.7
    )
)

# Ensure bidirectional edges show both weights with explicit labeling
edge_labels = {}
label_pos = {}
for u, v in G.edges():
    edge_labels[(u, v)] = f"{G[u][v]['weight']:.2f}"
    source_x, source_y = pos[u]
    target_x, target_y = pos[v]
    label_x = (0.7 * target_x + 0.3 * source_x)
    label_y = (0.7 * target_y + 0.3 * source_y)
    label_pos[(u, v)] = (label_x, label_y)
    # Label the reverse edge position if it exists
    if (v, u) in G.edges():
        rev_label_x = (0.3 * target_x + 0.7 * source_x)
        rev_label_y = (0.3 * target_y + 0.7 * source_y)
        label_pos[(v, u)] = (rev_label_x, rev_label_y)

# Draw edge labels for both directions to ensure non-overlapping
for (u, v), label in edge_labels.items():
    if (u, v) in label_pos:
        x, y = label_pos[(u, v)]
        # Set color based on direction
        label_color = '#0000FF' if u == 'CTHRC1+ mFib' else '#FF0000'
        plt.annotate(
            label,
            xy=(x, y),
            xytext=(0, 0),
            textcoords='offset points',
            ha='center',
            va='center',
            fontsize=10,  # Decreased font size
            color=label_color,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
        )
    if (v, u) in label_pos and (v, u) in edge_labels:
        rev_label = edge_labels[(v, u)]
        rev_x, rev_y = label_pos[(v, u)]
        # Set color for reverse direction
        label_color = '#0000FF' if v == 'CTHRC1+ mFib' else '#FF0000'
        plt.annotate(
            rev_label,
            xy=(rev_x, rev_y),
            xytext=(0, 0),
            textcoords='offset points',
            ha='center',
            va='center',
            fontsize=10,  # Decreased font size
            color=label_color,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
        )

# Final customization
plt.title("CTHRC1+ mFib Interaction Network (R=30)", 
         fontsize=18, pad=20, fontweight='bold')
plt.axis('off')

plt.tight_layout()

# Save the plot as a PDF
output_path = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Annotated/Stricture/Stricture_IF_Analysis/niche_prediction_linear/niche_interactions_R30.pdf"
plt.savefig(output_path, bbox_inches='tight', format='pdf')

plt.show()