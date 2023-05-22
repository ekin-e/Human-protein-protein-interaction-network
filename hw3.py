import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Read in the data from the CSV file
df = pd.read_csv('HURI.hgnc.csv')

# p_col1 = df[0].str.split(',', expand=True)
# p_col2 = df[1].str.split(',', expand=True)

# new_data = []
# for i in range(len(p_col1)):
#     if p_col1.iloc[i][0] != p_col2.iloc[i][0]:
#         new_data.append([p_col1[i][0], p_col2.iloc[i][0]])

# Filter out self-interactions
edges = [(row[0], row[1]) for row in df.values if row[0] != row[1]]

# Extract the edges as a list of tuples
# edges = [tuple(x) for x in df.values]

# Construct a graph object from the list of edges
G = nx.Graph()
G.add_edges_from(edges)

# Calculate the clustering coefficient for each node (excluding self-loops)
clustering_coeffs = nx.clustering(G)

# Calculate the degree of each protein
degrees = dict(G.degree())

# Plot the degree distribution as a histogram
plt.hist(list(degrees.values()), bins=50)
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.show()

# Q2
highest_degree_node = max(degrees, key=degrees.get)
highest_degree = degrees[highest_degree_node]
print(f"The highest degree protein is {highest_degree_node} with a degree of {highest_degree}.")

# Q3
# Create a list of (degree, clustering coefficient) tuples for each node
data = [(degrees[node], clustering_coeffs[node]) for node in degrees.keys()]
# Plot the data as a scatter plot
x, y = zip(*data)
plt.scatter(x, y, alpha=0.5)
plt.xlabel('Degree')
plt.ylabel('Clustering Coefficient')
plt.show()

#Q3c
# Find the degree and clustering coefficient for LSM6 and MAPK9, or print an error message if they are not found
try:
    lsm6_degree = G.degree['LSM6']
    lsm6_clustering = clustering_coeffs['LSM6']
    lsm6_neighbors = list(G.neighbors('LSM6'))
except KeyError:
    print("LSM6 not found in the graph")
    lsm6_degree = None
    lsm6_clustering = None
    lsm6_neighbors = []

try:
    mapk9_degree = G.degree['MAPK9']
    mapk9_clustering = clustering_coeffs['MAPK9']
    mapk9_neighbors = list(G.neighbors('MAPK9'))
except KeyError:
    print("MAPK9 not found in the graph")
    mapk9_degree = None
    mapk9_clustering = None
    mapk9_neighbors = []

# Print the results
if lsm6_degree is not None:
    print("LSM6:")
    print("Number of interaction partners:", lsm6_degree)
    print("Clustering coefficient:", lsm6_clustering)
    print("Interacting proteins:", lsm6_neighbors)
    print()
if mapk9_degree is not None:
    print("MAPK9:")
    print("Number of interaction partners:", mapk9_degree)
    print("Clustering coefficient:", mapk9_clustering)
    print("Interacting proteins:", mapk9_neighbors)

# Q4
lit_df = pd.read_csv('Lit_degrees.csv')
lit_degrees = dict(zip(lit_df['Protein\tDegree'].str.split('\t').str[0], lit_df['Protein\tDegree'].str.split('\t').str[1].astype(int)))

# Get the set of proteins with at least one interaction in both networks
common_proteins = set(G.nodes()).intersection(set(lit_degrees.keys()))

# Calculate the degrees of each protein in the two networks
sys_map_degrees = {protein: degree for protein, degree in G.degree(common_proteins)}
lit_degrees = {protein: degree for protein, degree in lit_degrees.items() if protein in common_proteins}


# Calculate the Pearson correlation between the two sets of degrees
corr, _ = pearsonr(list(sys_map_degrees.values()), list(lit_degrees.values()))

print(f"Pearson correlation between systematically mapped and literature-curated network degrees: {corr:.3f}")

#Q4b
protein_of_interest = None
for protein in G.nodes():
    if G.degree(protein) > 10 and nx.clustering(G, protein) > 0.2 and protein not in lit_degrees.keys():
        protein_of_interest = protein
        break
        
print(f"Protein with more than 10 interactions, clustering coefficient greater than 0.2, and no interactions in the literature-curated network: {protein_of_interest}")

# Interaction partners of the protein CELF5
try:
    celf5_neighbors = list(G.neighbors('CELF5'))
    print("CELF5 interaction partners:", celf5_neighbors)
except KeyError:
    print("CELF5 not found in the graph")



