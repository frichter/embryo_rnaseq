"""Process modules."""

import pandas as pd
import networkx as nx

mod_name = 'purple'  # salmon
mod_i = pd.read_csv('cs_edges_' + mod_name + '.txt', sep='\t')

mod_i = mod_i.loc[mod_i['weight'] > 0.05]

mod_i = mod_i.drop(['direction', 'fromAltName', 'toAltName'], axis=1)

out_loc = 'cs_edges_' + mod_name + '_min05.txt'
mod_i.to_csv(out_loc, index=False, sep="\t", float_format='%g')

node_l = mod_i.fromNode.tolist()
node_l.extend(mod_i.toNode.tolist())
len(set(node_l))
node_l = list(set(node_l))
with open('cs_edges_' + mod_name + '_min05_nodes.txt', 'w') as f:
    for node_i in node_l:
        _ = f.write(node_i + '\n')


mod_i = pd.read_csv('cs_edges_' + mod_name + '_min05_t5hubs.txt', sep='\t')
node_l = mod_i.fromNode.tolist()
node_l.extend(mod_i.toNode.tolist())
len(set(node_l))
node_l = list(set(node_l))
with open('cs_edges_' + mod_name + '_min05_nodes.txt', 'w') as f:
    for node_i in node_l:
        _ = f.write(node_i + '\n')

"""
import networkx as nx
import pandas as pd
"""

mod_name = 'purple'  # salmon
g = nx.Graph()

G = nx.Graph()
count = 0
net_loc = 'cs_edges_' + mod_name + '_min05.txt'
g_hubs = ["RPL7L1P8", "OR1F2P", "MINOS1P3", "AP006621.1", "RP11-585F1.10"][0]
with open(net_loc, 'r') as f:
    header = next(f)
    for line in f:
        line_wt = float(line.split("\t")[2])
        n1, n2 = line.split("\t")[0], line.split("\t")[1]
        # if (n1 in g_hubs) or (n2 in g_hubs):
        if line_wt >= 0.1:
            G.add_edge(n1, n2, weight=line_wt)


len(G.nodes())
len(G.edges())

with open('cs_edges_' + mod_name + '_min05_t1hub_nodes.txt', 'w') as f:
    for node_i in G.nodes():
        _ = f.write(node_i + '\n')


with open('cs_edges_' + mod_name + '_min1.txt', 'w') as f:
    _ = f.write('\t'.join(['fromNode', 'toNode', 'weight']))
    for u, v, d in G.edges(data=True):
        _ = f.write('\t'.join([u, v, str(d['weight'])]) + '\n')

g_hubs = ["RPL7L1P8", "OR1F2P", "MINOS1P3", "AP006621.1", "RP11-585F1.10"]
e_to_keep = [[u, v] for u, v in G.edges() if (u in g_hubs) or (v in g_hubs)]
# graph_list = list(nx.connected_component_subgraphs(G_sub))

with open('cs_edges_' + mod_name + '_min1_nodes.txt', 'w') as f:
    for node_i in G.nodes():
        _ = f.write(node_i + '\n')


#
#
#
#
#
