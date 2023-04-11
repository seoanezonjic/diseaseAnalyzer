#! /usr/bin/env python

import argparse
import os
import re
import json

##############
## METHODS
##############

def load_pairs(input_file, cols, types):
    pairs = []
    count = 0
    node_ids = {}
    nodes = []
    with open(input_file) as file:
        for line in file:
            fields = line.strip().split('\t')
            node_cols = [fields[c] for c in cols]
            node_cols = [nc.split(',') for nc in node_cols]
            for i, nodes_A in enumerate(node_cols):
                nodes_B = node_cols[i+1] if i < len(node_cols)-1 else None
                for node_A in nodes_A:
                    if node_A == '-':
                        continue
                    query_node = node_ids.get(node_A)
                    if query_node is None:
                        node_ids[node_A] = str(count)
                        add_node(nodes, node_A, str(count), types[i])
                        count += 1
                    if nodes_B is not None:
                        for node_B in nodes_B:
                            if node_B == '-':
                                continue
                            pairs.append([node_A, node_B])
    return pairs, node_ids, nodes, count

def add_node(nodes, name, count, node_type):
    node = {
        'data': {
            'id': count,
            'name': name,
        }
    }
    if node_type is not None:
        node['data']['type'] = node_type
    nodes.append(node)

def pairs2cys(pairs, node_ids, nodes, count):
    edges = []
    for node_A, node_B in pairs:
        edges.append({
            'data': {
                'id': str(count),
                'source': node_ids[node_A],
                'target': node_ids[node_B],
                'interaction': '-',
                'weight': 1.0
            }
        })
        count += 1
    cys = {
        'elements': {
            'nodes': nodes,
            'edges': edges
        }
    }
    return cys


##############
## OPTPARSE
##############

options = {}
parser = argparse.ArgumentParser(description=f'Usage: {__file__} [options]')

parser.add_argument('-i', '--input', dest='input', help='Path to input file')
parser.add_argument('-o', '--output', dest='output', help='Path to output file')
parser.add_argument('-c', '--cols', dest='cols', default='0,1', help='Columns with network nodes. 0 based. Default, 0,1')
parser.add_argument('-t', '--types', dest='types', default='', help='node types, one for each selected col. 0 based. Default, 0,1')

args = parser.parse_args()

options['input'] = args.input
options['output'] = args.output
options['cols'] = [int(col) for col in args.cols.split(',')]
options['types'] = args.types.split(',')


##############
## MAIN
##############

pairs, node_ids, nodes, count = load_pairs(args.input, args.cols, args.types)
q()
# Convert pairs to Cytoscape.js format
cys_net = pairs2cys(pairs, node_ids, nodes, count)

# Write Cytoscape.js JSON to file
with open(args.output + '.cyjs', 'w') as f:
    f.write(json.dumps(cys_net, indent=4))