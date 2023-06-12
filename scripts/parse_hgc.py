#! /usr/bin/env python

import argparse
import os
import re

##############
## METHODS
##############

def load_tabular(file):
    data = []
    with open(file) as f:
        for line in f:
            data.append(line.strip().split('\t'))
    return data

def load_clusters(file, nodes):
    clusters = {}
    for cl_id, node in load_tabular(file):
        if node in nodes:
            add_hash(cl_id, node, clusters)
        else:
            print(f"{node} node id from cluster {cl_id} does not exist in network")
    return clusters

def load_paths(paths_file, thr=None):
    path_data = {}
    nodes = {}
    for fields in load_tabular(paths_file):
        source, target = fields[:2]
        nodes[source] = True
        nodes[target] = True
        if thr is not None:
            sphere = int(fields[5])
            if sphere <= thr:
                fields[5] = sphere
            else:
                continue
        fields[0] = float(fields[0])  # Biological distance
        fields[6] = [g.split('[')[0] for g in fields[6].split(']')]  # Route
        add_hash(source, {target: fields}, path_data)
    return path_data, nodes

def load_paths_simple(paths_file):
    path_data = {}
    nodes = {}
    for fields in load_tabular(paths_file):
        source, target = fields[:2]
        nodes[source] = True
        nodes[target] = True
        fields[0] = float(fields[0]) #Biological distance
        fields[6] = fields[1].split(',') #Route
        add_hash(source, {target: fields}, path_data)
    return path_data, nodes

def write_tabular(data, file):
    with open(file, 'w') as f:
        for record in data:
            f.write('\t'.join(map(str, record)) + '\n')

def write_nested_tabular(data, file, col):
    with open(file, 'w') as f:
        for record in data:
            items = record[col]
            for item in items:
                record[col] = item
                f.write('\t'.join(map(str, record)) + '\n')

def add_hash(key, val, hash):
    query = hash.get(key)
    if query is None:
        if isinstance(val, dict):
            hash[key] = val
        else:
            hash[key] = [val]
    else:
        if isinstance(val, dict):
            query.update(val)
        else:
            query.append(val)

def compact_hgc(input_folder, gene_list=None):
    all_data = []
    for file_path in os.listdir(input_folder):
        if not file_path.endswith('.txt'):
            continue
        gene = os.path.splitext(file_path)[0]
        if gene_list is not None and gene not in gene_list:
            continue
        gene_data = load_tabular(os.path.join(input_folder, file_path))
        gene_data.pop(0)  # Remove header
        gene_data.pop(0)  # Remove gene self path
        #if gene_list is not None:
        #    gene_data = [g for g in gene_data if g[0] in gene_list]
        for g in gene_data:
            g.insert(0, gene)
        all_data.extend(gene_data)
    return all_data


def expand_clusters(clusters, paths, partial=False):
    expanded_clusters = []
    all_stats = []
    for id, genes in clusters.items():
        expanded_genes, avg_sh_path = expand_cluster(genes, paths, partial)
        expanded_clusters.append([id, expanded_genes])
        all_stats.append([id, avg_sh_path])
    return expanded_clusters, all_stats

def expand_cluster(genes, paths, partial=False):
	#Takes in a set of genes and a set of paths as inputs, and expands the set of genes by adding intermediate genes along the paths between each pair of genes in the original set
    genes = genes.copy()
    expanded_genes = []
    connected = True
    path_lengths = []
    while len(genes) > 1 and connected:
        source = genes.pop()
        for target in genes:
            pair_data = get_path_data(source, target, paths)
            if pair_data is not None:
                path = pair_data[6]
                path_lengths.append(len(path))
                expanded_genes.extend(path)
            else:
                print(f"{source} - {target} has no path", file=sys.stderr)
                path_lengths.append(None)
                if partial:
                    expanded_genes.extend([source, target])
                else:  # No path between source and target
                    connected = False
                    break
    expanded_genes = list(set(expanded_genes))
    n_genes = len(expanded_genes)
    avg_p_len = float('nan') if None in path_lengths or n_genes == 0 else sum(path_lengths) / n_genes
    return expanded_genes, avg_p_len

def get_network_distances(seed_groups, group_candidates, paths):
    distances_by_group = {}
    for id, candidates in group_candidates.items():
        seeds = seed_groups[id]
        candidates = candidates - seeds
        distances = {}
        for candidate in candidates:
            distance = get_average_path_distance(candidate, seeds, paths)
            distances[candidate] = distance
        if distances:
            distances_by_group[id] = distances
    return distances_by_group

def get_average_path_distance(candidate, seeds, paths):
    distances = []
    for target in seeds:
        path_data = get_path_data(candidate, target, paths)
        if path_data is not None:
            distances.append(path_data[0])
    if distances:
        return sum(distances) / len(distances)
    else:
        return float('nan')

def get_path_data(source, target, paths):
    pair_data = paths.get(source, {}).get(target)
    if pair_data is None:
        pair_data = paths.get(target, {}).get(source)
    return pair_data

##############
## OPTPARSE
##############

parser = argparse.ArgumentParser(description='Usage: ' + __file__ + ' [options]')

parser.add_argument('-i', '--input_file', dest='input_file', help='HGC aggregated data')
parser.add_argument('-f', '--input_folder', dest='input_folder', help='Folder with HGC specific connectome files files')
parser.add_argument('-g', '--gene_list', dest='gene_list', help='File with gene list to select pairs')
parser.add_argument('-p', '--partial', dest='partial', action='store_true', help='Allows clusters with incomplete paths')
parser.add_argument('-S', '--simple', dest='simple', action='store_true', help='Simple data path')
parser.add_argument('-o', '--output_file', dest='output_file', default='output.txt', help='Output combined file')
parser.add_argument('-s', '--stat_file', dest='stat_file', help='Output file with clusters statistics')
parser.add_argument('-t', '--threshold', dest='threshold', type=float, help='p-value threshold')

options = vars(parser.parse_args())

##############
## MAIN
##############

if "input_folder" in options:
    gene_list = None
    if "gene_list" in options and options["gene_list"]:
        gene_list = load_tabular(options["gene_list"])
        gene_list = [item for sublist in gene_list for item in sublist]
    all_data = compact_hgc(options["input_folder"], gene_list)
    write_tabular(all_data, options["output_file"])
elif "input_file" in options:
    if options["simple"]:
        paths, nodes = load_paths_simple(options["input_file"])
    else:
        paths, nodes = load_paths(options["input_file"], options["threshold"])
    clusters = load_clusters(options["gene_list"], nodes)
    expanded_clusters, all_stats = expand_clusters(clusters, paths, options["partial"])
    distances_by_cluster = get_network_distances(clusters, expanded_clusters, paths)
    with open(os.path.join(os.path.dirname(options["output_file"]), "distances.txt"), "w") as f:
        for cl_id, gene_distances in distances_by_cluster.items():
            gene_distances = sorted(gene_distances.items(), key=lambda x: x[1])
            count = 1
            for gene, distance in gene_distances:
                if not math.isnan(distance):
                    f.write("\t".join([cl_id, gene, str(distance), str(count / len(gene_distances))]) + "\n")
                    count += 1
    write_nested_tabular(expanded_clusters, options["output_file"], 1)
    if options["stat_file"]:
        write_tabular(all_stats, options["stat_file"])
