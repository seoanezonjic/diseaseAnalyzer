#!/usr/bin/env python

import argparse

#################################
##METHODS
#################################

def load_metrics_file(metrics_file):
    filtered_metrics = {}
    with open(metrics_file) as f:
        for line in f:
            line = line.rstrip()
            if 'group' in line:
                continue
            cluster_id, asp_value = line.split('\t')
            filtered_metrics[cluster_id] = float(asp_value)
    return filtered_metrics

def load_clusters_file(clusters_file):
    clustered_genes = {}
    with open(clusters_file) as f:
        for line in f:
            line = line.rstrip()
            cluster_id, gene = line.split('\t')
            if cluster_id not in clustered_genes:
                clustered_genes[cluster_id] = [gene]
            else:
                clustered_genes[cluster_id].append(gene)
    return clustered_genes

def combine_data(filtered_metrics, clustered_genes, contributing_diseases_per_cluster, output_file):
    storage = []
    for cluster, genes in clustered_genes.items():
        disease_per_cluster = contributing_diseases_per_cluster.get(cluster, 0)
        asp_value = filtered_metrics.get(cluster, None)
        if asp_value is not None:
            storage.append([len(genes), asp_value, disease_per_cluster, cluster])
    with open(output_file, 'w') as f:
        f.write('clusterSize\taspValue\tdiseasePerCluster\tclusterID\n')
        for genes, asp_value, disease_per_cluster, cluster in storage:
            f.write(f'{genes}\t{asp_value}\t{disease_per_cluster}\t{cluster}\n')

def load_diseases_file(cluster_diseases_file, disease_gene_file, clustered_genes_storage):
    cluster_diseases = {}
    with open(cluster_diseases_file) as f:
        for line in f:
            line = line.rstrip()
            disease_id, cluster_id = line.split('\t')
            if cluster_id not in cluster_diseases:
                cluster_diseases[cluster_id] = [disease_id]
            else:
                cluster_diseases[cluster_id].append(disease_id)

    disease_genes = {}
    with open(disease_gene_file) as f:
        for line in f:
            line = line.rstrip()
            disease_id, gene_id = line.split('\t')
            if disease_id not in disease_genes:
                disease_genes[disease_id] = [gene_id]
            else:
                disease_genes[disease_id].append(gene_id)

    contributing_diseases_per_cluster = {}
    for cluster_id, disease_ids in cluster_diseases.items():
        count = 0
        for disease in disease_ids:
            gene_ids = disease_genes.get(disease, [])
            if gene_ids:
                count += 1
        contributing_diseases_per_cluster[cluster_id] = count / len(disease_ids) * 100

    return contributing_diseases_per_cluster

#################################
##ARG-PARSER
#################################

parser = argparse.ArgumentParser(description=f'Usage: {__file__} [options]')
parser.add_argument('-d', '--diseases_file', required=True, help='Diseases file input')
parser.add_argument('-i', '--metrics_file', required=True, help='Metric file input')
parser.add_argument('-f', '--clusters_file', required=True, help='Clusters file input')
parser.add_argument('-g', '--disease_gene_file', required=True, help='Disease-gene file input')
parser.add_argument('-o', '--output_file', default='output.txt', help='Output combined file')
args = parser.parse_args()

#################################
##MAIN
#################################

filtered_metrics = load_metrics_file(args.metrics_file)
clustered_genes = load_clusters_file(args.clusters_file)
contributing_diseases_per_cluster = load_diseases_file(args.diseases_file, args.disease_gene_file, clustered_genes)
combine_data(filtered_metrics, clustered_genes, contributing_diseases_per_cluster, args.output_file)