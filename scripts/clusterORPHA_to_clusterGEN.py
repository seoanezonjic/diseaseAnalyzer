#! /usr/bin/env python

import argparse
import os
import re
from collections import defaultdict

##############
## METHODS
##############

def load_disease_clusters(file):
	disease_clusters = {}
	with open(file, 'r') as f:
		for line in f:
			line = line.strip()
			orpha_code, cluster_id = line.split("\t")
			disease_clusters[orpha_code] = cluster_id
	return disease_clusters

def get_cluster_genes(orpha_genes_file, disease_clusters):
	cluster_genes = {}
	with open(orpha_genes_file, 'r') as f:
			for line in f:
				line = line.strip()
				if 'DiseaseID' in line: continue
				orpha_code, genes = line.split("\t")
				genes_list = genes.split(',')
				cluster_id = disease_clusters.get(orpha_code)
				if cluster_id is not None:
					query = cluster_genes.get(cluster_id)
					if query is None:
						cluster_genes[cluster_id] = [genes_list]		
					else:
						query.append(genes_list)
	cluster_genes = {cluster: genes_by_disease for cluster, genes_by_disease in cluster_genes.items() if len(genes_by_disease) >= 2}			
	return cluster_genes

def apply_filters(cluster_genes, filter_value, min_groups):
	saved_genes = {}
	general_stats = defaultdict(int)
	for cluster, genes_by_disease in cluster_genes.items():
		stats = defaultdict(int)
		for gene_group in genes_by_disease:
			for gene in gene_group:
				stats[gene] += 1
				general_stats[gene] += 1
		filtered_genes = []
		for gene, number in stats.items():
			if number < min_groups: continue
			if number / len(genes_by_disease) >= filter_value:
				filtered_genes.append(gene)
		saved_genes[cluster] = filtered_genes
	return saved_genes, general_stats

##############
## OPTPARSE
##############

parser = argparse.ArgumentParser(description="Get disease clusters (cluster IDs and genes by disease")
parser.add_argument("-F", "--filter_value", type=float, help="Filter value", required=True, default=0.8)
parser.add_argument("-i", "--input_file", metavar="PATH", help="Input file", required=True)
parser.add_argument("-m", "--min_groups", type=int, help="Min diseases per gene", required=True, default=0)
parser.add_argument("-o", "--output_file", metavar="PATH", help="Output file", default="output.txt")
parser.add_argument("-t", "--orpha_genes", metavar="PATH", help="ORPHA-genes file (genes separated by commas)")
args = parser.parse_args()

##############
## MAIN
##############

disease_clusters = load_disease_clusters(args.input_file)
cluster_genes = get_cluster_genes(args.orpha_genes, disease_clusters)
saved_genes, gene_stats = apply_filters(cluster_genes, args.filter_value, args.min_groups)

stats_file = os.path.join(os.path.dirname(args.output_file), "gene_stats.txt")
with open(stats_file, "w") as f:
	f.write("gene\tfreq\n")
	for gene, number in gene_stats.items():
		f.write(gene + "\t" + str(number) + "\n")

with open(args.output_file, "w") as f:
	final_genes = []
	for clusterID, genes in saved_genes.items():
		if len(genes) > 0:
			for gene in genes:
				f.write(clusterID + "\t" + gene + "\n")