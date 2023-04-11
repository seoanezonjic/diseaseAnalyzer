#! /usr/bin/env python

import argparse
import csv

#################################
##METHODS
#################################

def load_orpha_mondo_genes_file(input_file):
    orpha_genes_storage = {}
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            orpha_code, mondo_code, gene_code = line.split("\t")
            query = orpha_genes_storage.get(orpha_code)
            if query is None:
                orpha_genes_storage[orpha_code] = [gene_code]
            else:
                query.append(gene_code)
    return orpha_genes_storage

def load_orpha_hpo_file(input_file):
    orpha_hpos_storage = {}
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if 'DiseaseID' in line:
                continue
            orpha_code, hpos = line.split("\t")
            orpha_hpos = hpos.split('|')
            orpha_hpos_storage[orpha_code] = orpha_hpos
    return orpha_hpos_storage

def load_orpha_clusters_file(input_file):
    orpha_clusters_storage = {}
    clusters_orpha_storage = {}
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            orpha_code, cluster_id = line.split("\t")
            orpha_clusters_storage[orpha_code] = cluster_id
            query = clusters_orpha_storage.get(cluster_id)
            if query is None:
                clusters_orpha_storage[cluster_id] = [orpha_code]
            else:
                query.append(orpha_code)
    return orpha_clusters_storage, clusters_orpha_storage

def combine_files(orpha_genes_storage, orpha_hpos_storage, orpha_clusters_storage, clusters_orpha_transl_storage):
    final_table = {}
    orpha_codes = set(orpha_genes_storage.keys()) | set(orpha_hpos_storage.keys())
    for orpha_code in orpha_codes:
        genes_ary = orpha_genes_storage.get(orpha_code)
        hpos_ary = orpha_hpos_storage.get(orpha_code)
        cluster_id = orpha_clusters_storage.get(orpha_code)
        tranls_dis = clusters_orpha_transl_storage.get(orpha_code)
        if hpos_ary is not None and cluster_id is not None:
            genes_ary = ['-'] if genes_ary is None else list(set(genes_ary))
            info_to_add = [tranls_dis, cluster_id, ",".join(genes_ary), ",".join(hpos_ary)]
            final_table[orpha_code] = info_to_add
    return final_table

def translate_orpha_diseases(orpha_dictionary_file, clusters_orpha_storage):
    clusters_orpha_transl_storage = {}
    with open(orpha_dictionary_file) as f:
        for line in f:
            disease_code, disease_name, hpo_code = line.strip().split("\t")
            clusters_orpha_transl_storage[disease_code] = disease_name
    translated_diseases = {}
    for cluster, diseases in clusters_orpha_storage.items():
        for disease_id in diseases:
            disease_name = clusters_orpha_transl_storage.get(disease_id)
            query = translated_diseases.get(cluster)
            if query is None:
                translated_diseases[cluster] = [disease_name]
            else:
                query.append(disease_name)
    return translated_diseases, clusters_orpha_transl_storage

def generate_output(output_file, string, data2write):
    with open(output_file, 'w') as f:
        f.write(string)
        for key, val in data2write.items():
            f.write(f"{key}\t{','.join(val)}\n")

#################################
##ARG-PARSER
#################################

parser = argparse.ArgumentParser(description='Parse info to create a table')

parser.add_argument('-c', '--orpha_clusters_file', type=str,
                    help='Input file with ORPHA codes & cluster IDs')
parser.add_argument('-d', '--orpha_dictionary', type=str,
                    help='Input file with ORPHA codes & description. Used for translation.')
parser.add_argument('-i', '--orpha_mondo_genes_file', type=str,
                    help='Input file with ORPHA, MONDO and gene codes')
parser.add_argument('-f', '--orpha_hpo_file', type=str,
                    help='Input file with ORPHA and HPO codes')
parser.add_argument('-o', '--output_file', type=str, default='output_table.txt',
                    help='Output file')
parser.add_argument('-O', '--output_clusters_orpha_table', type=str, default='output_clusters_orpha_table.txt',
                    help='Output file: table with clusters and diseases by cluster.')

args = parser.parse_args()

options = {
    'orpha_clusters_file': args.orpha_clusters_file,
    'orpha_dictionary': args.orpha_dictionary,
    'orpha_mondo_genes_file': args.orpha_mondo_genes_file,
    'orpha_hpo_file': args.orpha_hpo_file,
    'output_file': args.output_file,
    'output_clusters_orpha_table': args.output_clusters_orpha_table
}

#################################
##MAIN
#################################

orpha_genes_storage = load_orpha_mondo_genes_file(args.orpha_mondo_genes_file)
orpha_hpos_storage = load_orpha_hpo_file(args.orpha_hpo_file)
orpha_clusters_storage, clusters_orpha_storage = load_orpha_clusters_file(args.orpha_clusters_file)
translated_diseases, clusters_orpha_transl_storage = translate_orpha_diseases(args.orpha_dictionary, clusters_orpha_storage)

#diseases_num_by_cluster = list(map(lambda a: len(a), clusters_orpha_storage.values()))
#ave_diseases_by_cluster = sum(diseases_num_by_cluster) / len(diseases_num_by_cluster)

final_table = combine_files(orpha_genes_storage, orpha_hpos_storage, orpha_clusters_storage, clusters_orpha_transl_storage)

with open(options['output_file'], 'w') as f:
	f.write("DiseaseID\tDiseaseName\tClusterID\tGeneIDs\tHPOIDs\n")
	for orpha_code, info in final_table.items():
		f.write(orpha_code + "\t" + "\t".join(info) + "\n")

with open(options['output_clusters_orpha_table'], 'w') as f:
	f.write("ClusterID\tDiseaseIDs\n")
	for clusterID, diseases in translated_diseases.items():
		f.write(clusterID + "\t" + ",".join(diseases) + "\n")