#! /usr/bin/env python

import argparse
import re

##############
## METHODS
##############


def load_mondo_file(file):
    #mondo_genes = {MONDO: [geneA, geneB]}
    mondo_genes = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.rstrip()
            mondo_id, gene_hgnc, gene_symb = line.split("\t")
            query = mondo_genes.get(mondo_id)
            if query is None:
                mondo_genes[mondo_id] = [gene_symb]
            else:
                if gene_symb not in query: query.append(gene_symb)
    return mondo_genes

def find_genes(input_file, mondo_genes, output_file):
    orpha_mondo_genes = []
    with open(input_file) as f:
        for line in f:
            line = line.rstrip()
            orpha_id, mondo_id = line.split("\t")
            disease_genes = mondo_genes.get(mondo_id)
            if disease_genes is not None:
                for gene in disease_genes:
                    orpha_mondo_genes.append([orpha_id, mondo_id, gene])
    with open(output_file, 'w') as f:
        for info in orpha_mondo_genes:
            f.write('\t'.join(info) + '\n')

##############
## OPTPARSE
##############
parser = argparse.ArgumentParser(description="Find genes associated with diseases using ORPHA and MONDO codes.")
parser.add_argument("-i", "--input_file", metavar="PATH", help="List ORPHA and MONDO codes", required=True)
parser.add_argument("-m", "--mondo_file", metavar="PATH", help="MONDO file to find genes associated with disease", required=True)
parser.add_argument("-o", "--output_file", metavar="PATH", help="Output ORPHA, MONDO and genes file", default="output.txt")
args = parser.parse_args()

##############
## MAIN
##############

mondo_genes = load_mondo_file(args.mondo_file)
find_genes(args.input_file, mondo_genes, args.output_file)        