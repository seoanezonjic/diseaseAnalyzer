#! /usr/bin/env python


import argparse
import re

##############
## METHODS
##############


def load_dictionary(file):
	dictionary = {}
	with open(file, 'r') as f:
		for line in f:
			line = line.strip()
			if re.match(r'^#', line): continue
			string, genename, protein_size, annotation = line.split("\t")
			query = dictionary.get(string)
			if query is None:
				dictionary[string] = genename
			else:
				query.append(genename)
	return dictionary

def translate_file(file, dictionary):
	translations = []
	untranslated_prots = []
	with open(file, 'r') as f:
		for line in f:
			line = line.strip() # equivalent to chomp! in rb
			if re.search('protein1', line) is None:
				info = line.split(' ')
				comb_score = info[-1]
				geneID1 = dictionary[info[0]]
				geneID2 = dictionary[info[1]]
				if geneID1 is None:
					untranslated_prots.append(info[0])
				elif geneID2 is None:
					untranslated_prots.append(info[1])
				else:
					translations.append([geneID1, geneID2, comb_score])
			else:
				header = line.split("\t")

	return translations, list(set(untranslated_prots))


##############
## OPTPARSE
##############
parser = argparse.ArgumentParser(description="Translate cluster STRING Ids to gene clusters")
parser.add_argument("-d", "--dictionary_file", metavar="PATH", help="Input file with dictionary of terms", required=True)
parser.add_argument("-i", "--input_file", metavar="PATH", help="Input file to translate", required=True)
parser.add_argument("-n", "--untranslated_file", metavar="PATH", help="Output file with untranslated proteins", required=True)
parser.add_argument("-o", "--output_file", metavar="PATH", help="Output file", default="output.txt")
args = parser.parse_args()

##############
## MAIN
##############

dictionary = load_dictionary(args.dictionary_file)
translated_info, untranslated_prots = translate_file(args.input_file, dictionary)

with open(args.output_file, "w") as f:
	f.write("protein1\tprotein2\tcombined_score\n")
	for info in translated_info:
		f.write("\t".join(info) + "\n")

with open(args.untranslated_file, "w") as f:
	for protID in untranslated_prots:
		f.write(protID + "\n")