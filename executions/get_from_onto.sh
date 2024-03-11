#! /usr/bin/env bash

current=`pwd`
global_path=$current/..
scripts_path=$global_path/scripts
dataset_path=$global_path'/datasets'
list_path=$global_path'/lists'
temp_files=$global_path'/executions/temp_files'
ontology_file=$dataset_path'/mondo.obo'
export PATH=$scripts_path:$PATH
source ~soft_bio_267/initializes/init_python


if [ "$1" == "1" ]; then
	mkdir $list_path
	term=$2
	if [ "$3" == "orpha" ]; then
		codename='Orphanet'
		codetag='ORPHA'
	elif [[ "$3" == "omim" ]]; then
		codename='OMIM'
		codetag='OMIM'
	fi
	monarch_entities -d $term -r phenotype-disease -o $list_path/$term'_disease_results.tab'
	cut -f 1 $list_path/$term'_disease_results.tab' | sort -u > $list_path/$term.mondo
	semtools -i $list_path/$term.mondo --list -k $codename":[0-9]*" --xref_sense -O MONDO -o $list_path/$term'.dict'
	cut -f 2 $list_path/$term'.dict' | sort -u | sed "s/$codename/$codetag/g" > $dataset_path/$term.txt
fi

if [ "$1" == "2" ]; then

	wget https://data.monarchinitiative.org/latest/tsv/gene_associations/gene_disease.9606.tsv.gz -O $temp_files'/gene_disease.9606.tsv.gz'
	gunzip $temp_files'/gene_disease.9606.tsv.gz'
	# the preivous source seems deprecated. See that this source could replace iti (from https://hpo.jax.org/app/data/annotations) or if there is an alternative in the new monarch app:
	# wget 'https://objects.githubusercontent.com/github-production-release-asset-2e65be/41063438/545054dd-747d-4bba-ab19-4f18858f39af?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAVCODYLSA53PQK4ZA%2F20240311%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240311T180105Z&X-Amz-Expires=300&X-Amz-Signature=08e1f83f4ccb9b4f91d07302eec125a1ac646e596a66cef106cba8fd3b0f4edd&X-Amz-SignedHeaders=host&actor_id=4126037&key_id=0&repo_id=41063438&response-content-disposition=attachment%3B%20filename%3Dgenes_to_disease.txt&response-content-type=application%2Foctet-stream' -O genes_to_disease.txt
	# maybe we have to use the new api to retrieve this data (see module 1)


	wget http://geneontology.org/gene-associations/goa_human.gaf.gz -O $temp_files'/goa_human.gaf.gz'
	gunzip $temp_files'/goa_human.gaf.gz'

fi

if [ "$1" == "3" ]; then
	go=$2 #GO:0001525 => Angiogenesis
	source ~soft_bio_267/initializes/init_ruby
	cut -f 2,5 $temp_files'/gene_disease.9606.tsv' > $temp_files'/gene_MONDO4go.txt'
	grep -w $go $temp_files'/goa_human.gaf' | awk '{if($7 != "IEA" && $7 != "NAS" && $7 != "ISS") print $3}' | sort -u > $list_path/go_gene.lst
	grep -w -F -f $list_path/go_gene.lst $temp_files'/gene_MONDO4go.txt' > $list_path'/gene_MONDO4go_filt.txt'
	cut -f 1 $list_path'/gene_MONDO4go_filt.txt' | sort -u > $list_path'/go_genes.txt'
	cut -f 2 $list_path'/gene_MONDO4go_filt.txt' | sort -u > $list_path'/go_mondo.txt'
	get_disease_mondo.rb -i $list_path'/go_mondo.txt' -k 'Orphanet:[0-9]*' -s -f $ontology_file -S $temp_files/supp_mondo_orpha.txt -o $list_path/go_query.dict
	cut -f 2 $list_path/go_query.dict | sort -u | sed 's/Orphanet/ORPHA/g' > $dataset_path/$go.txt
fi
