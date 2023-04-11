#! /usr/bin/env bash

current=`pwd`
global_path=$current/..
scripts_path=$global_path/scripts
dataset_path=$global_path'/datasets'
list_path=$global_path'/lists'
temp_files=$global_path'/executions/temp_files'
ontology_file=$dataset_path'/mondo.obo'
export PATH=$scripts_path:$PATH
source ~soft_bio_267/initializes/init_pets
source ~soft_bio_267/initializes/init_python


if [ "$1" == "1" ]; then
	mkdir $list_path
	#MP:0000260 => MPO: abnormal angiogenesis
	#MP:0000259 => MPO: abnormal vascular development
	#HP:0002597 => HPO: Abnormality of the vasculature
	term=$2
	enc_term=`echo -e $term | sed 's/:/%3A/g'`
	wget 'https://api.monarchinitiative.org/api/bioentity/phenotype/'$enc_term'/diseases?rows=10000&facet=false&unselect_evidence=false&exclude_automatic_assertions=false&fetch_objects=false&use_compact_associations=false&direct_taxon=false' -O $list_path/$term.json
	get_items.py -i $list_path/$term.json > $list_path/$term.mondo 
	# get_items.rb -i $list_path/$term.json > $list_path/$term.mondo 
	get_disease_mondo.rb -i $list_path/$term.mondo -k 'Orphanet:[0-9]*' -s -f $ontology_file -S $temp_files/supp_mondo_orpha.txt -o $list_path/$term.dict
	cut -f 2 $list_path/$term.dict | sort -u | sed 's/Orphanet/ORPHA/g' > $dataset_path/$term.txt
fi

if [ "$1" == "2" ]; then

	wget https://archive.monarchinitiative.org/latest/tsv/gene_associations/gene_disease.9606.tsv.gz -O $temp_files'/gene_disease.9606.tsv.gz'
	gunzip $temp_files'/gene_disease.9606.tsv.gz'
	wget http://geneontology.org/gene-associations/goa_human.gaf.gz -O $temp_files'/goa_human.gaf.gz'
	gunzip $temp_files'/goa_human.gaf.gz'

fi

if [ "$1" == "3" ]; then
	go=$2 #GO:0001525 => Angiogenesis
	cut -f 2,5 $temp_files'/gene_disease.9606.tsv' > $temp_files'/gene_MONDO4go.txt'
	grep -w $go $temp_files'/goa_human.gaf' | awk '{if($7 != "IEA" && $7 != "NAS" && $7 != "ISS") print $3}' | sort -u > $list_path/go_gene.lst
	grep -w -F -f $list_path/go_gene.lst $temp_files'/gene_MONDO4go.txt' > $list_path'/gene_MONDO4go_filt.txt'
	cut -f 1 $list_path'/gene_MONDO4go_filt.txt' | sort -u > $list_path'/go_genes.txt'
	cut -f 2 $list_path'/gene_MONDO4go_filt.txt' | sort -u > $list_path'/go_mondo.txt'
	get_disease_mondo.rb -i $list_path'/go_mondo.txt' -k 'Orphanet:[0-9]*' -s -f $ontology_file -S $temp_files/supp_mondo_orpha.txt -o $list_path/go_query.dict
	cut -f 2 $list_path/go_query.dict | sort -u | sed 's/Orphanet/ORPHA/g' > $dataset_path/$go.txt
fi
