#! /usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem='20gb'
#SBATCH --time='01:00:00'
#SBATCH --constraint=sd
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


framework_dir=`dirname $0`
export global_path=$(readlink -f $framework_dir )
config=$1

scripts_path=$global_path/scripts
dataset_path=$global_path'/datasets'
execution_path=$global_path'/executions'
source $config
temp_files=$output_folder'/temp_files'
list_path=$output_folder'/lists'
#dis_white_list=$output_folder/dis_white_list

export PATH=$scripts_path:$PATH
source ~soft_bio_267/initializes/init_python

mkdir -p $output_folder $dataset_path $temp_files $list_path


if [[ "$modules" =~ "d" ]]; then
	echo 'Downloading and processing source files'
	wget https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz -O $dataset_path"/9606.protein.links.v11.5.txt.gz"
    wget https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz -O $dataset_path"/9606.protein.info.v11.5.txt.gz" # File with Human GeneID codes and STRING IDs (used as dictionary)
   	gunzip $dataset_path/*.gz

	standard_name_replacer -i $dataset_path'/9606.protein.links.v11.5.txt' -I $dataset_path"/9606.protein.info.v11.5.txt" -c 1,2 -s ' ' | tr ' ' '\t' > $dataset_path/string_transl_network.txt

    ### File with HPO phenotypes
	####wget http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt -O $temp_files"/genes_to_phenotype.txt"
	wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa -O $dataset_path'/phenotype.hpoa'
	grep -v '#' $dataset_path/phenotype.hpoa | grep -v -w 'NOT' | cut -f 1,2,4 > $dataset_path/dis_name_phen.txt
	touch $list_path'/white_list'
fi

if [[ "$modules" =~ "HPO" ]]; then # GEnerate disease list from HPO
	echo 'Searching HPO'
	# if [ "$disease_database" == "orpha" ]; then
	# 	codename='Orphanet'
	# 	codetag='ORPHA'
	# elif [[ "$disease_database" == "omim" ]]; then
	# 	codename='OMIM'
	# 	codetag='OMIM'
	# fi
	monarch_entities -d $keyword_term -r phenotype-disease -o $list_path/$keyword_term'_disease_results.tab'
	cut -f 1 $list_path/$keyword_term'_disease_results.tab' | sort -u > $list_path/$keyword_term.mondo
	#semtools -i $list_path/$keyword_term.mondo --list -k $codename":[0-9]*" --xref_sense -O MONDO -o $list_path/$keyword_term'.dict'
	#cut -f 2 $list_path/$keyword_term'.dict' | sort -u | sed "s/$codename/$codetag/g" > $list_path/$keyword_term.txt
fi


if [[ "$modules" =~ "1a" ]]; then # Filter disease list and convert to OMIM or Orphanet
	echo "Filtering disease list"
	if [ "$disease_database" == "orpha" ]; then
		codename='Orphanet'
		codetag='ORPHA'
	elif [[ "$disease_database" == "omim" ]]; then
		codename='OMIM'
		codetag='OMIM'
	fi
	echo 'Filtering list of phenotypes'
	semtools -i $list_path/$keyword_term.mondo --list -O MONDO -F "p($mondo_white_list)n($mondo_black_list)" > $temp_files/filtered_mondo_codes.txt
	semtools -i $temp_files/filtered_mondo_codes.txt --list -k $codename":[0-9]*" --xref_sense -O MONDO -o $temp_files/mondo_dis_codes.txt
	sed "s/$codename/$codetag/g" $temp_files/mondo_dis_codes.txt | cut -f 2 | sort -u > $temp_files/filtered_dis_codes.txt
fi

if [[ "$modules" =~ "1b" ]]; then # Create disease-gene network
	echo "create network"
	if [ "$disease_database" == "orpha" ]; then
		codename='Orphanet'
		codetag='ORPHA'
	elif [[ "$disease_database" == "omim" ]]; then
		codename='OMIM'
		codetag='OMIM'
	fi
	echo 'preparing files...'
	#intersect_columns -a $dis_codes -b $dis_white_list > $temp_files/filtered_dis_codes.txt	
	echo -e "DiseaseID\tHPOID" > $temp_files/disease_hpos.txt
	
	echo "$dataset_path/dis_name_phen.txt"
	grep -w -F -f $temp_files/filtered_dis_codes.txt $dataset_path/dis_name_phen.txt | cut -f 1,3 | sort -u | aggregate_column_data -i - -s '|' -x 1 -a 0 >> $temp_files/disease_hpos.txt	
	sed "s/$codetag/$codename/g" $temp_files/filtered_dis_codes.txt > $temp_files/disease_IDs
	echo 'executing semtools...'
	semtools -i $temp_files/disease_IDs --list -k "Orphanet:[0-9]*|OMIM:[0-9]*" -O MONDO -o $temp_files/disease_mondo_codes.txt # semtools with --list set is to get a dictionary
	sed -i "s/$codename/$codetag/g" $temp_files/disease_mondo_codes.txt
	cut -f 2 $temp_files/disease_mondo_codes.txt > $temp_files/mondo_codes.txt
	echo 'get genes for diseases from monarch'
	monarch_entities -i $temp_files/mondo_codes.txt -r disease-gene -a -o $temp_files/monarch_gene_disease_API # in progresss
	get_mondo_genes.py -i $temp_files/disease_mondo_codes.txt -m $temp_files/monarch_gene_disease_API -o $temp_files/disease_mondo_genes.txt
	cut -f 1,3 $temp_files/disease_mondo_genes.txt | sort -u > $temp_files/disease_genes.txt
fi

workflow_folder=$output_folder/workflow
if [[ "$modules" =~ "2" ]]; then
	source ~soft_bio_267/initializes/init_autoflow
	echo 'Launching analysis'
	var_info=`echo -e "\\$similarity_measure=$similarity_measure,
	\\$string_network=$dataset_path/string_transl_network.txt,
	\\$hub_zscore=2.5,
	\\$white_list=$list_path'/white_list',
	\\$string_dict=$dataset_path/9606.protein.info.v11.5.txt,
	\\$combined_score=$combined_score,
	\\$min_group=$min_group,
	\\$gene_filter=$gene_filter_value,
	\\$disease_gene_file=$temp_files/disease_genes.txt,
	\\$phenotype_annotation=$dataset_path'/dis_name_phen.txt',
	\\$disease_mondo_genes=$temp_files'/disease_mondo_genes.txt',
	\\$disease_hpo_file=$temp_files'/disease_hpos.txt',
	\\$all_diseases=$temp_files/filtered_dis_codes.txt,
	\\$report_template=$global_path/templates/report_template.txt,
	\\$disease_ids=$temp_files/disease_IDs,
	\\$scripts_path=$scripts_path" | tr -d '[:space:]' `
	AutoFlow -w templates/aRD_analysis.txt -t '7-00:00:00' -m '100gb' -c 4 -o $workflow_folder -n 'sr' -e -V $var_info $2
elif [[ "$modules" =~ "2b" ]]; then
	flow_logger -w -e $workflow_folder -r all $2
fi
