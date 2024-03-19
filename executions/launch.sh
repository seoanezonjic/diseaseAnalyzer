#! /usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem='20gb'
#SBATCH --time='01:00:00'
#SBATCH --constraint=sd
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


current=`pwd`
global_path=$current/..
scripts_path=$global_path/scripts
dataset_path=$global_path'/datasets'
temp_files=$global_path'/executions/temp_files'
output_folder=$global_path'/executions/aRD_workflow'
list_path=$global_path'/lists'
dis_codes=$dataset_path'/HP:0000365.txt'
#orpha_codes=$dataset_path'/HP:0410008.txt'
#orpha_codes=$dataset_path'/hearImp_orpha_codes
#orpha_codes=$dataset_path'/nmd_orpha_codes'
#orpha_codes=$dataset_path'/HP:0000729'

export PATH=$scripts_path:$PATH
source ~soft_bio_267/initializes/init_python

mkdir -p $output_folder $dataset_path $temp_files

string_network=$dataset_path'/string_data.txt'

if [ "$1" == "1" ]; then
	echo 'Downloading files'
	
	wget https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz -O $dataset_path"/9606.protein.links.v11.5.txt.gz"
	gunzip $dataset_path'/9606.protein.links.v11.5.txt.gz'
	mv $dataset_path'/9606.protein.links.v11.5.txt' $string_network # copied from original execution
	
	### File with Human GeneID codes and STRING IDs (used as dictionary)
    wget https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz -O $dataset_path"/9606.protein.info.v11.5.txt.gz"
    gunzip $dataset_path"/9606.protein.info.v11.5.txt.gz"

    ### File with HPO phenotypes
	####wget http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt -O $temp_files"/genes_to_phenotype.txt"
	wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa -O $dataset_path'/phenotype.hpoa'
	
	touch $temp_files/dis_white_list
fi

if [ "$1" == "1a" ]; then
	term=$2
	if [ "$3" == "orpha" ]; then
		codename='Orphanet'
		codetag='ORPHA'
	elif [[ "$3" == "omim" ]]; then
		codename='OMIM'
		codetag='OMIM'
	fi
	#MONDO:0003847 Mendelian disease (filter by those diseases with a genetic component)
	echo 'Filtering list of phenotypes'
	semtools -C $2 -O MONDO > $temp_files/filtered_mondo_codes.txt
	semtools -i $temp_files/filtered_mondo_codes.txt --list -k $codename":[0-9]*" --xref_sense -O MONDO -o $temp_files/mondo_dis_codes.txt
	sed "s/$codename/$codetag/g" $temp_files/mondo_dis_codes.txt | cut -f 2 > $temp_files/dis_white_list
fi

if [ "$1" == "1b" ]; then
	if [ "$2" == "orpha" ]; then
		codename='Orphanet'
		codetag='ORPHA'
	elif [[ "$2" == "omim" ]]; then
		codename='OMIM'
		codetag='OMIM'
	fi
	echo 'preparing files...'
	intersect_columns -a $dis_codes -b $temp_files/dis_white_list > $temp_files/filtered_dis_codes.txt	
	grep -v '#' $dataset_path/phenotype.hpoa | grep -v -w 'NOT' | cut -f 1,2,4 > $temp_files/dis_name_phen.txt
	echo -e "DiseaseID\tHPOID" > $temp_files/disease_hpos.txt
	grep -w -F -f $temp_files/filtered_dis_codes.txt $temp_files/dis_name_phen.txt | cut -f 1,3 | sort -u | aggregate_column_data -i - -s '|' -x 1 -a 0 >> $temp_files/disease_hpos.txt	
	sed "s/$codetag/$codename/g" $temp_files/filtered_dis_codes.txt > $temp_files/disease_IDs
	echo $temp_files/disease_IDs
	# semtools with --list set is to get a dictionary
	echo 'executing semtools...'
	semtools -i $temp_files/disease_IDs --list -k "Orphanet:[0-9]*|OMIM:[0-9]*" -O MONDO -o $temp_files/disease_mondo_codes.txt
	sed -i "s/$codename/$codetag/g" $temp_files/disease_mondo_codes.txt
	cut -f 2 $temp_files/disease_mondo_codes.txt > $temp_files/mondo_codes.txt
	echo 'get genes for diseases from monarch'
	monarch_entities -i $temp_files/mondo_codes.txt -r disease-gene -a -o $temp_files/monarch_gene_disease_API # in progresss
	get_mondo_genes.py -i $temp_files/disease_mondo_codes.txt -m $temp_files/monarch_gene_disease_API -o $temp_files/disease_mondo_genes.txt
	cut -f 1,3 $temp_files/disease_mondo_genes.txt > $temp_files/disease_genes.txt
	echo 'extract STRING data'
	standard_name_replacer -i $string_network -I $dataset_path"/9606.protein.info.v11.5.txt" -c 1,2 -s ' ' | tr ' ' '\t' > $temp_files/string_transl_network.txt
fi

gene_filter_values=( 0 )
combined_score_filts=( 700 )
similarity_measures=( 'lin' )
#similarity_measures=( 'lin' 'resnik' 'jiang_conrath' )
min_groups=( 0 )
if [ "$1" == "2" ]; then
	source ~soft_bio_267/initializes/init_autoflow
	echo 'Launching analysis'
	for similarity_measure in "${similarity_measures[@]}"
	do	
		for min_group in "${min_groups[@]}"
		do
			for combined_score in "${combined_score_filts[@]}"
			do
				for gene_filter_value in "${gene_filter_values[@]}"
				do
					execution_name=$similarity_measure"_"$min_group"_"$combined_score"_"$gene_filter_value
					var_info=`echo -e "\\$similarity_measure=$similarity_measure,
					\\$string_network=$temp_files/string_transl_network.txt,
					\\$hub_zscore=2.5,
					\\$white_list=$list_path'/white_list',
					\\$string_dict=$dataset_path/9606.protein.info.v11.5.txt,
					\\$combined_score=$combined_score,
					\\$min_group=$min_group,
					\\$gene_filter=$gene_filter_value,
					\\$disease_gene_file=$temp_files/disease_genes.txt,
					\\$phenotype_annotation=$temp_files'/dis_name_phen.txt',
					\\$disease_mondo_genes=$temp_files'/disease_mondo_genes.txt',
					\\$disease_hpo_file=$temp_files'/disease_hpos.txt',
					\\$all_diseases=$temp_files'/filtered_dis_codes.txt',
					\\$report_template=$global_path/executions/templates/report_template.txt,
					\\$scripts_path=$scripts_path" | tr -d '[:space:]' `
					AutoFlow -w templates/aRD_analysis.txt -t '7-00:00:00' -m '100gb' -c 4 -o $output_folder"/"$execution_name -n 'sr' -e -V $var_info $2
				done
			done
		done
	done
elif [ "$1" == "2b" ]; then
	source ~soft_bio_267/initializes/init_autoflow
	for similarity_measure in "${similarity_measures[@]}"
	do	
		for min_group in "${min_groups[@]}"
		do
			for combined_score in "${combined_score_filts[@]}"
			do
				for gene_filter_value in "${gene_filter_values[@]}"
				do
					execution_name=$similarity_measure"_"$min_group"_"$combined_score"_"$gene_filter_value
					flow_logger -w -e $output_folder"/"$execution_name -r all $2
				done
			done
		done
	done
fi

#orpha_codes=$dataset_path'/HP:0000407.txt'
#orpha_codes=$dataset_path'/orphas_2012.txt' #antiguo: raquel_aRD_orpha_codes.txt
#orpha_codes='/mnt/home/users/pab_001_uma/pedro/proyectos/angio/results/orpha_codes'
