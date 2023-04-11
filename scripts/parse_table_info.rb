#! /usr/bin/env ruby

##########################
## Script to parse info for supplementary tables
##########################

require 'optparse'
require 'csv'

##########################
#METHODS
##########################

def load_orpha_mondo_genes_file(input_file)
	orpha_genes_storage = {}
	File.open(input_file).each do |line|
		line.chomp!
		orpha_code, mondo_code, gene_code = line.split("\t")
		query = orpha_genes_storage[orpha_code]
		if query.nil?
			orpha_genes_storage[orpha_code] = [gene_code]
		else
			orpha_genes_storage[orpha_code] << gene_code
		end
	end
	return orpha_genes_storage
end

def load_orpha_hpo_file(input_file)
	orpha_hpos_storage = {}
	File.open(input_file).each do |line|
		line.chomp!
		next if line.include?('DiseaseID')
		orpha_code, hpos = line.split("\t")
		orpha_hpos = hpos.split('|')
		orpha_hpos_storage[orpha_code] = orpha_hpos
	end
	return orpha_hpos_storage
end

def load_orpha_clusters_file(input_file)
	orpha_clusters_storage = {}
	clusters_orpha_storage = {}
	File.open(input_file).each do |line|
		line.chomp!
		orpha_code, cluster_id = line.split("\t")
		orpha_clusters_storage[orpha_code] = cluster_id
		query = clusters_orpha_storage[cluster_id]
		if query.nil?
			clusters_orpha_storage[cluster_id] = [orpha_code]
		else
			query << orpha_code
		end
	end
	return orpha_clusters_storage, clusters_orpha_storage
end

def combine_files(orpha_genes_storage, orpha_hpos_storage, orpha_clusters_storage, clusters_orpha_transl_storage)
	final_table = {}
	orpha_codes = orpha_genes_storage.keys | orpha_hpos_storage.keys 
	orpha_codes.each do |orpha_code|
		genes_ary = orpha_genes_storage[orpha_code]
		hpos_ary = orpha_hpos_storage[orpha_code]
		cluster_id = orpha_clusters_storage[orpha_code]
		tranls_dis = clusters_orpha_transl_storage[orpha_code]
		unless hpos_ary.nil? || cluster_id.nil?
			genes_ary = ['-'] if genes_ary.nil?
			info_to_add = [tranls_dis, cluster_id, genes_ary.uniq.join(','), hpos_ary.uniq.join(',')]
			final_table[orpha_code] = info_to_add
		end
	end
	return final_table
end


def translate_orpha_diseases(orpha_dictionary_file, clusters_orpha_storage)
	clusters_orpha_transl_storage = {}
	# CSV.read(orpha_dictionary_file, col_sep: "\t").each do |line|
	# 	next if line.include?('#disease-db')
	# 	disease_name = line[2]
	# 	clusters_orpha_transl_storage[line[5]] = disease_name
	# end
	File.open(orpha_dictionary_file).each do |line|
		disease_code, disease_name, hpo_code = line.chomp.split("\t")
		clusters_orpha_transl_storage[disease_code] = disease_name
	end
	translated_diseases = {}
	clusters_orpha_storage.each do |cluster, diseases|
		diseases.each do |disease_id|
			disease_name = clusters_orpha_transl_storage[disease_id]
			query = translated_diseases[cluster]
			if query.nil?
				translated_diseases[cluster] = [disease_name]
			else
				query << disease_name
			end 
		end
	end
	return translated_diseases, clusters_orpha_transl_storage
end
##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:orpha_clusters_file] = nil
  opts.on("-c", "--orpha_clusters_file PATH", "Input file with ORPHA codes & cluster IDs") do |data|
    options[:orpha_clusters_file] = data
  end

  options[:orpha_dictionary] = nil
  opts.on("-d", "--orpha_dictionary PATH", "Input file with ORPHA codes & description. Used for translation.") do |data|
    options[:orpha_dictionary] = data
  end

  options[:orpha_mondo_genes_file] = nil
  opts.on("-i", "--orpha_mondo_genes_file PATH", "Input file with ORPHA, MONDO and gene codes") do |data|
    options[:orpha_mondo_genes_file] = data
  end

  options[:orpha_hpo_file] = nil
  opts.on("-f", "--orpha_hpo_file PATH", "Input file with ORPHA and HPO codes") do |data|
    options[:orpha_hpo_file] = data
  end

  options[:output_file] = 'output_table.txt'
  opts.on("-o", "--output_file PATH", "Output file") do |data|
    options[:output_file] = data
  end

  options[:output_clusters_orpha_table] = 'output_clusters_orpha_table.txt'
  opts.on("-O", "--output_clusters_orpha_table PATH", "Output file: table with clusters and diseases by cluster.") do |data|
    options[:output_clusters_orpha_table] = data
  end

end.parse!

##########################
#MAIN
##########################

orpha_genes_storage = load_orpha_mondo_genes_file(options[:orpha_mondo_genes_file])
orpha_hpos_storage = load_orpha_hpo_file(options[:orpha_hpo_file])
orpha_clusters_storage, clusters_orpha_storage = load_orpha_clusters_file(options[:orpha_clusters_file])
translated_diseases, clusters_orpha_transl_storage = translate_orpha_diseases(options[:orpha_dictionary], clusters_orpha_storage)

diseases_num_by_cluster = clusters_orpha_storage.values.map{|a| a.length}
ave_diseases_by_cluster = diseases_num_by_cluster.inject(0){|sum, a| sum + a}.fdiv(diseases_num_by_cluster.length)

final_table = combine_files(orpha_genes_storage, orpha_hpos_storage, orpha_clusters_storage, clusters_orpha_transl_storage)

File.open(options[:output_file], 'w') do |f|
	f.puts "DiseaseID\tDiseaseName\tClusterID\tGeneIDs\tHPOIDs"
	final_table.each do |orpha_code, info|
		f.puts "#{orpha_code}\t#{info.join("\t")}"
	end
end

File.open(options[:output_clusters_orpha_table], 'w') do |f|
	f.puts "ClusterID\tDiseaseIDs"
	translated_diseases.each do |clusterID, diseases|
		f.puts "#{clusterID}\t#{diseases.join(',')}"
	end
end

	
