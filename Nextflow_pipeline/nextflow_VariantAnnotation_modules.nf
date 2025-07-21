#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process coding_GWAS_effects{
	input:
		path cdGWAS_ProtineEffect
 	output:
		path "vars-annotated.json", emit: GWAS_variants_effect
	"""
 	# format input
	awk -F ',' '{print \$1, \$2, ".", \$3, \$4}' $cdGWAS_ProtineEffect | sort -h -k1,1 -k2,2 | sed -n '/^chr/!p' | uniq > input.tmp
	# This script run VEP with alphaMissense annotation 
	bash /container/scripts/annotate_vep_local.sh
	"""
}

process parse_coding_GWAS_effects{
	publishDir (
		path: "$HOME/$params.output_dir/fetched_data/variant_annotation/", 
		saveAs: {"coding_GWAS_variants_effects.csv"}, 
		mode: 'copy', 
		overwrite: true
	)
	input:
		path cdGWAS_ProtineEffect
 	output:
		path "parsed_vep.csv", emit: GWAS_variants_effect
	"""
	python /container/scripts/collate_ProteinInfo.py $cdGWAS_ProtineEffect
	"""
}

// Split Protvar info inro per gene files
process RareVars_AZ_protein_info_split {
	input:
		path rare_AZ_vars
    output:
    	path "rare_vars-*.csv", emit: RareVars_AZ_vars_PRE
	"""
	python /container/scripts/preprocess_get_ProteinInfo_AZ.py \
		$rare_AZ_vars \
		"Variant" \
		"-" 
	"""
}
// Fetch info from ProtVar and EBI
process RareVars_AZ_protein_info {
	input:
		path rare_AZ_vars
    output:
    	path "vars-annotated.json", emit: RareVars_AZ_vars_WORK
	"""
	# format input
	awk -F ',' '{print \$1, \$2, ".", \$3, \$4}' $rare_AZ_vars | sort -h -k1,1 -k2,2 | sed -n '/^chrom/!p' | uniq > input.tmp
	# This script run VEP with alphaMissense annotation 
	bash /container/scripts/annotate_vep_local.sh
	"""
}
// Put all Prot Var info together
process RareVars_AZ_protein_info_collate {
	publishDir (
		path: "$HOME/$params.output_dir/fetched_data/variant_annotation/", 
		saveAs: {"az_variants_effects.csv"}, 
		mode: 'copy', 
		overwrite: true
	)
	input:
		path rare_AZ_vars, stageAs: 'annotated_partition_??.tsv'
    output:
    	path "parsed_vep.csv", emit: AZ_variants_effect
	"""
	python /container/scripts/collate_ProteinInfo.py $rare_AZ_vars 
	"""
}

// ClinVar
// Split Protvar info inro per gene files
process RareVars_ClinVar_protein_info_split {
	input:
		path rare_ClinVar_vars
    output:
    	path "rare_vars-*.csv", emit: RareVars_ClinVar_vars_PRE
	"""
	python /container/scripts/preprocess_get_ProteinInfo_rare_ClinVar_vars.py \
		$rare_ClinVar_vars \
		"variantId" \
		"_" 
	"""
}

// Fetch info from ProtVar and EBI
process RareVars_ClinVar_protein_info {
	input:
		path rare_ClinVar_vars
    output:
    	path "vars-annotated.json", emit: RareVars_ClinVar_Protein_WORK
	"""
	# format input
	sort -h -k1,1 -k2,2 $rare_ClinVar_vars | uniq > input.tmp
	# This script run VEP with alphaMissense annotation 
	bash /container/scripts/annotate_vep_local.sh
	"""
}

// Put all Prot Var info together
process RareVars_ClinVar_protein_info_collate {
	publishDir (
		path: "$HOME/$params.output_dir/fetched_data/variant_annotation/", 
		saveAs: {"clinvar_variants_effects.csv"}, 
		mode: 'copy', 
		overwrite: true
	)
	input:
		path rare_ClinVar_vars, stageAs: 'annotated_partition_??.tsv'
    output:
    	path "parsed_vep.csv", emit: CV_variants_effect
	"""
	python /container/scripts/collate_ProteinInfo.py $rare_ClinVar_vars 
	"""
}
