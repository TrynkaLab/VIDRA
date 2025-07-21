#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// this is to get the clusters definition 
process cluster_to_genes {
    output:
		path "phenotypes.csv", emit: cluster_phenotypes
	"""
	python /container/scripts/extract_trait_id_from_cluster_file.py $params.cluster_file $params.subset_cluster
	"""
}

// Once one has the phenotypes, we need to extract genes associated with them.
// Genes extraction comes from colocalisation for common traits

// We get GWAS studies from this script
process link_pheno_to_study{
	publishDir "$HOME/Resources/", mode: 'copy'
	input:
		path cluster_phenotypes
 	output:
		path "gwas_study_info.csv"
		path "trait_to_gwas.csv"
	"""
	python /container/scripts/link_traits_to_GWAS.py $cluster_phenotypes
 	"""
}

// Get relevant variants using coloc evidences
// i.e. For the relevant studies what variants are colocalising
process coloc_variants{
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path relevant_gwas
 	output:
		path "colocalising_variants.csv", emit: colocalsing_variants
	"""
	python /container/scripts/get_coloc_vars.py $relevant_gwas $params.colocalisation_threshold
 	"""
}

// Per gene process
process coloc_variants_per_gene{
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path genes
 	output:
		path "colocalising_variants.csv", emit: colocalsing_variants
	"""
	python /container/scripts/get_coloc_vars_per_gene.py $genes $params.colocalisation_threshold
 	"""
}

// This script collect the varinats that are coding and have been identified in GWAS studies.
// Becasue these are coding variants, I can extract the protein info from ProtVar and EBI and do not need to do the same QTL effect on this vairiants
process coding_GWAS{ 
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy', pattern: "GWAS_coding_variants_NoColocalising.csv"
 	input:
		path genes_to_fetch
	output:
		path "GWAS_coding_variants_NoColocalising.csv", emit: coding_GWAS_variants_with_disease_effect
		path "GWAS_coding_variants_fromOT.csv", emit: coding_GWAS_variants_toFetchProtineEffect
	"""
	python /container/scripts/get_coding_GWAS_nonColoc.py $genes_to_fetch
 	"""
}


// Once one has the info on the colocalising variants, could extract the variant info from the respective GWAS and QTL experiments
// GWAS
process GWAS_variants_extraction{
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path colocalsing_variants
 	output:
		path "GWAS_variants.csv", emit: GWAS_variants
	"""
	python /container/scripts/get_GWAS_vars.py $colocalsing_variants 
 	"""
}
// QTL
process QTL_variants_extraction{
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path colocalsing_variants
 	output:
		path "QTL_variants.csv", emit: QTL_variants
	"""
	python /container/scripts/get_QTL_vars.py $colocalsing_variants 
 	"""
}

// In parallel to common variats we extract the variants that have been reported in clinvar or AZ. 
// Can be done in parallel becasue they don't need coloc evidences
// AZ
process Rare_AZvariants_extraction {
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path cluster_phenotypes
    output:
    	path 'az_variants.csv', emit: RareVars_AZ_vars
	"""
	python /container/scripts/get_AZ_vars.py $cluster_phenotypes $params.AZ_pheno_conversion_Table $params.AZ_data
	"""
}

process Rare_AZvariants_extraction_per_gene {
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path genes_symbols
    output:
    	path 'az_variants.csv', emit: RareVars_AZ_vars
	"""
	python /container/scripts/get_AZ_vars_per_gene.py $genes_symbols $params.AZ_pheno_conversion_Table $params.AZ_data
	"""
}

// AZ Burden tests
process AZburdens_extraction {
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path cluster_phenotypes
    output:
    	path 'burden_tests.csv', emit: burdens_AZ
	"""
	python /container/scripts/get_AZ_burden_tests.py $cluster_phenotypes
	"""
}

process AZburdens_extraction_per_gene {
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path genes_symbols
    output:
    	path 'burden_tests.csv', emit: burdens_AZ
	"""
	python /container/scripts/get_AZ_burden_tests_per_gene.py $genes_symbols
	"""
}
// ClinVar
process ClinVar_variants_extraction{
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path QTL_variants
 	output:
		path "clinvar_variants.csv", emit: clinvar_variants
	"""
	python /container/scripts/get_clinvar_vars.py $QTL_variants 
 	"""
}
process ClinVar_variants_extraction_per_gene{
	publishDir "$HOME/$params.output_dir/fetched_data/", mode: 'copy'
	input:
		path genes
 	output:
		path "clinvar_variants.csv", emit: clinvar_variants
	"""
	python /container/scripts/get_clinvar_vars_per_gene.py $genes 
 	"""
}
// Put together all the variants and split by gene-phenotype pair - pre-processing requiredd for the VIDRA workflow
// This script combines all the info fetched, and structure them in a shape that is useful for Stan script analyses
// output to a series of parquet files is convenient because it allow to easily split the combination gene-pheno to different files
process pre_processing_AS {
	publishDir "$HOME/pre_AS/$params.output_dir/bayes/", mode: 'move',  overwrite: true, pattern: "*_bayesDF/**"
	input:
		path coloc_anchor
		path gwas_vars
		path qtl_vars
		path pheno_link
		path AZdisease
		path AZprotein
		path ClinVarDisease
		path ClinVarProtein
		path az_burden
		path cdGWAS
		path cdGWASprotein
    output:
		path "${params.pre_AS_out_name}*/**", emit: merged_harmosised_data
	"""
	python /container/scripts/pre_processing_AS_per_gene_pheno.py \
		$coloc_anchor \
		$gwas_vars \
		$qtl_vars \
		$pheno_link \
		$AZdisease \
		$AZprotein \
		$ClinVarDisease \
		$ClinVarProtein \
		$params.pre_AS_out_name \
		$az_burden \
		$cdGWAS \
		$cdGWASprotein \
		$params.genes_to_fetch_complete \
		$params.convertOldEFOtoUpdateOnes
	"""
}
