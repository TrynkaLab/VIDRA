#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { cluster_to_genes; coding_GWAS; link_pheno_to_study; coloc_variants; coloc_variants_per_gene } from './nextflow_modules.nf'
include { GWAS_variants_extraction; QTL_variants_extraction; Rare_AZvariants_extraction; Rare_AZvariants_extraction_per_gene } from './nextflow_modules.nf'
include { AZburdens_extraction; AZburdens_extraction_per_gene } from './nextflow_modules.nf'
include { ClinVar_variants_extraction; ClinVar_variants_extraction_per_gene } from './nextflow_modules.nf'
// // Variant annotation modules
include { coding_GWAS_effects; parse_coding_GWAS_effects } from './nextflow_VariantAnnotation_modules.nf'
include { RareVars_AZ_protein_info_split; RareVars_AZ_protein_info; RareVars_AZ_protein_info_collate } from './nextflow_VariantAnnotation_modules.nf'
include { RareVars_ClinVar_protein_info_split; RareVars_ClinVar_protein_info; RareVars_ClinVar_protein_info_collate } from './nextflow_VariantAnnotation_modules.nf'
// Parse modules
include { pre_processing_AS } from './nextflow_modules.nf'
// Analysis modules
include { plotAS; create_scrambled_ref_distribution } from './nextflow_modules.nf'
include { compactCsvASstats; AS_stats_gene_phenotypes_linear; AS_stats_gene_phenotypes_sigmoid; AS_stats_test_hyperparameters; viewASmtx; compareASstats } from './nextflow_Analysis_modules.nf'
// Comparison modules
include { calculate_alternativeMR } from './nextflow_comparison_modules.nf'

// The entry points is getting information on all the variants starting from gene IDs
workflow get_data_from_genes{
	main:
		genes = channel.fromPath(params.genes_to_fetch)
		gene_symbols = channel.fromPath(params.gene_symbols_to_fetch)
		coding_GWAS(genes)
		coloc_variants_per_gene(genes)
		GWAS_variants_extraction(coloc_variants_per_gene.out.colocalsing_variants)
		QTL_variants_extraction(coloc_variants_per_gene.out.colocalsing_variants)
		Rare_AZvariants_extraction_per_gene(gene_symbols)
		AZburdens_extraction_per_gene(gene_symbols)
		ClinVar_variants_extraction_per_gene(genes)
		// link_pheno_to_study() // This needs to be run only once when need to link phenotypes to studies
}

// Annotate variants that are coding (i.e most do not have QTLs)
workflow annotate_variants{
	main:
		// Files to annotate
		gwasVar = channel.fromPath("$HOME/$params.output_dir/fetched_data/GWAS_coding_variants_NoColocalising.csv")
		AZrareVar = channel.fromPath("$HOME/$params.output_dir/fetched_data/az_variants.csv")
		ClinVar = channel.fromPath("$HOME/$params.output_dir/fetched_data/clinvar_variants.csv")
		// Annotation scripts
		// GWAS coding variants - not colocalising
		gwasVar \
		| coding_GWAS_effects \
		| parse_coding_GWAS_effects
		// AZ phewas variants
		AZrareVar \
		| RareVars_AZ_protein_info_split \
		| flatten \
		| RareVars_AZ_protein_info \
		| collect \
		| RareVars_AZ_protein_info_collate
		// ClinVar variants
		ClinVar \
		| RareVars_ClinVar_protein_info_split \
		| flatten \
		| RareVars_ClinVar_protein_info \
		| collect \
		| RareVars_ClinVar_protein_info_collate
}

// Collate all the information into a single file 
// Pre-processing step required to run the AS analysis stats
workflow parse_data{
	main:
		// Get all the files needed in this workflow module
		// These are fetched and generated in the previous workflow module (i.e. get_data and annotate_variants)
		// common non-coding variants
		coloc_anchor = channel.fromPath("$HOME/$params.output_dir/fetched_data/colocalising_variants.csv")
		gwas_vars = channel.fromPath("$HOME/$params.output_dir/fetched_data/GWAS_variants.csv")
		qtl_vars = channel.fromPath("$HOME/$params.output_dir/fetched_data/QTL_variants.csv")
		pheno_link = channel.fromPath("$HOME/Resources/trait_to_gwas.csv") 
		// rare variants
		// AZ
		AZdisease = channel.fromPath("$HOME/$params.output_dir/fetched_data/az_variants.csv")
		AZprotein = channel.fromPath("$HOME/$params.output_dir/fetched_data/variant_annotation/az_variants_effects.csv")
		// ClinVar
		ClinVarDisease = channel.fromPath("$HOME/$params.output_dir/fetched_data/clinvar_variants.csv")
		ClinVarProtein = channel.fromPath("$HOME/$params.output_dir/fetched_data/variant_annotation/clinvar_variants_effects.csv")
		// Burden tests from OT
		az_burden = channel.fromPath("$HOME/$params.output_dir/fetched_data/burden_tests.csv")
		// Common coding variants from GWAS
		cdGWAS = channel.fromPath("$HOME/$params.output_dir/fetched_data/GWAS_coding_variants_NoColocalising.csv")
		cdGWASprotein = channel.fromPath("$HOME/$params.output_dir/fetched_data/variant_annotation/coding_GWAS_variants_effects.csv")
		// Read in the info and split by gene-disease pair
		pre_processing_AS( 
			coloc_anchor,
			gwas_vars,
			qtl_vars,
			pheno_link,
			AZdisease,
			AZprotein,
			ClinVarDisease,
			ClinVarProtein,
			az_burden,
			cdGWAS,
			cdGWASprotein
			) 
}

// Stats analysis
// Compiling the stan models doens't work in the current version of nextflow pipeline,
// So when there is a new model to test, it has to be compiled manually before running the pipeline
// i.e. run the following commands manually: 
// path_to_model = '/container/stan_models/hierarchical_mix_models_robust_beta.stan'
// model = CmdStanModel(stan_file=path_to_model,cpp_options={'STAN_THREADS':'true'})

workflow analyse_data{
	main:
		// Define the list of genes I want to prioritise for manual curation
		allowed_genes = file(params.filter_ASanalysis_genes).readLines()
		scramble_vars = Channel.fromPath("$HOME/scrambled_variants/scrambles.csv")
		// Gene-trait pair file
		Channel.fromPath("$HOME/pre_AS/output_all_genes_20240320/bayes/preASwrangle_bayesDF_BLOOD/as_gene*", type: 'any', maxDepth: 3)
				.map { path -> 
						def components = path
											.toString() // convert to string
											.replace('as_gene=', '') // remove the prefixes
											.split('/') // split by the path
						return [components[-1], path] // return the gene, trait and path in a list
						}.set { AS_parsed_GT }
		// Gene-trait pair file - all are in the same gene file			
		AS_parsed_GT
			| AS_stats_gene_phenotypes_linear
}

// Test hyperparameters
workflow tets_hyperparameters{
	main:
		// Gene-trait pair file
		Channel.fromPath("$HOME/pre_AS/bayes/preASwrangle_bayesDF/**/*.parquet", type: 'any', maxDepth: 3)
				.map { path -> 
						def components = path
											.toString() // convert to string
											.replace('as_gene=', '') // remove the prefixes
											.replace('as_disease=', '')
											.split('/') // split by the path
						return [components[-3], components[-2], path] // return the gene, trait and path in a list
						}.set { AS_parsed_GT }
		// Gene-trait pair file		
		AS_parsed_GT.groupTuple( by: [0,1] ).set { as_gt }
		// AS object
		// Create hyperparameters references
		hyp1 = Channel.from(.01,.05,.1,.2)
		hyp1\
			.combine( as_gt )\
			.set{ hyperparameters }
		// Run the test
		hyperparameters
			| AS_stats_test_hyperparameters
}

workflow {
	main:
		get_data_from_genes() // alternative entry point to the previous one
		annotate_variants()
		parse_data()
		analyse_data()	
		stat_assess()
}