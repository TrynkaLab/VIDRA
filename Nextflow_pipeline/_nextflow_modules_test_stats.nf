#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process AS_stats_gene_phenotypes_linear_bayes{
	input:
		tuple val(phenotype),
				val(gene),
				path(filePathsGT)
    output:
		tuple val(phenotype),
				val(gene),
				env(LOSS),
				env(TYPE),
				path('*.npy'),
				emit: AS_estimated_exp,
				optional: true
		path('*.npy'),
			emit: AS_estimated_exp_slim,
			optional: true
	"""
	python /container/scripts/stat_test_scripts/ASestimate_single_gene_linear_bayes.py \
		$phenotype \
		$gene \
		$filePathsGT 
	"""
}