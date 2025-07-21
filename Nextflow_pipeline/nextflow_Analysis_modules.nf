#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process AS_stats_gene_phenotypes_linear{
	input:
		tuple val(gene), path(filePathsGT)
    output:
		path('*_Bayes_hierarchical.csv'),
			emit: dummy_emit,
			optional: true
	"""
	python /container/scripts/ASestimate_single_gene_linear_bayes.py \
		0.1 \
		$gene \
		$filePathsGT
	"""
}

process endsignal{
	input:
		val(dummy_input)
	output:
		path ('endSignal.txt')
	"""
	touch endSignal.txt
	"""
}

process compactCsvASstats{
	input:
    	path(csv_files)
    output:
    	file 'AScollated_${UUID.randomUUID()}.csv'
	"""
	cat $csv_files > ASCollated_${UUID.randomUUID()}.csv
	rm $csv_files
	"""
}

process AS_stats_test_hyperparameters{
	input:
		tuple val(h1),
				val(phenotype),
				val(gene),
				path(filePathsGT)
    output:
		path('*_Bayes_testHyperp.csv'),
			emit: AS_estimated_lin_slim,
			optional: true
	"""
	python /container/scripts/ASestimate_testHyperp_single_gene_linear_bayes.py \
		$h1 \
		$phenotype \
		$gene \
		$filePathsGT
	"""
}
