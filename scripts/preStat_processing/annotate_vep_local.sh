# annotate_vep_local.sh
/opt/vep/src/ensembl-vep/vep \
	--cache --species homo_sapiens --assembly GRCh38 \
	--dir_cache /opt/vep/.vep --dir_plugins /opt/vep/.vep/plugins \
	--json --canonical --everything --force_overwrite \
	-i input.tmp \
	-o vars-annotated.json \
	--plugin AlphaMissense,file=/plugin_resources/alphamissense/AlphaMissense_hg38.tsv.gz \
	--plugin CADD,file=/plugin_resources/CADD/whole_genome_SNVs.tsv.gz \
	--plugin mutfunc,db=/plugin_resources/mutfunc/mutfunc_data.db \
	--plugin EVE,file=/plugin_resources/EVE/eve_plugin/eve_merged.vcf.gz \
	--plugin MaveDB,file=/plugin_resources/maveDB/MaveDB_variants.tsv.gz \
	--plugin LoFtool,/plugin_resources/loftool/LoFtool_scores.txt \
	--plugin PrimateAI,/plugin_resources/primateAI/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz \
	--plugin REVEL,file=/plugin_resources/revel/new_tabbed_revel_grch38.tsv.gz \
	--plugin pLI,/plugin_resources/pli/pLI_values.txt \
	--plugin Blosum62 \
	--plugin Conservation,/plugin_resources/conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw