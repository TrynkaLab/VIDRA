# VIDRA: Variant-Informed Dose-Response Analysis

Code used in the manuscript: `Genetic dose-response modelling predicts drug mechanisms, dosing, and adverse events`

Key message of the manuscript:

The manuscript presents VIDRA (Variant-Informed Dose-Response Analysis), a statistical framework that models how genetic variation across the allele frequency and functional consequence spectrum informs gene-phenotype dose-response relationships. We integrated over 1.5 million associated variants from Open Targets data, to derive more than 88,000 dose-response-like relationships, increasing gene-phenotype associations 3-fold when incorporating rare variants. We benchmarked VIDRA against known drug targets to derive a Therapeutic Potential Score, identifying over 2,000 high-potential targets, most of which are currently untargeted. VIDRA supports target prioritisation, predicts direction of modulation, and extends to biomarker discovery and dose guidance, providing a versatile framework to inform multiple stages of drug development.

The nextflow pipeline perform the following analysis:

1. Ingest relevant information (i.e. common and rare variants)
2. Pre-process them to have a suitable structure for the Stan analyses
3. Perform VIDRA dose-response analysis per gene-phenotype
4. Analyse the results.

The Stan models are in the relevant directory: `stan_models`. To increase reproducibility, VIDRA pipeline was developed with a Docker container. it can be replicated with the standard command (locally after Docker installation). For more information on the methods we send to the manuscript. 

