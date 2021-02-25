# gene_to_phenotype_association

This is a simple wrapper around R packages <a href="https://cran.r-project.org/web/packages/ontologyIndex/index.html">ontologyIndex</a> and <a href="https://cran.r-project.org/web/packages/ontologySimilarity/index.html">ontologySimilarity</a>.

This script takes a table of genomic variants in some patients as an input and outputs 2 metrics for them:

1) raw similarity between a gene and patient's phenotype (0 = no similarity, 1 = perfect match)
2) phe-value - how well this particular gene explains the patient's phenotype, in comparison with other genes affected in this dataset (0 = this gene is the best match for this patient phenotype, 1 = this is absolutely random match and other genes explain the phenotype much better).

Requirements:
1) both abovementioned R packages installed,
2) genes_to_phenotype file from <a href="https://hpo.jax.org/app/download/annotation">OMIM</a> is located in the same folder as this script,
3) your input file "affected genes in samples" is tab-separated, each affected gene-sample takes a separate row, and your file has a header with the column names "disease_details_HPO_term_id" (here HPO terms OF PATIENTS should be described, separated by "; " split, for example - "HP:0000248 - Brachycephaly; HP:0000343 - Long philtrum;", without quotation marks!), "sample" (containing sample IDs - take care that same IDs belong to patients with the same HPO described), "gene" (gene ID used in OMIM).

You run the script as:

'Rscript match_genes_to_pheno.R input_file.txt'

Result file with the name 'phevalues.annot.tsv' will appear in the directory with the script.
