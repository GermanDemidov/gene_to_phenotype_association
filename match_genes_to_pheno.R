#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

file_name_with_variants = args[1]

### Master ontology created

library(ontologyIndex)
library(ontologySimilarity)

data(hpo)

samples_list <- read.csv(file_name_with_variants, stringsAsFactors = F, sep="\t")

samples_list = samples_list[which(samples_list$disease_details_HPO_term_id != ""),]

list_of_phenotypes = list()
list_of_submitters = list()

for (i in 1:nrow(samples_list)) {
  hpos_list <- strsplit(samples_list[i,]$disease_details_HPO_term_id, "; ")[[1]]
  if (length(hpos_list) > 0) {
    hpo_numbers = c()
    for (elem in hpos_list) {
      hpo_num = strsplit(elem, " - ")[[1]][1]
      if (hpo_num %in% hpo$id) {
        hpo_numbers = c(hpo_numbers, hpo_num)
      }
    }
    if (length(hpo_numbers) > 0) {
      list_of_phenotypes[[samples_list[i,]$sample]] = hpo_numbers
    }
  }
}

information_content <- descendants_IC(hpo)




### Gene to phenotype list
gene_to_pheno <- read.table("genes_to_phenotype.txt", sep="\t", quote="", stringsAsFactors = F)
# HP:0040283 - occasional (29-5%), HP:0040282 = frequent (30-79%), HP:0040284 = very rare (4-1%), HP:0040281 = very frequent (80-99%), HP:0040280 - obligate


### Dataset to analyze
gene_list <- samples_list$gene


background_similarities <- list()
all_samples = unique(samples_list$sample)
for (sample in all_samples) {
  print(sample)
  if (sample %in% names(list_of_phenotypes) & length(list_of_phenotypes[[sample]]) > 0) {
    vec_of_phenos = c()
    phenotype_patient_set = list_of_phenotypes[[sample]]
    if (length(phenotype_patient_set) > 0) {
      for (gene in unlist(gene_list)) {
        hpos_for_gene = list()
        phenotype_matching = -0.1
        diseases <- unique(gene_to_pheno[which(gene_to_pheno[,2] == gene), 9])
        for (disease in diseases) {
            hpos_for_gene[[disease]] = gene_to_pheno[which(gene_to_pheno[,2] == gene & gene_to_pheno[,9] == disease), 3]
            hpos_for_gene[[disease]] = hpos_for_gene[[disease]][which(hpos_for_gene[[disease]] %in% hpo$id)]
            
            set_of_hpos <- minimal_set(hpo, terms=hpos_for_gene[[disease]])
            lst_hpo = list(phenotype_patient_set, set_of_hpos)
            sim_mat <- get_sim_grid(ontology=hpo, information_content=information_content, term_sets=lst_hpo)
            vec_of_phenos <- c(vec_of_phenos, sim_mat[1,2])
        }
      }
      background_similarities[[sample]] = vec_of_phenos
    }
  }
}

best_gene_similarities = rep(-0.1, nrow(samples_list))
phevalues = rep(1.1, nrow(samples_list))
for (i in 1:nrow(samples_list)) {
  sample_ID = samples_list[i,]$sample
  if (sample_ID %in% names(list_of_phenotypes) & length(list_of_phenotypes[[sample_ID]]) > 0) {
    best_gene_similarities[i] = 0
    phenotype_patient_set = list_of_phenotypes[[sample_ID]]
    
    list_of_genes = unlist(samples_list[i,]$gene)
    whole_gene_set = c()
    for (gene in list_of_genes) {
      diseases <- unique(gene_to_pheno[which(gene_to_pheno[,2] == gene), 9])
      phenotype_matching = -0.1
      
      for (disease in diseases) {
        hpos_for_gene[[disease]] = gene_to_pheno[which(gene_to_pheno[,2] == gene & gene_to_pheno[,9] == disease), 3]
        hpos_for_gene[[disease]] = hpos_for_gene[[disease]][which(hpos_for_gene[[disease]] %in% hpo$id)]
        set_of_hpos <- minimal_set(hpo, terms=hpos_for_gene[[disease]])
        
        lst_hpo = list(phenotype_patient_set, set_of_hpos)
        sim_mat <- get_sim_grid(ontology=hpo, information_content=information_content, term_sets=lst_hpo)
        phenotype_matching = max(phenotype_matching, sim_mat[1,2])
      }
      
      best_gene_similarities[i] = max(best_gene_similarities[i], phenotype_matching)
      phevalues[i] = length(which(background_similarities[[sample_ID]] >= best_gene_similarities[i])) / length(background_similarities[[sample_ID]])
    }
  }
}

pdf("scores_of_variants.pdf")

plot(phevalues, best_gene_similarities, col=rgb(0,0,0,0.2), pch=19, main="phenotype matching vs phe-values", xlab="phe-value", ylab="Similarity" )

dev.off()





write.table(cbind(phevalues, best_gene_similarities, samples_list), quote=F, file="phevalues.annot.tsv", row.names=F, sep="\t")






