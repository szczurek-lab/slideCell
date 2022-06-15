#' perform Wilcoxon for potential geneA - geneB synthetic lethal pair
#' check if geneA altered cell lines are more dependent on geneB
#' @param sl_pair c(geneA_ENSEMBL, geneB_ENSEMBL) from a list of potential SL gene pairs
#' @param sens_matrix geneB dependency matrix in CCLE cell lines for a given knock-down experiment
#' @param alt_matrix geneA alteration matrix for CCLE cell lines
#' @param min_alt_no minimum no of cell lines which have geneA alteration, default 10
#' @return p-value of the Wilcoxon test
#' @export

spid <- function(sl_pair, sens_matrix1, alt_matrix1, plot=F, experiment_name="", min_alt_no = 10) {

  name_gene1 <- as.character(sl_pair[1])
  name_gene2 <- as.character(sl_pair[2])

  # checking if geneA is present in alteration matrix
  if (!(name_gene1 %in% colnames(alt_matrix1))) {
    print(paste(name_gene1, "is not present in alteration matrix for", experiment_name, "result is NA"))
    return(NA)

  # checking if geneB is present in sensitivity matrix
  } else if (!(name_gene2 %in% colnames(sens_matrix1))) {
    print(paste(name_gene2, "is not present in sensitivity matrix for ", experiment_name, "result is NA"))
    return(NA)
  }

  # list of cell lines suitable for testing - present in both matrices
  common_cell_lines <- intersect(rownames(sens_matrix1), rownames(alt_matrix1))

  # geneA alteration list
  alt_gene1 <- alt_matrix1[common_cell_lines, name_gene1]
  names(alt_gene1) <- rownames(alt_matrix1[common_cell_lines,])
  # geneB sensitivity list
  sens_gene2 = sens_matrix1[common_cell_lines,name_gene2]
  names(sens_gene2) <- rownames(sens_matrix1[common_cell_lines,])

  # list of cell line names with geneA alteration
  alt_cls <- names(alt_gene1[alt_gene1==1])

  # geneB sensitivity list in geneA altered cell lines
  sens_gene2_alt_gene1_1 = sens_gene2[alt_cls]

  if (length(sens_gene2_alt_gene1_1) < min_alt_no) {
    print("number of geneA altered cell lines is lower than required min")
    return(NA)
  }

  # list of cell line names withOUT geneA alteration
  not_alt_cls <- names(alt_gene1[alt_gene1==0])

  # geneB sensitivity list in geneA NOT altered cell lines
  sens_gene2_alt_gene1_0 = sens_gene2[not_alt_cls]

  res <- wilcox.test(sens_gene2_alt_gene1_1, sens_gene2_alt_gene1_0, alternative = c("less"))
  if (plot) {
    pl <- slideCell::plot_Wilcoxon_rank_sum_test(sl_pair, sens_matrix1=sens_matrix1, alt_matrix1=alt_matrix1, pval=signif(res$p.value,3), experiment=experiment_name)
    print(pl)
  }
  return(res$p.value)


}

#' plot boxplot for Wilcoxon for potential geneA - geneB synthetic lethal pair
#' check if geneA altered cell lines are more dependent on geneB
#' @param sl_pair c(geneA_ENSEMBL, geneB_ENSEMBL) from a list of potential SL gene pairs
#' @param sens_matrix1 geneB dependency matrix in CCLE cell lines for a given knock-down experiment
#' @param alt_matrix1 geneA alteration matrix for CCLE cell lines
#' @param experiment name of the experiment
#' @param pval p-value from the Wilcoxon test
#' @return boxplot of the Wilcoxon test
#' @export

plot_Wilcoxon_rank_sum_test <- function(sl_pair, sens_matrix1, alt_matrix1, experiment="", pval=F) {

    name_gene1=as.character(sl_pair[1])
    name_gene2=as.character(sl_pair[2])

    if (pval) {
        add_to_title = paste(" p_val =", pval)
    }
    else {
        add_to_title = ""
    }

    # list of cell lines suitable for testing - present in both matrices
    common_cell_lines <- intersect(rownames(sens_matrix1), rownames(alt_matrix1))

    # geneA alteration list
    alt_gene1 <- alt_matrix1[common_cell_lines, name_gene1]
    names(alt_gene1) <- rownames(alt_matrix1[common_cell_lines,])
    # geneB sensitivity list
    sens_gene2 = sens_matrix1[common_cell_lines,name_gene2]
    names(sens_gene2) <- rownames(sens_matrix1[common_cell_lines,])


    filter_G1=as.data.frame(alt_gene1==1)

    sens_1=as.data.frame(sens_gene2)
    sens_2=cbind(sens_1,filter_G1)
    colnames(sens_2)=c("sens","is_alt")

    filter2=is.na(sens_2$sens)
    sens_2=sens_2[!filter2,]

    #sens_wf=sens_2[order(sens_2$sens),]

    ggpubr::ggboxplot(sens_2, x = "is_alt", y = "sens",
    ylab = paste(name_gene2, "dependency score"), xlab = paste(name_gene1, "alteration"), add = "jitter", color = "is_alt", palette = c("grey", "red"), title = paste("Wilcoxon rank-sum test for", name_gene1, name_gene2, experiment, add_to_title)) + ggplot2::theme(legend.position="none") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, , size = 9))
}



#' verify list of potential SL gene pairs
#' performs Wilcoxon test for each pair
#' where geneA which is altered in the minimum no of cell lines for Wilcoxon test
#' and geneB is present in the chosen knock-down experiment
#' @param sl_pair_list list of geneA geneB pairs for verifing
#' @param sens_matrix list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package
#' @param data_set_names names for the dependency data sets - will be visible in results
#' @param alt_matrix1 geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package
#' @param min_alt_no minimum no of cell lines which have geneA alteration, default 10
#' @return list of Wilcoxon rank-sum test results separately for each dependency data set
#' @export

verify_sl_pairs <- function(sl_pair_list,
  sens_matrix=list(slideCell::sens_Drive_RSA, slideCell::sens_Drive_Ataris, slideCell::sens_Achilles_AC, slideCell::sens_Achilles_Demeter),
  data_set_names=c("Drive_RSA", "Drive_Ataris", "Achilles_AC", "Achilles_Demeter"),
  alt_matrix1 = slideCell::gene_cl_mut_matrix,
  min_alt_no = 10) {

  RSA = apply(sl_pair_list, 1, slideCell::spid,
    sens_matrix[[1]], alt_matrix1, min_alt_no = min_alt_no, experiment_name=data_set_names[1])

  Ataris = apply(sl_pair_list, 1, slideCell::spid,
    sens_matrix[[2]], alt_matrix1, min_alt_no = min_alt_no, experiment_name=data_set_names[2])

  AC = apply(sl_pair_list, 1, slideCell::spid,
    sens_matrix[[3]], alt_matrix1, min_alt_no = min_alt_no, experiment_name=data_set_names[3])

  Demeter = apply(sl_pair_list, 1, slideCell::spid,
    sens_matrix[[4]], alt_matrix1, min_alt_no = min_alt_no, experiment_name=data_set_names[4])

  results = list(sl_pair_list, RSA, Ataris, AC, Demeter)
  names(results) = c("sl_pair_list", data_set_names)

  return(results)


}
