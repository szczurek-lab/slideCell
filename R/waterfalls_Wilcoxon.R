#' checks primary site of a cell line
#' where geneA which is altered in the minimum no of cell lines for Wilcoxon test
#' and geneB is present in the chosen knock-down experiment
#' @param cl_id cell line id
#' @param cl_ann cell line annotations
#' @return primary site
#' @export

primary_site <- function(cl_id, cl_ann) {
  as.character(cl_ann[cl_ann[,1]==cl_id,5])
}

#' arranges single legend
#' @return grid.grab object
#' @export

grid_arrange_shared_legend <- function(..., title = "title", ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplot2::ggplotGrob(plots[[1]] +
                    ggplot2::theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- grid::unit(0.2, "npc") # percentage of screen for the legend
  gl <- lapply(plots, function(x) x +
                 ggplot2::theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)



  combined <- switch(position,
                     "bottom" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                            legend,ncol = 1,
                                            heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight)),
                     "right" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth)))


  grid::grid.newpage()

  grid::grid.draw(combined)
  #grid::grid.text(title, x=0.9, y=0.9, gp = grid::gpar(fontsize =7))
  #vwhb

  grid::grid.grab()


}

#' creates one waterfall plot
#' where geneA altered cell lines (for Wilcoxon test subpopulation) are colored
#' and geneB depenency is plotted
#' @param sl_pair geneA geneB pair for verifing
#' @param sens list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package
#' @param alt geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package
#' @param experiment name of the dependency data set
#' @return waterfall plot
#' @export

wf <- function (sl_pair, sens, alt, experiment, pval=F, pallete=slideCell::sites_colors, cl_annotations=slideCell::cl_annotations) {

  name_gene1=as.character(sl_pair[1])
  name_gene2=as.character(sl_pair[2])

  if (pval) {
      add_to_title = paste(" p_val =", pval)
  }
  else {
      add_to_title = ""
  }

  # list of cell lines suitable for testing - present in both matrices
  common_cell_lines <- intersect(rownames(sens), rownames(alt))

  # geneA alteration list
  alt_gene1 <- alt[common_cell_lines, name_gene1]
  names(alt_gene1) <- rownames(alt[common_cell_lines,])
  # geneB sensitivity list
  sens_gene2 = sens[common_cell_lines,name_gene2]
  names(sens_gene2) <- rownames(sens[common_cell_lines,])


  filter_G1=as.data.frame(alt_gene1==1)

  sens_1=as.data.frame(sens_gene2)
  sens_2=cbind(sens_1,filter_G1)
  colnames(sens_2)=c("sens","is_alt")

  filter2=is.na(sens_2$sens)
  sens_2=sens_2[!filter2,]

  sens_wf=sens_2[order(sens_2$sens),]

  filter3 = sens_wf$is_alt == T
  primary_site=rep("cl_not_altered",nrow(sens_wf))
  primary_site[filter3] = sapply(rownames(sens_wf[filter3,]), slideCell::primary_site, cl_annotations)

  sens_wf <- as.data.frame(cbind(cell_lines=seq(1,nrow(sens_wf),1),cl=rownames(sens_wf),sens_wf,primary_site))

  return(
  ggplot2::ggplot(data=sens_wf, ggplot2::aes(x=cell_lines, y=sens, fill = primary_site)) +
    ggplot2::geom_bar(stat="identity")+ ggplot2::scale_fill_manual(values = pallete[as.character(sens_wf$primary_site)])+
    ggplot2::ggtitle(paste(experiment, add_to_title)) +
    ggplot2::labs(x=paste0("cell lines, coloured when ", name_gene1, " altered"), y = paste0("dependency score ", name_gene2)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 8), axis.text.x = ggplot2::element_text(size=6),
          axis.title=ggplot2::element_text(size=7), legend.title = ggplot2::element_text(size =7), legend.key.size = grid::unit(0.25, "cm"),
          legend.text = ggplot2::element_text(size = 6))
          )

  #dev.off()
}



#' creates 4 waterfall plots
#' where geneA altered cell lines (for Wilcoxon test subpopulation) are colored
#' and geneB depenency is plotted
#' @param sl_pair list of geneA geneB pairs for verifing
#' @param sens_matrix list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package
#' @param data_set_names names for the dependency data sets - will be visible in results
#' @param alt_matrix1 geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package
#' @return 4 waterfall plots
#' @export
plot_many_waterfalls <- function(sl_pair,
sens_matrix=list(slideCell::sens_Drive_RSA, slideCell::sens_Drive_Ataris, slideCell::sens_Achilles_AC, slideCell::sens_Achilles_Demeter),
data_set_names=c("Drive_RSA", "Drive_Ataris", "Achilles_AC", "Achilles_Demeter"),
alt_matrix1 = slideCell::gene_cl_mut_matrix) {

    is_x1=F
    is_x2=F
    is_x3=F
    is_x4=F
    G1 = as.character(sl_pair[1])
    G2 = as.character(sl_pair[2])


    if (G1 %in% colnames(alt_matrix1) & G2 %in% colnames(sens_matrix[[1]])) {
        pv1 = signif(slideCell::spid(sl_pair, sens_matrix1 = sens_matrix[[1]], alt_matrix1=alt_matrix1), digits=2)
        x1 = slideCell::wf(sl_pair, sens = sens_matrix[[1]], alt = alt_matrix1, experiment=data_set_names[1], pval=pv1)
        is_x1=T

    }

    if (G1 %in% colnames(alt_matrix1) & G2 %in% colnames(sens_matrix[[2]])) {
        pv2 = signif(slideCell::spid(sl_pair, sens_matrix1 = sens_matrix[[2]], alt_matrix1=alt_matrix1), digits=2)
        x2 = slideCell::wf(sl_pair, sens = sens_matrix[[2]], alt = alt_matrix1, experiment=data_set_names[2], pval=pv2)
        is_x2=T
    }

    if (G1 %in% colnames(alt_matrix1) & G2 %in% colnames(sens_matrix[[3]])) {
        pv3 = signif(slideCell::spid(sl_pair, sens_matrix1 = sens_matrix[[3]], alt_matrix1=alt_matrix1), digits=2)
        x3 = slideCell::wf(sl_pair, sens = sens_matrix[[3]], alt = alt_matrix1, experiment=data_set_names[3], pval=pv3)
        is_x3=T
    }


    if (G1 %in% colnames(alt_matrix1) & G2 %in% colnames(sens_matrix[[4]])) {
        pv4 = signif(slideCell::spid(sl_pair, sens_matrix1 = sens_matrix[[4]], alt_matrix1=alt_matrix1), digits=2)
        x4 = slideCell::wf(sl_pair, sens = sens_matrix[[4]], alt = alt_matrix1, experiment=data_set_names[4], pval=pv4)
        is_x4=T
    }

    if(is_x4){
        x = grid_arrange_shared_legend(x1,x2,x3,x4, ncol = 2, nrow = 2, position = "right")

    }

    if(!is_x4){
        x = grid_arrange_shared_legend(x1,x2,x3, ncol = 2, nrow = 2, position = "right")

    }


    dev.off()

}

#' creates 4 waterfall plots
#' where geneA altered cell lines (for Wilcoxon test subpopulation) are colored
#' and geneB depenency is plotted
#' @param sl_pair list of geneA geneB pairs for verifing
#' @param sens_matrix list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package
#' @param data_set_names names for the dependency data sets - will be visible in results
#' @param alt_matrix1 geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package
#' @return 4 waterfall plots
#' @export
plot_many_waterfalls_old <- function(sl_pair,
sens_matrix=list(slideCell::sens_Drive_RSA, slideCell::sens_Drive_Ataris, slideCell::sens_Achilles_AC, slideCell::sens_Achilles_Demeter),
data_set_names=c("Drive_RSA", "Drive_Ataris", "Achilles_AC", "Achilles_Demeter"),
alt_matrix1 = slideCell::gene_cl_mut_matrix) {

    e1 <- slideCell::wf(sl_pair, sens = sens_matrix[[1]], alt = alt_matrix1, experiment=data_set_names[1])
    e2 <- slideCell::wf(sl_pair, sens = sens_matrix[[2]], alt = alt_matrix1, experiment=data_set_names[2])
    e3 <- slideCell::wf(sl_pair, sens = sens_matrix[[3]], alt = alt_matrix1, experiment=data_set_names[3])
    e4 <- slideCell::wf(sl_pair, sens = sens_matrix[[4]], alt = alt_matrix1, experiment=data_set_names[4])
    gridExtra::grid.arrange(e1, e2, e3, e4, ncol=2)


}
