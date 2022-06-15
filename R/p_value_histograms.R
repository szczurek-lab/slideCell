#' plots p-value histogram for a list of p-values
#' enrichment set - geneA altered cell lines
#' geneB sensitivity is a ranking criterion
#' @param p_values p-values vector
#' @param name  name of the data set used to produce p-values
#' @param b_param the threshold parameter - the probablity that the number of p-values in each bin will exceed the plotted level, default is 1%
#' @return p-values histogram
#' @export

p_value_one_histogram <- function(p_values, name, b_param=0.01) {
    #plot p-values histogram with bin width equal to b_param
    hist(p_values, main = NA, breaks = seq(0, 1, b_param), xlab="p-values")
    title(paste("Histogram for", name), line = 0,  cex.main = 0.8, )
    abline(h=qbinom((1-b_param), length(p_values), b_param), col="red")
}


#' plots 4 p-value histograms
#' @param p_values_list list of p-values vectors from 4 different testing setups
#' @param data_set_names names for each testing setup
#' @param b_param the threshold parameter - the probablity that the number of p-values in each bin will exceed the plotted level, default is 1%
#' @return p-values histograms
#' @export

plot_many_histograms <- function(p_values_list, data_set_names, analysis_type, b_param=0.01) {

    par(mfrow = c(2, 2))
    slideCell::p_value_one_histogram(na.omit(p_values_list    ), name = data_set_names[1], b_param=b_param)
    slideCell::p_value_one_histogram(na.omit(p_values_list[2][[1]]), name = data_set_names[2], b_param=b_param)
    slideCell::p_value_one_histogram(na.omit(p_values_list[3][[1]]), name = data_set_names[3], b_param=b_param)
    slideCell::p_value_one_histogram(na.omit(p_values_list[4][[1]]), name = data_set_names[4], b_param=b_param)
    title(paste("p-values histograms for", analysis_type), line = 0, outer = TRUE)
    par(mfrow = c(1, 1))

}


