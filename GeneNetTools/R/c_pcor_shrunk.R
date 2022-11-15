#' Validate that a  json file meets a required schema
#'
#' @param params a json file to be validated against a schema
#' @param pschema a predefined json schema
#'
#' @export
validate_parameters <- function(params,pschema="pca_projection_schema.json"){
  schemafile <- system.file("extdata", pschema, package = "GeneNetTools")
  jsonvalidate::json_validate(params,schemafile,verbose=TRUE,error=TRUE)
}

#' Reads columns from a file in table format
#'
#' Additional info
#' @param filename a string data file name including the relative path
#' @param select_columns a vector including the column names to be read from the data file
#'
#' @return if succeed, this function returns a data table
#' @export
#' @keywords internal
#'
# #' @examples
read_data <- function(filename,select_columns){
  # Check for empty list
  if (length(select_columns)<1) {
    select_columns = NULL
    cat("Reading all columns\n")
  } else {
    cat("Selected columns",select_columns,"\n")
  }
  tryCatch(cols <- data.table::fread(filename,select = select_columns),
           error = function(c) {
             c$message <- paste0(c$message, " (in ", filename, ")")
             stop(c)
           }
           # ,warning = function(c) {
           #   c
           #   }
  )
  # print(cols)
}

#' This function validates a json structure
#'
#' @param fileparams a json structure stored in a file name to be validated
#'
#' @return an R list of parameters extracted from the json structure
#' @export
#'
# #' @examples
validate_json_file <- function(fileparams) {
  if (file.exists(fileparams)){
    tryCatch(lp <-  jsonlite::fromJSON(fileparams),
             error = function(c) {
               c$message <- paste0(c$message, " (in ", fileparams, ")")
               stop(c)
             }
    )
  } else {
    message <- paste0("Parameter file '", fileparams, "' not found.")
    stop(message)
  }
}


#' Partial correlation shrunk
#'
#' This function computes confidence intervals for the partial correlation with shrinkage.
#' @param lparams a list of parameters created using a JSON file.
#' This file should contain the following name/value pairs.
#'
#' "filename": <string, required>
#'
#' "variables": <array, strings representing column names>
#'
#' "cutoff": <number, required threshold for the p-value of the partial correlation>
#'
#' "verbose": <boolean, required to display detailed description on the terminal>
#'
#' @return Forest plot of partial correlations in Rplot.pdf
#' @export
#'
# #' @examples
c_pcor_shrunk <- function(lparams){
  # read file

  # read file
  dt <- read_data(lparams$filename,lparams$variables)

  # estimate partial correlations with shrinkage
  estimated.pcor.ecoli <- GGM.shrunk(x = as.matrix(dt),
                                     verbose = lparams$verbose)

  # Confidence intervals
  ci.ecoli <- confint.GGM(x = estimated.pcor.ecoli , alpha = 0.05)
  ci.ecoli <- ci.ecoli[order(ci.ecoli[,1], decreasing = T),]

  # CI forest plot
  ci.ecoli.df <- data.frame(ci.ecoli[ci.ecoli[,2]>0, ])

  # un-comment for the complete list
  ci.ecoli.df <- ci.ecoli.df[1:20,]

  p <- ggplot2::ggplot(data=ci.ecoli.df, ggplot2::aes(y= 1:nrow(ci.ecoli.df),
                                                      x=ci.ecoli.df[,1], xmin=ci.ecoli.df[,2],
                                                      xmax=ci.ecoli.df[,3])) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(height=.4) +
    ggplot2::labs(title='', x='scaled pcors (95% CI)', y = '') +
    ggplot2::geom_vline(xintercept=c(0,.1,.3), color='black', linetype='dashed', alpha=.5) +
    ggplot2::scale_y_continuous(name = "", breaks=1:nrow(ci.ecoli.df), labels=row.names(ci.ecoli.df))+
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = "transparent", color = "black"),
                   panel.background = ggplot2::element_rect(fill = "white"))

  print(p)

}

#==============================
#' pval_pcor_shrunk
#'
#' TODO: add proper title and description
#'
#' @param lparams a list of parameters created using a JSON file.
#' This file should contain the following name/value pairs.
#'
#' "filename": <string, required>
#'
#' "variables": <array, strings representing column names>
#'
#' "cutoff": <number, required threshold for the p-value of the partial correlation>
#'
#' "verbose": <boolean, required to display detailed description on the terminal>
#'
#' @return a plot
#' @export
#'
# #' @examples
c_pval_pcor_shrunk <- function(lparams){

  # read file
  dt <- read_data(lparams$filename,lparams$variables)

  # estimate partial correlations with shrinkage
  estimated.pcor <- GGM.shrunk(x = as.matrix(dt),
                                     verbose = lparams$verbose)

  # p-values using the new implementation Equation 6
  pval.ttest <- ttest.shrunk(x = estimated.pcor)

  GGM <- vec2sm(pval.ttest)
  GGM <- abs(GGM) < lparams$cutoff
  diag(GGM) <- 0


  g1 <-
    igraph::graph_from_adjacency_matrix(
      abs(GGM),
      mode = c("undirected"),
      diag = FALSE,
      weighted = TRUE
    )

  plot(g1,
       edge.width = 3,
       edge.color = "#fc8d62",
       vertex.color = "#8da0cb", vertex.size = 3,
       vertex.frame.color = "#8da0cb",
       vertex.label.color = "black",
       vertex.label.cex = 1,
       vertex.label.dist = 1,
       vertex.label = colnames(dt)# ,
       #layout = igraph::layout
  )

  #print(p)
}

#==============================
#' c_zscore_shrunk
#'
#' This function compares two networks using two data files
#'
#' TODO: add proper title and description
#'
#' @param lparams a list of parameters created using a JSON file.
#' This file should contain the following name/value pairs.
#'
#' "filename": <string, required data file 1>
#'
#' "filename2": <string, required data file 2>
#'
#' "variables": <array, strings representing column names>
#'
#' "cutoff": <number, required threshold for the p-value of the partial correlation>
#'
#' "verbose": <boolean, required to display detailed description on the terminal>
#'
#' @return a plot
#' @export
# #' @examples
c_zscore_shrunk <- function(lparams){
  # se definen dos archivos de entrada en el JSON schema y se agregan en el JSON file
  # read file
  data1 <- read_data(lparams$filename,lparams$variables)
  data2 <- GeneNet::ggm.simulate.data(50, diag(ncol(data1)))
  #read_data(lparams$filename2,lparams$variables)

  estimated.pcor.1 <- GGM.shrunk(as.matrix(data1))
  estimated.pcor.2 <- GGM.shrunk(as.matrix(data2))

  # compare
  z <- compare.GGM(x1 = estimated.pcor.1 ,
                   x2 = estimated.pcor.2  )
  cut = 1.96
  id = ( z[ , 'z-score'] > cut) + 2*( z[ , 'z-score'] < -cut)

  #---------------------
  # Bland Altman plot
  #---------------------
  # Matrix index of the top pcors
  id0 = order(abs(z[ , 'z-score']), decreasing = T)[1:9]

  idxs = matrix(data = 1:(ncol(data1)*ncol(data1)), nrow = ncol(data1),
                ncol = ncol(data1))
  idxs2 = matrix(data = idxs %in% id0,
                 nrow = ncol(data1), ncol = ncol(data1))
  idxs3 = which(idxs2, arr.ind = T)

  # 9 top names
  #res = data.frame( 'probe1'= all_new_gene$external_gene_name[idxs3[,1]],
  #                  'probe2' = all_new_gene$external_gene_name[idxs3[,2]] )
  # plot
  df = data.frame('diff' = (estimated.pcor.1 - estimated.pcor.2),
                  'ave' = 0.5*(estimated.pcor.1 + estimated.pcor.2))
  colnames(df) <- c('diff' , 'ave' )

  plot( x = df$ave, y = df$diff,
        pch = 20, cex= 0.5 ,
        cex.lab = 1,  cex.axis = 1,
        ylim = c(-0.1,0.1),
        col = c('black'),
        xlab = 'ave pcors B6 D2',
        ylab = 'diff pcors B6 D2', las = 2)
  abline(a = 0, b = 0)
  #text(x = df$ave[id0] ,
  #     y = df$diff[id0]*1.1 ,
  #     labels = paste(res[,1],res[,2],sep = '-'),
  #     cex = 0.75, adj = 1)
}
