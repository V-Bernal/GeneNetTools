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


#' ci_pcor_shrunk
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
ci_pcor_shrunk <- function(lparams){
  # read file
  dt <- read_data(lparams$filename,lparams$variables)

  estimated.pcor.ecoli <- GGM.shrunk(x = as.matrix(dt),
                       verbose = lparams$verbose)

  # GGM <- vec2sm(
  #   result
  # )

  # p-values
  # pval.tts.ecoli <- ttest.shrunk(x = result)
  #
  # # GGM[abs(GGM) < lparams$cutoff] <- 0
  # # print(GGM)
  #
  # GGM <- vec2sm(
  #   pval.tts.ecoli
  # )
  #
  # GGM <- abs(GGM) < lparams$cutoff
  # diag(GGM) <- 0
  #
  # g1 <-
  #   igraph::graph_from_adjacency_matrix(
  #     abs(GGM),
  #     mode = c("undirected"),
  #     diag = FALSE,
  #     weighted = TRUE
  #   )
  #
  # plot(g1,
  #      edge.width = 3,
  #      edge.color = "#fc8d62",
  #      vertex.color = "#8da0cb", vertex.size = 3,
  #      vertex.frame.color = "#8da0cb",
  #      vertex.label.color = "black",
  #      vertex.label.cex = 1,
  #      vertex.label.dist = 1,
  #      vertex.label = colnames(dt),
  #      layout = igraph::layout
  # )

  #--------------------
  # Confidence intervals
  ci.ecoli <- confint.GGM(x = estimated.pcor.ecoli , alpha = 0.05)
  ci.ecoli

  ci.ecoli <- ci.ecoli[order(ci.ecoli[,1], decreasing = T),]
  ci.ecoli

  # CI forest plot
  ci.ecoli.df = data.frame(ci.ecoli[ci.ecoli[,2]>0, ])

  # un-comment for the complete list
  ci.ecoli.df = ci.ecoli.df[1:20,]

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
