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
#' @return if succeds, this function returns a data table
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


#' Title
#'
#' @param lparams
#'
#' @return
#' @export
#'
#' @examples
pcor_shrunk_api <- function(lparams){
  # read file

  dt <- read_data(lparams$filename,lparams$variables)

  result <- GGM.shrunk(x = as.matrix(dt), lambda = lparams$lambda, verbose = lparams$verbose)
  # pl <- hist(result)
  #print(pl)

  print(class(result))

  result <- abs(result)>0.3

  GGM <- GeneNetTools::vec2sm.tools(
    result
  )

  diag(GGM) <- 0

  g1 <-
    igraph::graph_from_adjacency_matrix(
      abs(GGM),
      mode = c("undirected"),
      diag = FALSE,
      weighted = TRUE
    )

  plot(g1 #,
#       layout = layout_nicely(g1)
  )

}
