#' Reshape Folder
#' 
#' Disaggregates all cellular MAgPIE output in a folder to 0.5 degree files
#' based on spam matrices.
#' 
#' 
#' @usage reshape_folder(folder=".")
#' @param folder folder in which the files should be disaggregated
#' @author Jan Philipp Dietrich
#' @export
#' @importFrom magclass write.magpie
#' @seealso \code{\link{reshape_file}}, \code{\link{speed_aggregate}}
#' @examples
#' 
#'  \dontrun{reshape_folder()}
#' 
reshape_folder <- function(folder=".") {
  files <- grep("^cell\\.",dir(folder),value=TRUE)
  for(f in files) {
    out <- reshape_file(path(folder,f))
    if(!is.null(out)) write.magpie(out,path(folder,gsub("\\.[^\\.]+$","_0.5.mz",f)))
  }
}
