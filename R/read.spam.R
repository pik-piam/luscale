#' Read spam files
#' 
#' Read transformation matrix from file. Data has the following structure:
#' \enumerate{ \item one 4byte Integer which contains the number of values
#' saved (n) \item n*2 matrix (integer) which contains the indices of the
#' spam-object (extracted with triplet(spam\_object)) \item vector with n
#' elements (integer) containing the values of the spam-object \item nrows,
#' ncols (integer) }
#' 
#' 
#' @usage read.spam(file_name, file_type=NULL)
#' @param file_name name of the spam-file
#' @param file_type File type, either "spam" or "sz" (for "spam-zipped"). File
#' ending is used if file_type is not mentioned
#' @return \item{input_spam}{spam-object containing the transformation matrix
#' in spam-format}
#' @note Wildcards are supported in \code{file_name}
#' @author Jan Philipp Dietrich
#' @export
#' @importFrom utils tail
#' @seealso \code{\link{write.spam}}
#' @examples
#' 
#' \dontrun{
#' require(spam)
#' a <- matrix(c(0,1,1,0,0,1),3,2)
#' b <- as.spam(a)
#' 
#' print(b)
#' #      [,1] [,2]
#' # [1,]    0    0
#' # [2,]    1    0
#' # [3,]    1    1
#' # Class 'spam'
#' 
#' write.spam(b,"test.spam")
#' 
#' c <- read.spam("test.spam")
#' 
#' print(c)
#' #      [,1] [,2]
#' # [1,]    0    0
#' # [2,]    1    0
#' # [3,]    1    1
#' # Class 'spam'
#' 
#' }
#spam_in_out provides functions to read and write spam objects from/into files.
#Version 1.10 - Jan Philipp Dietrich
# 1.01: Added ability to save matrices with numeric instead of integer values
# 1.02: removed compression for compatibility with python
# 1.03: fixed bug that caused huge file sizes because 0 were also saved
# 1.04: added "require(spam)"
# 1.05: inserted return-commands
# 1.06: read.spam/write.spam can now handle wildcards
# 1.07: Error is now raised if file in read.spam does not exist
# 1.08: included unexpanded file_name in warning message for ambiguous file names
# 1.09: matrix dimensions are now also saved
# 1.10: added compressed format "spam-zipped" .sz

#Write transformation matrix into file. Data has the following structure:
#  1. one 4byte Integer which contains the number of values saved (n)
#  2. n*2 matrix (integer) which contains the indices of the spam-object (extracted with triplet(spam_object))
#  3. vector with n elements (integer) containing the values of the spam-object 
#  4. (nrow,ncol)-vector (integer)

#Read transformation matrix.

read.spam <- function(file_name,file_type=NULL) {  
  if(length(Sys.glob(file_name))==0) {
	  stop(paste("file",file_name,"does not exist"))
	}

  #expand wildcards
  file_name_unexpanded <- file_name	
  file_name <- Sys.glob(file_name)
  if(length(file_name)>1) {
    file_name <- file_name[1]
    warning(paste("file name",file_name_unexpanded,"is ambiguous, only first alternative is used!"))
  }
  if(is.null(file_type)) {
    file_type <- tail(strsplit(file_name,'\\.')[[1]],1)
  }
  if(file_type == "gz" | file_type == "sz") {
    zz <- gzfile(file_name,"rb")
  } else {  
    zz <- file(file_name,"rb")
  }
  n_input <- readBin(zz,integer(),size=4,n=1)
  input_triplet <- list()
  input_triplet$indices <- matrix(readBin(zz,integer(),n=n_input*2),n_input,2)
  input_triplet$values <- readBin(zz,numeric(),n=n_input)
  out <- spam::spam(input_triplet)
  tmp <- readBin(zz,integer(),n=2) 
  if(length(tmp)==2) dim(out) <- tmp
  close(zz)
  return(out)
}