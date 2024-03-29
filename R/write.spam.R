#' Write file in spam format
#'
#' Write transformation matrix into file. Data has the following structure:
#' \enumerate{ \item one 4byte Integer which contains the number of values
#' saved (n) \item n*2 matrix (integer) which contains the indices of the
#' spam-object (extracted with \code{triplet(spam_object)}) \item vector with n
#' elements (integer) containing the values of the spam-object \item nrows,
#' ncols (integer) }
#'
#'
#' @usage write.spam(input_spam, file_name, file_type=NULL)
#' @param input_spam Spam-object, that should be written to file
#' @param file_name file name the object should be written to
#' @param file_type File type, either "spam" or "sz" (for "spam-zipped"). File
#' ending is used if file_type is not mentioned
#' @note Wildcards are supported in \code{file_name}
#' @author Jan Philipp Dietrich
#' @export
#' @importFrom utils tail
#' @seealso \code{\link{read.spam}}
#' @examples
#'
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
#' file.remove("test.spam")
#'
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



write.spam <- function(input_spam,file_name,file_type=NULL) {
  #expand wildcards
  file_name <- paste(Sys.glob(dirname(file_name)),basename(file_name),sep="/")
  if(length(file_name)>1) {
    file_name <- file_name[1]
    warning("file name is ambiguous, only first alternative is used!")
  }
  input_spam <- spam::as.spam(input_spam)

  #if(format(Sys.time(),"%d%m") == "0104") schneeball()
  #if file-type is not mentioned file-ending is used as file-type
  if(is.null(file_type)) {
    file_type <- tail(strsplit(file_name,'\\.')[[1]],1)
  }
  if(file_type == "gz" | file_type == "sz") {
    zz <- gzfile(file_name,"wb")
  } else {
    zz <- file(file_name,"wb")
  }
  writeBin(as.integer(length(spam::triplet(input_spam)$values)),zz,size=4)
  writeBin(as.integer(as.vector(spam::triplet(input_spam)$indices)),zz)
  writeBin(as.numeric(as.vector(spam::triplet(input_spam)$values)),zz)
  writeBin(as.integer(dim(input_spam)),zz)
  close(zz)
}
