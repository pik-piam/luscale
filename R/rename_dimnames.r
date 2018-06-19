#' rename_dimnames
#' 
#' Renames the dimnames of an array or MAgPIE object after a query.
#' 
#' 
#' @usage rename_dimnames(data,dim=1,query=NULL,from=NULL,to=NULL)
#' @param data Array
#' @param dim The dimension to be renamed.
#' @param query If NULL, query is automatically searched for. Otherwhise an
#' array, data.frame or the path of a csv with at least two columns. One column
#' has to have the name of "from", the other one the name of "to". Some queries
#' can be found in the svn-folder tools/queries.
#' @param from Only required if query is not NULL. Column of the query with
#' original dimnames of the incoming dataset
#' @param to Only required if query is not NULL. Column of the query with the
#' target dimnames of the outcoming dataset
#' @return An array with different dimnames
#' @note translate.with.query has the same functionality, is more efficient,
#' yet more complicated to use.
#' @author Benjamin Bodirsky, Ulrich Kreidenweis
#' @export
#' @importFrom magclass is.magpie unwrap as.magpie
#' @keywords array
#################################
#### rename by query         ####
#################################
                                  
# Version 1.00 - Benjamin Bodirsky
# Version 1.01: minor change: also allow data.frames as query - Ulrich Kreidenweis

rename_dimnames<-function(data,dim=1,query=NULL,from=NULL,to=NULL) {
  if(dim-floor(dim)>0){ dim=old_dim_convention(dim) }  # check whether old naming convention is active
  if(is.data.frame(query)) {
    query <- sapply(query,as.character)
  }
  
  ismagpie<-FALSE
  if(!is.array(data)){stop("data has to be an array")}  
  if(dim>2){
    if (is.magpie(data)){
      data<-unwrap(data)
      ismagpie<-TRUE
    }
  }
  
  tmp=FALSE
  if (!is.null(to)) { 
    if((dim==1)&(to%in%(c("glo","GLO")))){
      tmp=TRUE
    }
  }
    
  if (tmp==TRUE) {
    dimnames(data)[[dim]]<-paste("GLO",1:dim(data)[dim],sep=".")
  } else {
    #select the adequate mapping
    
    dimnames(data)[[dim]]<-iconv(dimnames(data)[[dim]],from="ANSI_X3.4-1986",to="ASCII")
    data_categories<-as.vector(dimnames(data)[[dim]])
    
    mapping <- AutomaticMapping(x=data_categories,mapping=query,from=from,to=to)
    query<-mapping
    from<-colnames(mapping)[1]
    to<-colnames(mapping)[2]
    
    #start actual function  
    query_categories_from <- as.vector(query[,from])
    query_categories_to <- as.vector(query[,to])
    names(query_categories_to) <- query_categories_from
    if (any ((data_categories %in% query_categories_from)==FALSE)) {
      warning("Following categories in data had no query entry: ",paste(unique(data_categories[which((data_categories %in% query_categories_from)==FALSE)]),collapse = ' '))
    }    
    result <- query_categories_to[data_categories]
    if (any(is.na(result))){result[which(is.na(result))]<-"XXX"}
    dimnames(data)[[dim]]<-result
  
  }   
  
  
  if (ismagpie==TRUE){
    data<-as.magpie(as.array(data))
  }
  
  
  return(data)
}
