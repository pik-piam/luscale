#' AutomaticMapping
#' 
#' Function automatically finds fitting mapping for provided MAgPIE object for a given target aggregation.
#' 
#' 
#' @usage AutomaticMapping(x,mapping=NULL,from=NULL,to=NULL)
#' @param x MAgPIE object
#' @param mapping A array or data.frame containing a mapping query or mapping name
#' @param from Only required if query is not NULL. Column of the query with
#' original dimnames of the incoming dataset
#' @param to Only required if query is not NULL. Column of the query with the
#' target dimnames of the outcoming dataset
#' @author Benjamin Leon Bodirsky
#' @export
#' @importFrom utils data read.csv
#' @export
AutomaticMapping<-function(x,mapping=NULL,from=NULL,to=NULL){

  if (is.null(to)){to="NULL"}
  #if (is.null(from)){from="NULL"}
  
  luqueries <- NULL
  data("luqueries", envir=environment() ,package = "luscale")
  if(!is.vector(x)){stop("x has to be a vector")}
  x<-iconv(x,from="ANSI_X3.4-1986",to="ASCII")
  
  out<-matrix(data = c(x,x),ncol = 2,nrow = length(x),dimnames = list(NULL,c("from","to")))
  out[,2]<-"XXX"
  
  if (to=="glo"|to=="GLO") { 
    out[,2]<-"GLO"
  } else {
    
    
    #select mapping if none specified:
    if (is.null(mapping)) {
      # use standard mapping
        if (all(x %in% luqueries$spatial$fao_iso[,"fao"])) {      
          mapping<-"fao_iso"
          message("using mapping fao_iso")
        } else if (all(x %in% luqueries$spatial$fbs_iso[,"fbs"])) {      
          mapping<-"fbs_iso" 
          message("using mapping fbs_iso")   
        } else if (all(substr(x,1,3) %in% luqueries$spatial$iso_reg[,"iso"])) {          
          mapping<-"iso_reg"       
          message("using mapping iso_reg")
        } else if (all(substr(x,1,3) %in% luqueries$spatial$wb_iso[,"wb"])) {      
          mapping<-"wb_iso"         
          message("using mapping wb_iso")               
        } else if (all(x %in% luqueries$spatial$cell_iso[,"cell"])&(to=="iso")) {          
          mapping<-"cell_iso"       
          message("using mapping cell_iso")            
        } else if (all((x %in% luqueries$goods$faostat_kcr[,"fao"]))) {
          mapping<-"faostat_kcr"
          message("using mapping faostat_kcr")
        } else if (all(x %in% luqueries$goods$faostat_kli[,"fao"])) {
          mapping<-"faostat_kli"
          message("using mapping faostat_kli")            
        } else {stop("No standard mapping exists for these names. Please specify mapping.")}    
    } 
    
    #section to define standard queries  
    if (is.array(mapping) | is.data.frame(mapping)) {
      mapping<-mapping
      cn <- colnames(mapping)
      if(is.null(from) & to=="NULL"){
        from <- cn[1]
        to <- cn[2]
      } else if(is.null(from)) {
        from <- setdiff(cn,to)
      } else if(to=="NULL") {
        to <- setdiff(cn,from)
      }
    } else if (mapping=="fao_iso") {
      mapping <- luqueries$spatial$fao_iso
      if (to=="reg"){again=TRUE}
      if (is.null(from)){from<-"fao"}
      if (to=="NULL"){to<-"iso"}
    } else if (mapping=="fbs_iso") {
      mapping <- luqueries$spatial$fbs_iso  
      if (to=="reg"){again=TRUE}
      if (is.null(from)){from<-"fbs"}          
      if (to=="NULL"){to<-"iso"}   
    } else if (mapping=="iso_reg") {
      x<-substr(x,1,3)
      mapping <- luqueries$spatial$iso_reg                
      if (is.null(from)){from<-"iso"}         
      if (to=="NULL"){to<-"reg"}   
    } else if (mapping=="wb_iso") {  
      x<-substr(x,1,3)
      mapping <- luqueries$spatial$wb_iso  
      if (is.null(from)){from<-"wb"}
      if (!(to%in%c("NULL","iso"))){stop("wb_iso incompatible with selection of to")}
      if (to=="NULL"){to<-"iso"}  
    } else if (mapping=="cell_iso") {  
      mapping <- luqueries$spatial$cell_iso         
      if (is.null(from)){from<-"cell"}
      if (to=="NULL"){to<-"iso"}             
    } else if (mapping%in%c("cell_reg","cluster_reg")) {  
      mapping <- array(
        dim=c(length(x),2),
        dimnames=list(NULL,c("cell","reg")),
        data=c(x,
               substr(x,1,3))
        )
      if (is.null(from)){from<-"cell"}
      if (to=="NULL"){to<-"reg"}                                    
    } else if (mapping=="faostat_kcr") {
      mapping<-luqueries$goods$faostat_kcr
      if (is.null(from)){from<-"fao"}
      if (to=="NULL"){to<-"MAgPIE_item"}  
    } else if (mapping=="faostat_kli") {
      mapping<-luqueries$goods$faostat_kli
      if (is.null(from)){from<-"fao"}
      if (to=="NULL"){to<-"MAgPIE_item"} 
    } else if (is.character(mapping)) {
      if(file.exists(mapping)){
        mapping<-read.csv(mapping,header=TRUE,sep=";")  
      } else {
        stop("no mapping or mapping file with the name indicated in mapping")
      }
      mapping<-array(dim=c(dim(mapping)[1],2),data=c(as.vector(mapping[,from]),as.vector(mapping[,to])),dimnames=list(as.vector(mapping[,to]),c(from,to)))            
    } else {stop("mapping has to be an array, data.frame or the path to an csv file")}      
    
    if (is.null(from)|(to=="NULL")){stop("from and to needed")}
    if ((from %in% dimnames(mapping)[[2]])==FALSE) {stop("from hast to be a colname of the mapping")}
    if ((to %in% dimnames(mapping)[[2]])==FALSE) {stop("to has to be a colname of the mapping")}     
    
    mapping<-mapping[,c(which(dimnames(mapping)[[2]]==from),which(dimnames(mapping)[[2]]==to))]
    
    #which outs exist in mapping
    matching<-which(out[,1]%in%mapping[,1])
    #which dont exist in mapping
    unmatch<-which(!out[,1]%in%mapping[,1])
    if (length(unmatch)>0) {
      warning("Following categories in data had no mapping entry: ",paste(unique(out[unmatch,1]),". mapping set to XXX",collapse = ' '))
    }
    #replace outs that exist by mapping entry
    out[matching,2]<-mapping[match(out[matching,1],mapping[,1]),2]
  
  }
  colnames(out)=c(from,to)
  return(out)
}