#' @title GetData_CITMIC
#' @description Get the example data
#' @param Data  A character should be one of "GEP", "matrix_cell_go_inter", "matrix_cell_go_jaccard"
#' @return Data
#' @export

GetData_CITMIC<-function(Data){
  if(!exists("CITMIC_Data")) {
    utils::data("CITMIC_Data",package="CITMIC")
  }
  if (Data=="GEP")
  {
    dataset<- get("GEP",envir=CITMIC_Data)
    return(dataset)
  }
  if (Data=="matrix_cell_go_inter")
  {
    dataset<- get("matrix_cell_go_inter",envir=CITMIC_Data)
    return(dataset)
  }
  if (Data=="matrix_cell_go_jaccard")
  {
    dataset<- get("matrix_cell_go_jaccard",envir=CITMIC_Data)
    return(dataset)
  }
}
