#'@title CITMIC
#' @description The function "CITMIC" is used to identify cell infiltration in tumor microenvironment by calculating intercellular crosstalk.
#'@param GEP An example gene expression profile.
#'@param weighted This parameter specifies whether to create a weighted graph for the cell crosstalk network. If null, an unweighted graph is created, and the elements of the adjacency matrix indicate the number of edges between vertices. If true, a weighted graph is created(default: TRUE).
#'@param base Standardized log base of data for improving data distribution(default: 10).
#'@param damping Restart the probability of the random-walk algorithm (default: 0.9).
#'@param cl.cores The number of CPU cores applied to this task(default:1).
#'@param cell.type Preset the relevant cell type (e.g. if the solid tumor tissue does not contain 'HSC', it is better to remove it when we preset it.). When cell.type = NULL, all 86 cell types from the evaluation project are used by default.
#'@return A data frame containing the cell infiltration score for each sample.
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph V
#' @importFrom igraph page_rank
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats pt
#' @importFrom stats sd
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom stats na.omit
#' @importFrom fastmatch fmatch
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterEvalQ
#'@usage CITMIC(GEP,weighted = TRUE,base = 10,damping=0.90,cl.cores=1,cell.type=NULL)
#'@export
#'@examples
#' # Obtain the example data
#' GEP<-GetData_CITMIC("GEP")
#' # Run the function
#' \donttest{lnScore<-CITMIC(GEP,weighted = TRUE,base = 10,damping=0.90,cl.cores=1,cell.type=NULL)}



CITMIC<-function(GEP,weighted = TRUE,base = 10,damping=0.90,cl.cores=1,cell.type=NULL){

  GEP<-GEP
  cl<-makeCluster(cl.cores)
  clusterExport(cl,varlist = "GEP",envir = environment())
  clusterEvalQ(cl, {
    matrix_cell_go_inter<- CITMIC::GetData_CITMIC("matrix_cell_go_inter")
    matrix_cell_go_jaccard<-CITMIC::GetData_CITMIC("matrix_cell_go_jaccard")
    go_cell_score_row<-function(a,table,del,c){
      score_row<-rep(0,c)

      for(j in 1:length(a)){
        if(a[j]!=""){
          gene<-unlist(strsplit(a[j], split = ","))
          location<-fastmatch::fmatch(gene, table)

          dell<- na.omit(del[location])
          de_score1<-median(as.numeric(dell))
          if (!is.na(de_score1)) {
            score_row[j]<-de_score1
          }

        }

      }
      return(score_row)
    }
    median_inter<-function(matrix_cell_go_inter,matrix_cell_go_jaccard,GEP){

      GEPscore<-cbind(rownames(GEP),GEP[,1])
      table <- GEPscore[,1]
      del <- GEPscore[, 2]
      median_score<-matrix(0,nrow=length(rownames(matrix_cell_go_inter)),ncol=length(colnames(matrix_cell_go_inter)))
      for(k in 1:length(rownames(matrix_cell_go_inter))){

        Genes_vector<-matrix_cell_go_inter[k,]
        row<-go_cell_score_row(Genes_vector,table,del,length(colnames(matrix_cell_go_inter)))
        median_score[k,]<-row

      }
      matrix_median_genes<-median_score*matrix_cell_go_jaccard
      colnames(matrix_median_genes)<-colnames(matrix_cell_go_inter)
      rownames(matrix_median_genes)<-rownames(matrix_cell_go_inter)
      matrix_cell_score<-t(matrix_median_genes)%*%matrix_median_genes
      matrix_cell_score[is.na(matrix_cell_score)]<-0
      diag(matrix_cell_score)<-0
      return(matrix_cell_score)
    }
  })
  result_cell<-parLapply(cl=cl,X=colnames(GEP),function(x){
                           Zvalue<-data.frame(cbind(rownames(GEP),GEP[,x]))
                           Zvalue <- data.frame(Zvalue[!is.infinite(Zvalue[,2]),])
                           Zvalue <- data.frame(Zvalue[!is.na(Zvalue[,2]),])
                           Zvalue <- data.frame(Zvalue[!is.nan(Zvalue[,2]),])
                           Zvalue[,2]<-as.numeric(Zvalue[,2])
                           Z<-data.frame(Zvalue[,2])
                           rownames(Z)<-Zvalue[,1]
                           score<-median_inter(matrix_cell_go_inter,matrix_cell_go_jaccard,Z)
                           return(score)
                         })
  names(result_cell)<-colnames(GEP)
  stopCluster(cl)
  random_crosstalk<-function(result_cell,damping=damping){



    adj.final<-as.matrix(result_cell)
    graph = graph_from_adjacency_matrix(adj.final,mode=c("undirected"),weighted=weighted)
    temp = page_rank(graph, vids=V(graph), directed=FALSE, damping=damping, weights=NULL)
    rank = temp$vector
    rank1 = as.matrix(rank)


    return(rank1)
  }
  Score_rankwalk<-data.frame(row.names=rownames(result_cell[[1]]))
  for(i in names(result_cell)){
    score<-random_crosstalk(result_cell[[i]],damping)
    Score_rankwalk<-cbind(
      Score_rankwalk,
      score[rownames(Score_rankwalk),]
    )

  }
  colnames(Score_rankwalk)<-names(result_cell)
  cell_num<-apply(Score_rankwalk,2,function(x){
    y=log(x,base=base)
    return(y)
  })
  lnscore<-(cell_num-min(cell_num))/(max(cell_num)-min(cell_num))
  if(is.null(cell.type)){
    lnscore<-lnscore
  }else{
    lnscore<-lnscore[cell.type,]
  }
  return(lnscore)
}

