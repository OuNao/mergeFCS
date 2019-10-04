normalize_flowFrames<-function(ff = list(), channels=c(1,2,3,7,10,11), verbose = TRUE){
  ff2<-ff
  for (marker in channels){
    row_max<-max(sapply(ff, nrow))
    col_mat<-matrix(data=NA_real_, nrow=row_max, ncol=length(ff))
    for (frame in 1:length(ff)) {
      col_mat[,frame]<-c(ff[[frame]]@exprs[,marker], rep(NA_real_, (row_max-nrow(ff[[frame]]))))
    }
    #norm_col_mat<-limma::normalizeQuantiles(col_mat)
    norm_col_mat<-preprocessCore::normalize.quantiles(col_mat)
    #mode(norm_col_mat)<-"integer"
    for (frame in 1:length(ff)) {
      ff2[[frame]]@exprs[,marker]<-norm_col_mat[!is.na(norm_col_mat[,frame]),frame]
    }
  }
  return(ff2)
}

#' @title A function to merge flow cytometry data files (.FCS).
#'
#' @description A function to merge flow cytometry data files (.FCS) and impute variable parameters using kNN on the common parameters.
#' @param files Vector of file names of the .FCS files to merge
#' @param common Vector of common columns number
#' @param variable Vector of variable columns number
#' @param normalize Quantile normalize the common parameters before imputation
#' @param verbose Verbosity
#' @return Merged flowFrame
#' @keywords FCS
#' @export
merge_FCS<-function(files = c(), common=c(1,2,3,7,10,11), variable=c(4,5,6,8,9), normalize = TRUE, verbose = TRUE) {

  if (!(is.vector(files) && is.vector(common) && is.vector(variable) && is.character(files) && is.numeric(common) && is.numeric(variable))) stop("Wrong parameter types.")

  num_frames<-length(files)

  for (fil in 1:num_frames) {
    if (!file.exists(files[fil])) stop(paste("File number", fil, "does not exist."))
  }

  # read files
  ff<-list()
  for (fil in 1:num_frames) {
    ff[[fil]]<-flowCore::read.FCS(files[fil])
  }

  # compensate
  for (fil in 1:num_frames) {
    ff[[fil]]<-flowCore::compensate(ff[[fil]], flowCore::keyword(ff[[fil]], "SPILL")$SPILL)
  }

  # normalize
  if (normalize) ff<-normalize_flowFrames(ff, channels = common, verbose = verbose)

  # unmerged compensated flowSet
  unmerged<-flowCore::flowSet(ff)

  # create merge matrix
  numero_cols<-length(common)+ num_frames*length(variable) + 1
  numero_linhas<-0
  for (i in 1:num_frames) numero_linhas<-numero_linhas+nrow(ff[[i]])
  nomes_linhas<-NULL
  nomes_colunas<-c()
  for (i in 1:length(common)) nomes_colunas<-c(nomes_colunas, gsub("NA ","",paste(parameters(unmerged[[num_frames]])@data$desc[common[i]],parameters(unmerged[[num_frames]])@data$name[common[i]])))
  for (frame in 1:num_frames) {
    for (i in 1:length(variable)) nomes_colunas<-c(nomes_colunas, gsub("NA ","",paste(parameters(unmerged[[frame]])@data$desc[variable[i]],parameters(unmerged[[frame]])@data$name[variable[i]])))
  }
  nomes_colunas[numero_cols]<-"tube"
  nomes<-list(NULL,nomes_colunas)
  if (verbose) {
    cat("\033[34;1mFCS files to merge:\033[0m", files, sep="\n")
    cat("\n\033[34;1mflowFrames number:\033[0m", num_frames, " / ncol:", numero_cols, " / nrow:", format(numero_linhas, scientific=FALSE), "\n")
    cat("\n\033[34;1mCommon parameters:\033[0m\n\t\t")
    for (i in 1:num_frames) cat("\033[32mflowFrame:\033[31m", i, "\t")
    cat("\n")
    for (i in 1:length(common)){
      cat("\033[32mcol number:\033[31m",common[i],"\033[0m\t")
      for (j in 1:num_frames) {
        cat(unmerged[[j]]@parameters@data$desc[common[i]], unmerged[[j]]@parameters@data$name[common[i]], "\t")
      }
      cat("\n")
    }
    cat("\n\033[34;1mVariable parameters:\033[0m\n\t\t")
    for (i in 1:num_frames) cat("\033[32mflowFrame:\033[31m", i, "\t")
    cat("\n")
    for (i in 1:length(variable)){
      cat("\033[32mcol number:\033[31m",variable[i],"\033[0m\t")
      for (j in 1:num_frames) {
        cat(unmerged[[j]]@parameters@data$desc[variable[i]], unmerged[[j]]@parameters@data$name[variable[i]], "\t")
      }
      cat("\n")
    }
    cat("\n\n\033[34;1mFinal parameter names:\033[0m\n")
    print(nomes_colunas)
  }
  if (verbose) cat("\n\033[33;1mCreating merged matrix...")
  matriz<-matrix(data = NA, nrow = numero_linhas, ncol = numero_cols, dimnames=nomes)

  # common markers
  if (verbose) cat("\n\n\033[33;1mSubstituting values of common parameters...")
  indice<-0
  coluna<-0
  for (col in 1:length(common)){
    val_comum<-c()
    for (i in 1:num_frames) {
      val_comum<-c(val_comum, exprs(unmerged[[i]])[,common[col]])
    }
    matriz[,col]<-val_comum
  }
  coluna<-coluna+length(common)

  # impute all values of variable markers
  if (verbose) cat("\n\n\033[33;1mImputing values of variable parameters...\033[0m")
  a<-Sys.time()
  coluna2<-coluna
  if (verbose) cat("\nProcessed column number:", coluna2, "/", (numero_cols-1))
  for (frame in 1:num_frames) {
    Z<-FNN::get.knnx(exprs(unmerged[[frame]])[,common], matriz[,1:length(common)], 1)
    for (col in 1:length(variable)){
      #prev<-FNN::knn.reg(exprs(unmerged[[frame]])[,common], matriz[,1:length(common)], exprs(unmerged[[frame]])[,variable[col]], 1)
      #pred<-prev$pred
      pred<-exprs(unmerged[[frame]])[,variable[col]][Z$nn.index]
      matriz[,coluna2+col]<-pred
      if (verbose) cat("\rProcessed column number:", (coluna2+col), "/", (numero_cols-1))
    }
    coluna2<-coluna2+length(variable)
  }
  time_impute<-Sys.time()-a
  if (verbose) cat("\n\033[0mTotal time to impute values of variable parameters:\033[31m", time_impute, attributes(time_impute)$units)

  # substitute imputed values with real when exist
  if (verbose) cat("\n\n\033[33;1mSubstituting values of real variable parameters...")
  for (frame in 1:num_frames) {
    for (col in 1:length(variable)){
      matriz[(indice+1):(indice+nrow(unmerged[[frame]])),coluna+col]<-exprs(unmerged[[frame]])[,variable[col]]
    }
    matriz[(indice+1):(indice+nrow(unmerged[[frame]])),numero_cols]<-frame
    indice<-indice+nrow(unmerged[[frame]])
    coluna<-coluna+length(variable)
  }

  # create new merged flowFrame
  if (verbose) cat("\n\n\033[33;1mCreating new flowFrame...")
  meta <- data.frame(name=nomes_colunas, desc=NA, range=262144L)
  meta$minRange<-apply(matriz, 2, min)
  meta$maxRange<-apply(matriz, 2, max)
  rownames(meta)<-sprintf("$P%s", 1:numero_cols)
  keys<-list()
  for (col in 1:length(nomes_colunas)) {
    keys[[paste0("$P",col,"N")]]<-as.character(nomes_colunas[col])
    if (col %in% 1:length(common)) {
      keys[[paste0("$P",col,"R")]]<-unmerged[[num_frames]]@description[[paste0("$P",common[col],"R")]]
    } else if (col<length(nomes_colunas)) {
      frame<-as.integer(((col-length(common)-1)/length(variable))+1)
      colunm<-(col-((frame-1)*length(variable)))-length(common)
      keys[[paste0("$P",col,"R")]]<-unmerged[[frame]]@description[[paste0("$P",variable[colunm],"R")]]
    } else {
      keys[[paste0("$P",col,"R")]]<-num_frames
    }
  }
  merged<-flowFrame(exprs = matriz, description = keys, parameters = Biobase::AnnotatedDataFrame(meta))

  if (verbose) cat("\n\n\033[32;1mDone...")
  return(merged)
}

#' @title A function to concatenate flow cytometry data files (.FCS).
#'
#' @description A function to concatenate flow cytometry data files (.FCS) with the same parameters.
#' @param fs FlowSet with .FCS files with identical parameters to train the model.
#' @param ev_max Max number of events from each flowFrame.
#' @param ev_min Max number of events from each flowFrame.
#' @param verbose Verbosity.
#' @return List with concatenated flowFrame, files, samples, ev_max, ev_min and verbose values if verbose = TRUE or only the concatenated flowFrame if verbose = FALSE.
#' @keywords FCS
#' @export
concat_FCS<-function(fs, ev_max=5000, ev_min=0, verbose = TRUE) {
  if (class(fs)[1] != "flowSet") stop("fs is not a flowSet!")
  if (!(is.numeric(ev_max)) || !(is.numeric(ev_min))) stop("ev_max and ev_min must be numeric!")
  is.equal <- function(mylist) {
    check.eq <- sapply(mylist[-1], function(x) {x == mylist[[1]]})
    as.logical(apply(check.eq, 1, prod))
  }
  param_names<-lapply(1:length(fs), function(x) {names(fs[[x]])})
  if (sum(!is.equal(param_names)) > 0) stop("The parameters of all flowFrames in fs must be identical")
  num_frames<-length(fs)
  samples<-list()
  for (i in 1:num_frames) {
    if (nrow(fs[[i]])>ev_max) samples[[i]]<-sample(1:nrow(fs[[i]]), ev_max, replace=FALSE) else if (nrow(fs[[i]])<ev_min) samples[[i]]<-0 else samples[[i]]<-sample(1:nrow(fs[[i]]), nrow(fs[[i]]), replace=FALSE)
  }
  ff<-fs[[1]][samples[[1]]]
  for (i in 2:num_frames) ff@exprs<-rbind(ff@exprs, fs[[i]][samples[[i]]]@exprs)
  ff@description$GUID<-"Concatenated flowFrame"
  files<-ls(fs@frames)
  return<-list("frame"=ff, "files"=files, "samples"=samples, "ev_max"=ev_max, "ev_min"=ev_min, "verbose"=verbose)
  if (verbose) return(return) else return(ff)
}

### todo
# translate all to english
