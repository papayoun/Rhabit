checkList <- function(my_list, name = "my_list"){
  if (!inherits(my_list, "list")){
    stop(paste(name, "should be a list"))
  }
  else{
    test_null <- sapply(my_list, is.null)
    if (any(test_null))
      stop(paste("Element(s):", paste(which(test_null), collapse = " and "),
                 "has/have null element in both cov_list and grad_fun"))
  }
  return(NULL)
}

check2Lists <- function(list1, list2){
  if (!(inherits(list1, "list") & inherits(list2, "list"))){
    stop("cov_list and grad_fun should be either NULL or lists")
  }
  else{
    test_null <- mapply(function(x, y){
        is.null(x) & is.null(y)
      },
      list1, list2)
    if (any(test_null))
      stop(paste("Element(s):", paste(which(test_null), collapse = " and "),
                 "has/have null element in both cov_list and grad_fun"))
  }
  return(NULL)
}
checkCovGrad <- function(cov_list, grad_fun){
  if (is.null(cov_list) & is.null(grad_fun)){
    stop("Either cov_list of grad_fun must be non NULL")
  }
  else if (is.null(cov_list)){
    checkList(grad_fun, "grad_fun")
  }
  else if (is.null(grad_fun)){
    checkList(cov_list, "cov_list")
  }
  else{
    check2Lists(cov_list, grad_fun)
  }
}
