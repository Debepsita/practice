

lrd=function(matrices,lower=30,upper=100)
{
  library(r.jive)
  library(PMA)
  # Check if the input is a list of matrices
  if (!all(sapply(matrices, is.matrix))) {
    stop("All inputs must be matrices.")
  }
  #function to give the updated matrices containing only the common columns
  common_columns <- function(matrices) {


    # Find the common column names
    com_cols <- Reduce(intersect, lapply(matrices, colnames))

    # If no common columns are found
    if (length(com_cols) == 0) {
      message("No common columns found.")
      return(NULL)
    }

    # Subset matrices to include only the common columns
    updated_matrices <- lapply(matrices, function(mat) {
      mat[, com_cols, drop = FALSE]
    })

    return(updated_matrices)
  }
  comd=common_columns(matrices)
  if (is.null(comd)) {
    stop("No common columns found.")
  }
  rr=jive(comd,maxiter=20,showProgress="False")
  jm=rr$joint
  h2=do.call(rbind, jm) #rbind the list of matrices.
  h_sub2=h2
  pro.list <- vector("list", length = rr$rankJ)
  n1=colnames(comd[[1]])
  c1=lower
  c2=upper
  for(i in 1:(rr$rankJ))
  {
    pp=length(n1)
    l=104
    z=1

    while (l<c1 || l>c2 && z< pp^(0.25)) {



      tt1=SPC(h_sub2,sumabsv=z, K=1)
      l=length(tt1$v[tt1$v!=0])
      z=z+0.1


    }

    ind=which(tt1$v!=0) # indices of protein previously selected

    pro.list[[i]]=n1[ind]
    x=(tt1$d)*(tt1$u)%*%t((tt1$v))
    hs1=h_sub2-x
    h_sub2=hs1[,-ind]
    n1=n1[-ind]

  }
  return(pro.list)
}






