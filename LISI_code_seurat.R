LSI <- function(obj,data_list,var.names,reduct.dim,reduct.dim2,center=FALSE,method = "irlba"){
  ##Decomposite the reference
  
  ##scaling
  if(!center) dat <- apply(data_list[[1]]@assays$RNA@data[var.names,],2,scale)
  if(center) dat <- apply(data_list[[1]]@assay$RNA@data[var.names,],2,scale,center=TRUE,scale=FALSE)
  if(method == "rARPACK") ref.svd <- rARPACK::svds(dat,k=reduct.dim)
  if(method == "irlba") ref.svd <- irlba::irlba(dat, nv = reduct.dim,center = TRUE)
  
  ##define encoder function
  encode.m <- diag(1/ref.svd$d) %*% t(ref.svd$u)
  
  ##perform scaling to all matrix, calculate cor after projection and concatenate together
  cor <- NULL
  for(i in 1:length(data_list)){
    if(!center) scale.data <- apply(data_list[[i]]@assay$RNA@data[var.names,],2,scale)
	if(center){scale.data <- apply(data_list[[i]]@assay$RNA@data[var.names,],2,scale,center=TRUE,scale=FALSE);scale.data <- scale.data[var.names,]}
	proj <- encode.m %*% scale.data
	cor <- cbind(cor, ref.svd$v%*%proj)
  }
  
  ##perform SVD and concatenate vector
  if(method == "rARPACK") embed <- rARPACK::svds(cor,k=reduct.dim2)$v
  if(method == "irlba") embed <- irlba::irlba(cor, nv =reduct.dim2)$v
  
  ###integrate
  rownames(embed) <- colnames(obj)
  colnames(embed) <- paste0("PC_",1:ncol(embed))
  obj[['lsi_reduct']] <- CreateDimReducObject(embed,key="PC_")
}

