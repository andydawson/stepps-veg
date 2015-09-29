
other_build <- function(y, counts, other.idx){
  
  sum.other = sum(other.idx, na.rm=TRUE)
  
  if (sum.other>1){
    if (!is.null(y$other.hardwood))
      other.vec = y$other.hardwood + rowSums(counts[, other.idx])
    else 
      other.vec = rowSums(counts[, other.idx])
  } else if  (sum.other==1) {
    if (!is.null(y$other.hardwood))
      other.vec = y$other.hardwood + counts[, other.idx]
    else 
      other.vec = counts[, other.idx]
  } else {
    if (!is.null(y$other.hardwood))
      other.vec = y$other.hardwood 
    else 
      other.vec = rep(0, nrow(y))
  }
  
  return(other.vec)
}

# split the data to separate michigan upper and lower peninsula
split_mi <- function(meta, longlat){
  
  if (any(colnames(meta)=='region')){
    meta$state = meta$region
  } 
  
  if (longlat){
    centers_ll = data.frame(x=meta$long, y=meta$lat)
  } else {
    centers = data.frame(x=meta$x, y=meta$y)
    
    coordinates(centers) <- ~x + y
    proj4string(centers) <- CRS('+init=epsg:3175')
    
    centers_ll <- spTransform(centers, CRS('+proj=longlat +ellps=WGS84'))
    centers_ll <- as.matrix(data.frame(centers_ll))
  }
  
  if (!any(colnames(meta) == 'state')){
    meta = data.frame(meta, state=rep(NA, nrow(meta)))
    meta$state = map.where(database='state', centers_ll[,1], centers_ll[,2])
    meta[((meta$state == 'michigan:south') & (!is.na(meta$state))), 'state'] = 'michigan_south'
    meta[((meta$state == 'michigan:north') & (!is.na(meta$state))), 'state'] = 'michigan_north'
  }
  
  idx.mi = which(meta$state=='michigan_north')
  meta$state2 = as.vector(meta$state)
  meta$state2[idx.mi] = map.where(database="state", centers_ll[idx.mi,1], centers_ll[idx.mi,2])
  
  if (any(is.na(meta$state2))){
    idx.na = which(is.na(meta$state2))
  }
  idx.not.na = which(!is.na(meta$state2))
  
  idx.mi.s = which(meta$state=='michigan_south')
  meta$state2[idx.mi.s] = 'michigan:south'
  
  if (length(idx.na)>0){
    for (i in 1:length(idx.na)){
      idx = idx.na[i]
      
      centers = centers_ll[idx.not.na,]
      dmat    = rdist(matrix(centers_ll[idx,], nrow=1) , matrix(centers, ncol=2))
      min.val = dmat[1,which.min(dmat[which(dmat>1e-10)])]
      
      idx_close = which(dmat == min.val)
      state     = map.where(database="state", centers[idx_close,1], centers[idx_close,2])
      
      meta$state2[idx] = state
    }
  }
  
  meta$state2[which(meta$state2[idx.mi]=='minnesota')] = 'michigan:north'
  
  idx.bad = which((meta$state2=='michigan:north') & (meta$y<8e5))
  meta$state2[idx.bad] = 'michigan:south'
  
  idx.umw = which(meta$state2 %in% c('michigan:north', 'michigan:south', 'wisconsin', 'minnesota'))
  idx.not.umw = which(!(meta$state2 %in% c('michigan:north', 'michigan:south', 'wisconsin', 'minnesota')))
  if (length(idx.not.umw)>0){
    for (i in 1:length(idx.not.umw)){
      idx = idx.not.umw[i]
      
      centers = centers_ll[idx.umw,]
      dmat    = rdist(matrix(centers_ll[idx,], nrow=1) , matrix(centers, ncol=2))
      min.val = dmat[1,which.min(dmat[which(dmat>1e-10)])]
      
      idx_close = which(dmat == min.val)
      state     = map.where(database="state", centers[idx_close,1], centers[idx_close,2])
      
      meta$state2[idx] = state
    }
  }
  
  
#   meta$state2[which(meta$state2 == 'illinois')] = 'wisconsin'
#   meta$state2[which(meta$state2 == 'north dakota')] = 'minnesota'
#   meta$state2[which(meta$state2 == 'south dakota')] = 'minnesota'
#   meta$state2[which(meta$state2 == 'ohio')] = 'michigan:south'
#   meta$state2[which(meta$state2 == 'indiana')] = 'michigan:south'
  
  return(meta)
  
}


y_build <- function(counts, taxa_sub){ 
  
  three.p = data.frame(read.table(file='data/level3P_v0.2.csv', sep=",", row.names=NULL, header=TRUE, stringsAsFactors=FALSE))
  three.p = rbind(three.p, c('Other conifer', 'Other conifer', 'TRUE'))
  three.p[,1:2] = as.data.frame(apply(three.p[,1:2],2,function(x)gsub('\\s+', '.',x)))
  three.p[,1:2] = as.data.frame(apply(three.p[,1:2],2,function(x)gsub('/', '.',x)))
  
  
  taxa_veg = colnames(counts)
  taxa_use = taxa_sub
  taxa_3p = tolower(three.p[,2])
  
  con = as.logical(three.p[match(taxa_veg, taxa_3p),3])
  
  other.hw.idx = !(taxa_veg %in% taxa_use) & !con
  other.con.idx = !(taxa_veg %in% taxa_use) & con
  
  if (sum(taxa_veg %in% taxa_use) < length(taxa_sub)) print('One or more of the taxa in the provided list is not in 
                                                            the pls data or appears under a different name.')
  
  y       = counts[, taxa_veg %in% taxa_use]
  
  # other is super annoying
  
  y$other.hardwood = other_build(y, counts, other.hw.idx)
  y$other.conifer  = other_build(y, counts, other.con.idx)
  
  y       = unname(round(as.matrix(y)))
  
  return(y)
}
