library(fields)
library(mvtnorm)
library(sp)
library(rgdal)

source('r/utils/veg_pp_helpers.r')
source('r/utils/veg_build_data_funs.r')

# suff = 'wi'
# suff='wi'
# suff='rescaled_dist'
#suff="v0.3"
suff='v0.4'

one.eta = FALSE
bt      = TRUE
  
# state can be: mi, wi, mn
# states_pls = c('minnesota', 'wisconsin','michigan_north', 'michigan_south')
states_pls = c('minnesota', 'wisconsin','michigan:north')
# states_pls = 'michigan:south'
#states_pls = 'minnesota'
cells = -1
# cells = seq(1001, 2500)
#cells = seq(1001, 1050)

clust.rat = 120
#nclust    = 120
#nclust    = 80
nclust = nclust
# clust.rat = 10

# specify the taxa to use
# must be from the list: taxa_sub = c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech', 'elm', 'spruce', 'ash', 'hemlock')
# always have: 'other.hardwood' and 'other.conifer'

taxa_sub = c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech', 'elm', 'spruce', 'ash', 'hemlock')
#taxa_sub = c('oak', 'pine') #, 'maple', 'birch')
#taxa = c(taxa_sub, 'other')

# conversion table
tree_type = read.table('data/assign_HW_CON.csv', sep=',', row.names=1, header=TRUE)

##########################################################################################################################
## chunk 1: read in and subset pls and pollen data
##########################################################################################################################
# pls.raw     = data.frame(read.table(file='data/western_comp_stepps_v0.3-1.csv', sep=",", row.names=NULL, header=TRUE))
# pls.raw = data.frame(read.table(file='data/composition_v0.3.csv' , sep=",", row.names=NULL, header=TRUE))
# pls.raw = data.frame(read.table(file='../stepps-data/data/composition/composition_v0.4.csv' , sep=",", row.names=NULL, header=TRUE))
# colnames(pls.raw) = tolower(colnames(pls.raw))

pls.raw = data.frame(read.table(file='../stepps-data/data/pls_umw_v0.5.csv' , sep=",", row.names=NULL, header=TRUE))
colnames(pls.raw) = tolower(colnames(pls.raw))

# figure out why chestnut and atlantic.white.cedar are all NA!
pls.raw = pls.raw[,!(colnames(pls.raw) %in% c('chestnut', 'atlantic.white.cedar'))]

# conversion table
convert = read.table('data/dict-comp2stepps.csv', sep=',', row.names=1, header=TRUE)

# pull the subset of proportions
taxa.start.col = min(match(tolower(rownames(convert)), colnames(pls.raw)), na.rm=TRUE)

pls_dat  = pls.raw[,taxa.start.col:ncol(pls.raw)]
colnames(pls_dat) = as.vector(convert[match(colnames(pls_dat), tolower(rownames(convert))),1])
pls_dat_collapse  = sapply(unique(colnames(pls_dat)), 
                            function(x) rowSums( pls_dat[ , grep(x, names(pls_dat)), drop=FALSE]) )
counts = data.frame(pls_dat_collapse[,sort(colnames(pls_dat_collapse))])
meta   = pls.raw[,1:(taxa.start.col-1)]
# kilometers
# pls$X = pls$X/1000
# pls$Y = pls$Y/1000

meta_tmp = split_mi(meta, longlat=FALSE)

plot(meta_tmp$x, meta_tmp$y)

# filter
counts = counts[which(meta_tmp$state2 %in% states_pls),]
meta   = meta_tmp[which(meta_tmp$state2 %in% states_pls),]
if (length(cells) > 1){
  counts = counts[cells,]
  meta   = meta[cells,]
}

points(meta$x, meta$y, col='blue')

# megameters!
meta$x = meta$x/1000000
meta$y = meta$y/1000000

# tree_type_names = tolower(rownames(tree_type))

# left = rownames(tree_type)[!(rownames(tree_type) %in% taxa_sub)]
left = rownames(tree_type)[!(rownames(tree_type) %in% toupper(taxa_sub))]
y = data.frame(counts[, colnames(counts) %in% toupper(taxa_sub)])

taxa_other_hw = rownames(tree_type)[which(tree_type$type == 'HW')]
taxa_other_con = rownames(tree_type)[which(tree_type$type == 'CON')]

if (sum(left %in% taxa_other_hw)>1){
  y$OTHER.HARDWOOD = rowSums(counts[,left[left %in% taxa_other_hw]])
} else {
  y$OTHER.HARDWOOD = counts[,left[left %in% taxa_other_hw]]
}

if (sum(left %in% taxa_other_con)>1){
  y$OTHER.CONIFER = rowSums(counts[,left[left %in% taxa_other_con]])
} else {
  y$OTHER.CONIFER = counts[,left[left %in% taxa_other_con]]
}

y = y[,sort(colnames(y))]

taxa = colnames(y)
y    = unname(round(as.matrix(y)))
#rownames(y) = NULL
# y = y_build(counts, taxa_sub) # fix this if we want to use a subset of taxa

K = as.integer(ncol(y))
W = K-1
N = nrow(y)

centers_pls = cbind(meta$x, meta$y)

# distance matrix
d = rdist(as.matrix(centers_pls), as.matrix(centers_pls))
diag(d) <- 0

##########################################################################################################################
## determine knot locations
##########################################################################################################################
# source('r/utils/simDataFuns.r')
# source('r/data/utils/dataPlotFuns.r')

if (is.na(nclust)){
  nclust = ceiling(N/clust.rat)
}
  
knot_coords = kmeans(centers_pls, nclust)$centers
knot_coords = unname(knot_coords)
N_knots     = dim(knot_coords)[1]

# subgrid = regular_subgrid(centers_pls, dx=60,dy=60)
# knot_coords = knots_in_domain4(subgrid, centers_pls, cell_width = 8)
# N_knots = dim(knot_coords)[1]
# par(mfrow=c(1,1))
# # png(filename='domain_cells.png')
# pdf(file='domain_cells.pdf', width=8, height=8)
# plot(centers_pls[,1], centers_pls[,2], asp=1, axes=F, col='antiquewhite4',xlab='',ylab='', pch=19, cex=0.2)
# dev.off()
# 
# par(mfrow=c(1,1))
# # png(filename='domain_knots.png')
# pdf(file='domain_knots.pdf', width=8, height=8)
plot(centers_pls[,1], centers_pls[,2], asp=1, axes=F,  col='antiquewhite4', xlab='',ylab='', pch=19, cex=0.2)
points(knot_coords[,1], knot_coords[,2], col='red', pch=19, cex=0.5)
# dev.off()

# distance matrix
d_knots = rdist(knot_coords, knot_coords)
diag(d_knots) <- 0

# distance matrix
d_inter = rdist(centers_pls, knot_coords)
d_inter[which(d_inter<1e-8)]=0

##########################################################################################################################
## chunk: qr decompose X
##########################################################################################################################
x = matrix(1, nrow=N, ncol=1)
N_p = N

temp = qr(x)
Q = qr.Q(temp)
R = qr.R(temp)

P = Q %*% t(Q)
# M = diag(N) - P

if (all(P-P[1,1]<1.0e-12)){
  P = P[1,1]
  N_p = 1
}
##########################################################################################################################
## chunk: save data to file
##########################################################################################################################
if (suff != "") {suff = paste('_', suff, sep='')}

dump(c('K', 'N', 'N_knots', 
       'y', 
       'd_knots', 'd_inter',
       'P', 'N_p'), 
     file=paste('r/dump/veg_data_', K, 'taxa_', N, 'cells_', N_knots, 'knots', suff, '.dump',sep=""))

save(K, N, N_knots, 
     y, 
     d, d_knots, d_inter, 
     P, N_p,
     centers_pls, knot_coords, taxa,
     file=paste('r/dump/veg_data_', K, 'taxa_', N, 'cells_', N_knots, 'knots', suff, '.rdata',sep=""))

##########################################################################################################################
## chunk: generate and save initial conditions
##########################################################################################################################

if (bt){
  W = K-1
} else {
  W = K
  suff=paste(suff,'_nbt', sep='')
}

if (one.eta){
  eta   = c(1) 
  sqrt_eta   = c(1)
  
  suff=paste(suff,'_1eta', sep='')
  
} else {
  eta   = rep(1, W) 
  sqrt_eta   = rep(1, W)
}

rho   = rep(0.2, W)
mu    = rep(0,W) 

sigma = rep(0.5, W)

init_out = veg_build_inits(K, N, N_knots, eta, rho, mu, d_knots, d_inter)
alpha = init_out$alpha_init
alpha_raw = alpha
g = init_out$g_init

dump(c('eta', 'sqrt_eta', 'rho', 'mu', 'sigma', 'alpha', 'alpha_raw', 'g'), 
     file=paste('r/dump/veg_data_', K, 'taxa_', N, 'cells_', N_knots, 'knots', suff, '_inits.dump',sep=""))