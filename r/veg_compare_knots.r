# compare veg preds for different knot numbers

library(ggplot2)
library(reshape2)
library(fields)
source('r/read_stanbin.r')

save_plots = TRUE

run_40knots_c1 = list(suff_data = '12taxa_6341cells_40knots_v0.4',
                      suff_fit  = '12taxa_6341cells_40knots_v0.4_c1')
run_80knots_c1 = list(suff_data = '12taxa_6341cells_80knots_v0.4',
                      suff_fit  = '12taxa_6341cells_80knots_v0.4_c1')
run_120knots_c1 = list(suff_data = '12taxa_6341cells_120knots_v0.4',
                       suff_fit  = '12taxa_6341cells_120knots_v0.4_c1')
run_160knots_c1 = list(suff_data = '12taxa_6341cells_160knots_v0.4',
                       suff_fit  = '12taxa_6341cells_160knots_v0.4_c1')
run_200knots_c1 = list(suff_data = '12taxa_6341cells_200knots_v0.4',
                       suff_fit  = '12taxa_6341cells_200knots_v0.4_c1')
run_220knots_c1 = list(suff_data = '12taxa_6341cells_220knots_v0.4',
                       suff_fit  = '12taxa_6341cells_220knots_v0.4_c1')
run_240knots_c1 = list(suff_data = '12taxa_6341cells_240knots_v0.4',
                       suff_fit  = '12taxa_6341cells_240knots_v0.4_c1')
run_260knots_c1 = list(suff_data = '12taxa_6341cells_260knots_v0.4',
                       suff_fit  = '12taxa_6341cells_260knots_v0.4_c1')


runs = list(run_40knots_c1,
            run_80knots_c1,
            run_120knots_c1,
            run_160knots_c1,
            run_200knots_c1,
            run_220knots_c1,
            run_240knots_c1,
            run_260knots_c1)


eta = data.frame(eta=numeric(0), taxon=character(0), knots=numeric(0))
rho = data.frame(rho=numeric(0), taxon=character(0), knots=numeric(0))
mu  = data.frame(mu=numeric(0), taxon=character(0), knots=numeric(0))
g   = data.frame(g=numeric(0), cell=numeric(0), taxon=character(0), knots=numeric(0))

for (run in runs) {

  suff_data = run$suff_data
  suff_fit  = run$suff_fit
  
  # load the data
  load(paste0('r/dump/veg_data_', suff_data, '.rdata'))
  W = K-1
  
#   # load the data
#   load(paste0('figures/', suff_fit, '/veg_pars_', N_knots, 'knots.rdata'))
  if (!file.exists(sprintf('output/%s.bin', suff_fit))) {
    next
  }
  
  post_dat  = read_stanbin(sprintf('output/%s.bin', suff_fit))
  par_names = post_dat$par_names[7:length(post_dat$par_names)]
  post      = post_dat$samples[,7:ncol(post_dat$samples)]
  
  eta = rbind(eta, data.frame(eta   = matrix(colMeans(post[,which(par_names == 'eta')])), 
                              taxon = matrix(taxa[1:W]), 
                              knots = rep(N_knots, W)))
  
  rho = rbind(rho, data.frame(rho   = matrix(colMeans(post[,which(par_names == 'rho')])), 
                              taxon = matrix(taxa[1:W]), 
                              knots = rep(N_knots, W)))
  
  mu = rbind(mu, data.frame(mu    = matrix(colMeans(post[,which(par_names == 'mu')])), 
                              taxon = matrix(taxa[1:W]), 
                              knots = rep(N_knots, W)))
  
  for (k in 1:K){
    g_cols = seq((3*W + k + W*N_knots), (3*W + W*N_knots + W*N), by=W) 
    g_mean = rowMeans(t(post[,g_cols]))
    cell   = sapply(strsplit(names(g_mean), '\\.'), function(x) x[[3]])
    g = rbind(g, data.frame(g    = g_mean,
                            cell = cell,
                            taxon = rep(taxa[k], length(cell)),
                            knots = rep(N_knots, length(cell))))
  }
  
}

p <- ggplot(data=eta) + geom_point(data=eta, size=4, aes(x=taxon, y=eta, colour=factor(knots))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p, file='figures/eta.pdf')

p <- ggplot(data=rho) + geom_point(data=rho, size=4, aes(x=taxon, y=rho, colour=factor(knots))) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p, file='figures/rho.pdf')

p <- ggplot(data=mu) + geom_point(data=mu, size=4, aes(x=taxon, y=mu, colour=factor(knots))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p, file='figures/mu.pdf')

knot_vals = sort(unique(g$knots))
n_vals = length(knot_vals)

# for (i in 1:(length(knot_vals)-1) {
#   for (k in 1:K) {
#     taxon = taxa[k]
#     g[which(g$knots == knot_vals[i]) ,]
#     foo=g[which(g$knots == knot_vals[i]) ,]
#   }
# }

g_wide = dcast(g, cell + taxon ~ knots, mean, value.var='g')

# # dis = data.frame(cell=numeric(0), taxon=character(0), dis=numeric(0))
# dis = g_wide[,1:2]
# 
# for (i in 1:(n_vals-1)) {
#   knots = knot_vals[i]
#   
#   diff = data.frame((g_wide[,as.character(knots)] - g_wide[,as.character(knot_vals[n_vals])])^2)
#   colnames(diff) = paste0('diff_', knots)
#   
#   dis = cbind(dis, diff)
# }

dis = data.frame(matrix(NA, nrow=N, ncol=n_vals))
colnames(dis) = c('cell', knot_vals[1:(n_vals-1)])

for (j in 1:N) {
  dat = g_wide[g_wide$cell == j,]
  
  dis[j,1] = j
    
  for (i in 1:(n_vals-1)) {
    knots = knot_vals[i]
    
    d = sqrt(sum((dat[, as.character(knots)] - dat[,as.character(knot_vals[n_vals])])^2))
    
    dis[j,i+1] = d
    
  }
}

colnames(centers_pls) = c('x', 'y')
dis = cbind(dis, centers_pls)
dis=melt(dis, id.vars=c('cell', 'x', 'y'))
dis=data.frame(dis)

p <- ggplot(data=dis) + geom_raster(data=dis, aes(x=x, y=y, fill=value)) + scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
p <- p + facet_grid(variable~.)
p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5))) 
# print(d)
ggsave(p, file='figures/knot_diss.pdf')

dis_tot = aggregate(value ~ variable, dis, sum)
p<- ggplot() + geom_point(data=dis_tot, aes(x=variable, y=value))
ggsave(p, file='figures/dis_total.pdf')
