library(ggplot2)
library(rstan)
library(reshape2) # do i need reshape?
library(fields)

subDir = 'figures'

source('r/utils/veg_process_funs.r')
source('r/utils/veg_plot_funs.r')
source('r/read_stanbin.r')

#####################################################################################################################

suff_data = run$suff_data
suff_fit  = run$suff_fit

# suff_data = '12taxa_6341cells_120knots_v0.4'
# suff_fit  = '12taxa_6341cells_120knots_v0.4_c1'

decomp     = TRUE
bt         = TRUE
mpp        = TRUE
save_plots = TRUE


fname = sprintf('output/%s.bin', suff_fit)


subDir1 = paste(subDir, '/', suff_fit, sep='')
if (!file.exists(subDir1)){
  dir.create(file.path(getwd(), subDir1))
}

# load the data and posterior draws
load(paste0('r/dump/veg_data_', suff_data, '.rdata'))
post_dat = read_stanbin(fname)
post_dat$par_names = post_dat$par_names[7:length(post_dat$par_names)]
post_dat$samples = post_dat$samples[,7:ncol(post_dat$samples)]
#post_dat = load_stan_output(suff_fit)

post = post_dat$samples
 
W = K-1

# compute ess and write to file
sink(sprintf('%s/ess.txt', subDir1), type='output')
print('Effective samples sizes : ')
# ess_all = rowSums(sapply(post_dat$samples[,1:(W*3)], function(y) apply(y, 2, function(x) ess_rfun(x))))
ess_all = apply(post_dat$samples[,1:(W*3)], 2, function(x) ess_rfun(x))
print(as.matrix(ess_all[!(sapply(strsplit(names(ess_all),'\\['), function(x) x[[1]]) == 'log_lik')]))
sink()

# compute par summary and write to file
sink(sprintf('%s/summary.txt', subDir1), type='output')
print('Effective samples sizes : ')
# ess_all = rowSums(sapply(post_dat$samples[,1:(W*3)], function(y) apply(y, 2, function(x) ess_rfun(x))))
mean_pars = t(t(colMeans(post)))
colnames(mean_pars) = "mean"
quants    = t(apply(post, 2, function(x) quantile(x, probs=c(0.5, 0.025, 0.975))))
print(cbind(mean_pars, quants))
sink()


post=post[,1:(3*W)]
save(post, file=paste0(subDir1, '/veg_pars_', N_knots, 'knots.rdata'))

# ess(fit, W)

# Rprof('compute_prop_chains_Halpha.out')
#out = compute_prop_chains(fit, d_knots, d_inter, K, P, decomp=decomp, bt=bt, mpp=mpp) 
out = compute_prop_chains_Halpha(post_dat, d_knots, d_inter, K, P, decomp=decomp, bt=bt, mpp=mpp) 
# Rprof(NULL)  
# summaryRprof("compute_prop_chains_Halpha.out")

# benchmark(compute_prop_chains_Halpha(fit, d_knots, d_inter, K, P, decomp=decomp, bt=bt, mpp=mpp) , compute_prop_chains_HalphaR(fit, d_knots, d_inter, K, P, decomp=decomp, bt=bt, mpp=mpp) )

r = out$r
g = out$g
Halpha= out$Halpha

sumHalpha = array(NA, c(W, dim(post_dat$samples)[1]))

for (k in 1:W){
  sumHalpha[k,] = colSums(Halpha[k,,])
}

mu = get_mu(post_dat, W)

pdf(file=paste0(subDir1, '/compare_mu_Halpha.pdf'), width=8, height=8)
for (k in 1:W){
  par(mfrow=c(2,1))
  par(oma=c(0,0,2,0))
  plot(mu[,k], type="l", ylab=paste0('mu[', k, ']'))
  #   lines(mu_t[t+1,1,], col="blue")
  plot(sumHalpha[k,], type="l", ylab=paste0('sum_Halpha[', k, ']'))
  title(main=taxa[k], outer=TRUE)
}
dev.off()

# save(r, file=paste('veg/output/', suff_fit, '_r.rdata', sep=''))
# 
# load(file=paste('veg/output/', suff_fit, '_r.rdata', sep=''))

r_mean = compute_prop_means(r)
g_mean = compute_prop_means(g)
Halpha_mean = compute_prop_means(Halpha)

rm(Halpha)
rm(g)
gc()

####################################################################################################
# chunk: trace plots
####################################################################################################

suff_fig = paste0(N_knots, 'knots')
trace_plots(post_dat, K, suff=suff_fig, save_plots=save_plots, fpath=subDir1)

suff_fig = paste0(N_knots, 'knots')
trace_plots_props(r, suff=suff_fig, save_plots, fpath=subDir1)
rm(r)

# trace_plots_transform(fit, K, taxa, burnin=25, suff=suff_fit, save_plots=TRUE, fpath=subDir)

####################################################################################################
# chunk: plot observed proportions
####################################################################################################
suff_fig = paste0(N_knots, 'knots')
plot_process_maps(g_mean, centers_pls, taxa, K, thresh=1, suff_fig, save_plots, fpath=subDir1)

suff_fig = paste0(N_knots, 'knots_Halpha')
plot_process_maps(Halpha_mean, centers_pls, taxa, K, thresh=1, suff_fig, save_plots, fpath=subDir1)

####################################################################################################
# chunk: plot observed proportions
####################################################################################################

suff_fig = paste0(N_knots, 'knots')
plot_data_maps(y, centers_pls, taxa, K, thresh=1, suff_fig, save_plots, fpath=subDir1)

####################################################################################################
# chunk: plot predicted proportions
####################################################################################################

suff_fig = paste0(N_knots, 'knots')
plot_pred_maps(r_mean, centers_pls, taxa, K, thresh=1, suff_fit, save_plots, fpath=subDir1)

limits <- get_limits(centers_pls)

suff_fig = paste0('props_binned_', N_knots, 'knots')
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# plot_data_maps_binned(r_mean, centers=centers_pls, taxa=taxa, K, T, breaks, suff=suff4, save_plots=save_plots)
plot_data_maps_binned(r_mean, centers_pls, taxa, K, breaks, limits, suff=suff_fig, save_plots, fpath=subDir1)

# }
