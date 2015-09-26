library(gridExtra)
library(maptools)

us.shp <- readShapeLines('data/map_data/us_alb.shp',
                         proj4string=CRS('+init=epsg:3175'))
us.shp@data$id <- rownames(us.shp@data)
us.fort <- fortify(us.shp, region='id') 

# # trace plots
# # fit is a stanfit object
# trace_plots <- function(fit, K, suff=suff, save_plots=TRUE){
#   
#   W = K-1 
#   post = extract(fit, permuted=FALSE, inc_warmup=FALSE)
#   
# #   plotname = make_plotname(type_pars, dim_pars, suff)
#   
#   labels = colnames(post[,1,])
#   for (i in 1:length(labels)){
#     labels[i]=sub("[","",labels[i],fixed=TRUE) 
#     labels[i]=sub("]","",labels[i],fixed=TRUE)
#   }
#     
#   par(mfrow=c(4,2))
#     
#   if (save_plots){
#     pdf(paste('veg/figures/trace_pars_', suff, '.pdf', sep=""), width=6, height=4)
#     par(mfrow=c(1,1))
#   }
#     
#   for (i in 1:(4*W)){
#       #set the y limits!!!
#     plot(post[,1,i], type="l", ylab=labels[i], xlab="iter")
#     if (dim(post)[2] >= 2){
#       lines(post[,2,i], col="blue")
#     }
#   }
# #   mtext(suff,outer=TRUE,line=1) 
#     
#   if (save_plots){
#     dev.off()
#   }
# 
# }


# trace plots
# fit is a stanfit object
trace_plots <- function(post_dat, K, suff=suff, save_plots=TRUE, fpath){
  
  one.eta = FALSE
  W = K-1 
  
  par_names = post_dat$par_names
  post = post_dat$samples
  
  if (length(which(par_names == 'eta')) == 1){
    one.eta = TRUE
  } 
  
  if (length(which(par_names == 'rho')) == K){
    W = K
  }
  
  labels = colnames(post)
  for (i in 1:length(labels)){
    labels[i]=sub("[","",labels[i],fixed=TRUE) 
    labels[i]=sub("]","",labels[i],fixed=TRUE)
  }
  
  par(mfrow=c(4,2))
  
  if (save_plots){
    pdf(paste(fpath, '/trace_pars_', suff, '.pdf', sep=""), width=6, height=8)
    par(mfrow=c(3,1))
  }
  
  if (!one.eta){
    for (i in 1:W){
      for (j in 1:3){
        Sys.sleep(0.1)
        plot(post[,(j-1)*W + i], type="l", ylab=labels[(j-1)*W + i], xlab="iter")
#         if (dim(post)[2] >= 2){
#           lines(post[,i], col="blue")
#         }
      }
    }
  } else {
    plot(post[,1,1], type="l", ylab=labels[1], xlab="iter")
    for (i in 1:W){
      for (j in 1:2){
        Sys.sleep(0.1)
        plot(post[,(j-1)*W + i + 1], type="l", ylab=labels[(j-1)*W + i + 1], xlab="iter")
#         if (dim(post)[2] >= 2){
#           lines(post[,2,i], col="blue")
#         }
      }
    }
  }
  #   mtext(suff,outer=TRUE,line=1) 
  
  if (save_plots){
    dev.off()
  }
  
}

# trace plots
# fit is a stanfit object
trace_plots_transform <- function(fit, K, taxa, burnin=25, suff=suff_fit, save_plots=TRUE, fpath){
  
  W = K-1 
  post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  niters = dim(post)[1]
  burnin = burnin
  
  col_substr = substr(colnames(post[,1,]),1,2)
  
  #sqrt_eta = summary(fit)$summary[,'mean'][which(col_substr == 'sq')]
  idx_eta = which(col_substr == 'et')
  idx_rho = which(col_substr == 'rh')
#   sigma = summary(fit)$summary[,'mean'][which(col_substr == 'si')]
#   mu = summary(fit)$summary[,'mean'][which(col_substr == 'mu')]

  #   plotname = make_plotname(type_pars, dim_pars, suff)
  
#   labels = colnames(post[,1,])
#   for (i in 1:length(labels)){
#     labels[i]=sub("[","",labels[i],fixed=TRUE) 
#     labels[i]=sub("]","",labels[i],fixed=TRUE)
#   }
  
  
  cex_lab=1
  if (save_plots){
    pdf(paste(fpath, '/trace_pars_transform_', suff, '.pdf', sep=""), width=12, height=8)
    par(cex.axis=0.8, cex.lab=0.8, cex.main=0.8, cex.sub=1)
  }
  
  for (i in 1:W){
    
    eta_post = post[,1,idx_eta[i]][burnin:niters]
    rho_post = post[,1,idx_rho[i]][burnin:niters]
    
    #set the y limits!!!
    par(mar=c(5,4,4,2)+0.2)
#     par(oma=c(0,0,1,0))
    par(mfcol=c(2,2))
    plot(eta_post, type="l", ylab=paste('eta', i, sep=''), xlab="iter")
    plot(rho_post, type="l", ylab=paste('rho', i, sep=''), xlab="iter")
    plot(log(eta_post/rho_post), type="l", ylab=paste('log(eta', i, '/rho', i, ')', sep=''),, xlab="iter")
    plot(log(eta_post) + log(rho_post), type="l", ylab=paste('log(eta', i, ') + log(rho', i, ')', sep=''), xlab="iter")
#     plot(log(eta_post) - log(rho_post), type="l", ylab=paste('log(eta) - log(rho)', i, sep=' '), xlab="iter")
    title(main=paste('Taxa', i, ';', taxa[i], sep=' '), outer=TRUE, line=-2)
    
  }
  #   mtext(suff,outer=TRUE,line=1) 
  
  if (save_plots){
    dev.off()
  }
  
}


trace_plots_props <- function(r, suff=suff_fit, save_plots, fpath=subDir){
  
  K = dim(r)[1] # n taxa
  
  cols = rep('black', K)#rainbow(K)
  
  par(mfrow=c(4,2))
  if (save_plots){
    pdf(paste(fpath, "/trace_props_", suff, ".pdf", sep=""), width=8, height=4)
    par(mfrow=c(1,1))
  }
  
  for (i in 1:10){
    for (k in 1:K){
      plot(r[k,i,], type="l", col=cols[k], ylab=paste('r',i, 'taxon', k, sep=' '), 
           xlab="iter", ylim=c(min(r[k,i,]),max(r[k,i,])))
    }
  }
  
  
  if (save_plots){
    dev.off()
  }
}

theme_clean <- function(plot_obj){
  plot_obj <- plot_obj + theme(axis.ticks = element_blank(), 
                               axis.text.y = element_blank(), 
                               axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank())
  
  return(plot_obj)
}

# plot_data_maps(y, centers_pls, taxa, K, thresh=1, suff_fit, save_plots)
plot_data_maps <- function(y, centers, taxa, K, thresh, suff, save_plots, fpath){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  prop_dat = t(apply(y, 1, function(x) x/sum(x)))
  colnames(prop_dat) = taxa
  
  prop_dat = data.frame(prop_dat, x = centers[,1], y = centers[,2])
  prop_dat = melt(prop_dat, id.vars = c('x','y'))
  
  if (!is.na(thresh)){
    prop_dat$value[which(prop_dat$value > thresh)] = thresh
  }
  
#   d <- ggplot(data=elm, aes(x=x, y=y)) + geom_tile(aes(x=x, y=y, fill=value)) #+ scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
#   d <- d + facet_wrap(~variable, ncol=5)
# #   d <- theme_clean(d) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5))) +
# #     ggtitle("Observed")
#   
  d <- ggplot() + geom_raster(data=prop_dat, aes(x=x, y=y, fill=value)) + scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
  d <- d + facet_wrap(~variable, ncol=5)
  d <- theme_clean(d) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5))) +
    ggtitle("Observed")
  print(d)
  
  Sys.sleep(1)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_data_maps_', suff, '.pdf', sep=''), scale=1)
  }
  #   print(p)
}

plot_pred_maps <- function(r_mean, centers, taxa, K, thresh, suff, save_plots, fpath){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  prop_dat = t(r_mean)
  colnames(prop_dat) = taxa
  
  prop_dat = data.frame(prop_dat, x = centers[,1], y = centers[,2])
  prop_dat = melt(prop_dat, c('x','y'))
  
  if (!is.na(thresh)){
    prop_dat$value[which(prop_dat$value > thresh)] = thresh
  }
  
  d <- ggplot() + geom_raster(data=prop_dat, aes(x=x, y=y, fill=value)) + scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
  d <- d + facet_wrap(~variable, ncol=5)
  d <- theme_clean(d) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5))) +
    ggtitle("Predicted")
  print(d)
  
  Sys.sleep(1)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_pred_maps_', suff, '.pdf', sep=''), scale=1)
  }
  #   print(p)
}

plot_process_maps <- function(g, centers, taxa, K, thresh, suff, save_plots, fpath){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  W = dim(g)[1]
  
  dat = t(g)
  colnames(dat) = taxa[1:W]
  
  dat = data.frame(dat, x = centers[,1], y = centers[,2])
  dat = melt(dat, c('x','y'))
  
  if (!is.na(thresh)){
    dat$value[which(dat$value > thresh)] = thresh
  }
  
  d <- ggplot() + geom_raster(data=dat, aes(x=x, y=y, fill=value)) + scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
  d <- d + facet_wrap(~variable, ncol=5)
  d <- theme_clean(d) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5))) +
    ggtitle("Process")
  print(d)
  
  Sys.sleep(1)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_process_maps_', suff, '.pdf', sep=''), scale=1)
  }
  #   print(p)
}

plot_data_maps_binned <- function(y, centers, taxa, K, breaks, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  y = t(y)
  
  props_data = t(apply(y, 1, function(x) if (sum(x) > 0) {x/sum(x)} else {x}))
  #colnames(props_data) = taxa
  
  props_data_binned = matrix(0, nrow=nrow(props_data), ncol=ncol(props_data))
  colnames(props_data_binned) <- colnames(props_data)
  
  for (i in 1:ncol(props_data)){
    props_data_binned[,i] = cut(props_data[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data_binned[,k], 
                                          x     = centers[,1]*rescale, 
                                          y     = centers[,2]*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
  prop_dat$type = rep('PLS', nrow(prop_dat))
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props))) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits)
  p <- p + facet_wrap(~taxon, ncol=6)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
#   p <- p + theme(strip.text.x = element_blank(),
#                  strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  #print(p)
  #Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.eps', sep=''), scale=1)
    #     dev.off()
  }
  return(p)
}

add_map_albers <- function(plot_obj, map_data=us.fort, limits){
  p <- plot_obj + geom_path(data=map_data, aes(x=long, y=lat, group=group),  colour='grey55') + 
    #     scale_x_continuous(limits = c(min(umw.coord$x, na.rm=TRUE), max(umw.coord$x, na.rm=TRUE))) +
    #     scale_y_continuous(limits = c(min(umw.coord$y, na.rm=TRUE), max(umw.coord$y, na.rm=TRUE)))#, colour = "black", size = 1, fill = "white", aes(x=long, y=lat, group = group))
    # #   
    #     scale_x_continuous(limits = c(min(dat[,1], na.rm=TRUE), max(dat[,1], na.rm=TRUE))) +
    #     scale_y_continuous(limits = c(min(dat[,2], na.rm=TRUE), max(dat[,2], na.rm=TRUE)))
    scale_x_continuous(limits = limits$xlims*1000000) +
    scale_y_continuous(limits = limits$ylims*1000000) #+ coord_map("albers")
  return(p)
  
}

get_limits <- function(centers){
  xlo = min(centers[,1])
  xhi = max(centers[,1])
  
  ylo = min(centers[,2])
  yhi = max(centers[,2])
  
  return(list(xlims=c(xlo,xhi),ylims=c(ylo, yhi)))
}  

# 
# plot_Halpha_maps <- function(Halpha, centers, taxa, K, thresh, suff, save_plots, fpath){
#   
#   if (is.null(taxa)){taxa=seq(1,K)}
#   
#   W = dim(Halpha)[1]
#   
#   dat = t(g)
#   colnames(dat) = taxa[1:W]
#   
#   dat = data.frame(dat, x = centers[,1], y = centers[,2])
#   dat = melt(dat, c('x','y'))
#   
#   if (!is.na(thresh)){
#     dat$value[which(dat$value > thresh)] = thresh
#   }
#   
#   d <- ggplot() + geom_raster(data=dat, aes(x=x, y=y, fill=value)) + scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
#   d <- d + facet_wrap(~variable, ncol=5)
#   d <- theme_clean(d) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5))) +
#     ggtitle("Process")
#   print(d)
#   
#   Sys.sleep(1)
#   if (save_plots){
#     ggsave(file=paste(fpath, '/veg_process_maps_', suff, '.pdf', sep=''), scale=1)
#   }
#   #   print(p)
# }