### SAR model with independent normal-gamma prior and beta prior for rho ###
###  implements the sigma prior of Simpson, Rue et al.
rm(list=ls())

require(Matrix)
require(wbstats)
require(truncnorm)
require(ggplot2)
require(MASS)
require(dplyr)
require(sf)

# Fast logdet approximation function
source("./Estimation/lndetPaceBarry.R")

load("./Estimation/xy.RData")

#########--------------------------------------------------------------
W.list = W.list[c("trade","Queen","flight","FTA")]
#########--------------------------------------------------------------
#correct data in shpfile
y <- shp %>% 
  dplyr::select(starts_with("2020-")) %>% 
  st_drop_geometry %>%
  as.matrix()
#per capita
wbpop = wb(indicator = "SP.POP.TOTL",startdate = "2018",enddate = "2018")
wbpop$value = wbpop$value / 100000
wbpop = wbpop[match(shp$iso3c,wbpop$iso3c),]
for (j in 1:ncol(y)) {
  y[,j] = y[,j] / wbpop$value
}
N <- nrow(y)
y[is.na(y)] <- 0
y <- t(y)
#y <- log(y + .01)
y <- log(y + min(y[y>0])/100)
#y <- diff(y)
X <- apply(y, 2, lag)
X <- na.omit(X)
y <- y[2:(nrow(y)),]
T <- nrow(y)

X <- t(X)
y2 = t(y)
y <- as.vector(y2)
Xraw <- X

X1 <- as.vector(t(Xraw))

countries = rep(as.factor(shp$iso3c),T)
country_dummies = model.matrix(~. + 0,data = data.frame(countries))

#summary(lm(y~X1 + country_dummies + 0))
WY <- y2
mix <- rep(1,T)
for(jj in 1:ncol(y2)){
  WY[,jj] <- W.list[[mix[jj]]] %*% WY[,jj]
}

#summary(lm(y~X1 + country_dummies + as.vector(WY) + 0))

as.bdiag = function (x) {bdiag(split(x, rep(1:ncol(x), each = nrow(x))))}

#summary(lm(y~X1 + country_dummies + as.matrix(as.bdiag(WY)) + 0))




X <- cbind(as.vector(Xraw), country_dummies)

Y = as.matrix(y)

ddd = length(W.list)

smalln = nrow(y2)
k = ncol(X)
n = nrow(X)

### let us assign prior values
# beta mean and variance are proper, but with high variance (= low prior information)
beta_prior_mean = matrix(0,k,1)
beta_prior_var = diag(k) * 10^8
# rho prior is beta, with 
beta_prob = function(rho,a) 1/beta(a,a) * ((1+rho)^(a-1) *(1-rho)^(a-1))/(2^(2*a - 1))
rho_a =1.01
# sigma prior is accrding to Simpson, 
sigma_prob = function(sig,lambda) (lambda/2*sig^(-3/2))*exp(-lambda*sig^(-1/2))
lambda_sig = 1
sigma_prop = "truncnorm"


### calibration parameters for rho sampling
ccs = rep(1,T) #scaling of rho proposals
c_adjust = 1.1 #proposal distribution adjustment
rhos_accept = seq(0,T) #counter for rho acceptance rates

### set-up the gibbs sampler
# total number of draws
niter = 2000
# retain only S2 and discard S1, with S = S1 + S2
nretain = 1000
ndiscard = niter - nretain
# save the posterior draws here
postb = matrix(0,k,nretain)
posts = matrix(0,1,nretain)
postr = matrix(0,T,nretain)
rtmp = matrix(0,T,niter) # rho acceptance rates
# post.direct = matrix(0,k,nretain)
# post.indirect = matrix(0,k,nretain)
# post.total = matrix(0,k,nretain)
sig_accept = 0
acc.rates <- matrix(0,ncol=(T), nrow=nretain) # W mix accept rate

if (sigma_prop == "truncnorm") {
  c_sig = 10
  c_sig_adjust = 1.1
  sigma_accept_rate = matrix(0,niter, 1)
}

# set-up for griddy gibbs
griddy_n = 500
all_logdets = matrix(0,griddy_n,ddd)
rrhos = seq(-1,1,length.out = griddy_n + 2)[-c(1,griddy_n+2)]
cat("Pre-calculate griddy GIBBS...")
for (d in 1:ddd) {
  for (ii in 1:griddy_n) {
    all_logdets[ii,d] = log(det(diag(smalln) - rrhos[ii] * W.list[[d]]))
  }
}
cat("Done!\n")
# AYs = array(0,c(n,griddy_n))
# ## storage for efficient partial derivatives
# ai_diags = rep(0,griddy_n)
# ai_tots = rep(0,griddy_n)
# cat("Pre-calculate griddy GIBBS...")
# for (ii in 1:griddy_n) {
#   A = (.sparseDiagonal(n) - rrhos[ii] * W)
#   AYs[,ii] = as.matrix(A %*% Y)
#   AI = solve(A)
#   ai_diags[ii] = sum(diag(AI))
#   ai_tots[ii] = sum(AI)
# }
# cat("Done!\n")

# starting values (won't matter after sufficient draws)
#curr.beta = mvrnorm(1,beta_prior_mean,beta_prior_var)
curr.beta = as.matrix(solve(t(X) %*% X) %*% t(X) %*% Y)
curr.sigma = 1
curr.inds = rep(round(griddy_n/2),T)
curr.rhos = rrhos[curr.inds]
curr.mix = sample(1:ddd,T,replace = T)
curr.logdet = sum(all_logdets[cbind(curr.inds, curr.mix)])

# pre-calculate some terms for faster draws
beta_prior_var_inv = solve(beta_prior_var)
XpX = t(X) %*% X
curr.WY = curr.Ay = Y
for (tt in 1:T) {
  ind1 = ((tt-1) *smalln + 1) : (tt*smalln)
  curr.WY[ind1] = W.list[[curr.mix[tt]]] %*% Y[ind1] 
  curr.Ay[ind1] = Y[ind1] - curr.rhos[tt]*curr.WY[ind1]
}


### Gibbs sampling
for (iter in 1:niter) {
  cat("iter:",iter,"rho: min",min(curr.rhos),"mean",mean(curr.rhos),"max",max(curr.rhos),"sigma",curr.sigma,"\n")
  
  # draw beta
  V = solve(beta_prior_var_inv + 1/curr.sigma * XpX )
  b = V %*% (beta_prior_var_inv%*%beta_prior_mean + 1/curr.sigma * t(X) %*% curr.Ay )
  curr.beta = mvrnorm(1,b,V)
  
  # draw sigma using M-H step
  curr.xb = X %*% curr.beta
  curr.ESS = as.double(crossprod(curr.Ay - curr.xb))
  curr.llh = log(1/(sqrt(curr.sigma)*sqrt(2*pi)))*n - curr.ESS / (2*curr.sigma) + log(sigma_prob(curr.sigma,lambda_sig))
  # candidate draw
  if (sigma_prop == "uniform") {
    # uniform proposal would be an alternative
    Um = 1.1 * (2 * runif(1) -1)
    prop.sigma = curr.sigma * exp(Um)
  } else if (sigma_prop == "truncnorm") {
    prop.sigma = rtruncnorm(1,curr.sigma,c_sig,a = 0,b = Inf)
  }
  prop.llh = log(1/(sqrt(prop.sigma)*sqrt(2*pi)))*n - curr.ESS / (2*prop.sigma) + log(sigma_prob(prop.sigma,lambda_sig))
  acc = prop.llh - curr.llh
  if (log(runif(1)) < acc) {
    curr.sigma = prop.sigma
    sig_accept = sig_accept + 1
  }
  if (sigma_prop == "truncnorm") {
    sigma_accept_rate[iter] = sig_accept / iter
    if (iter < ndiscard / 2) {
      if (sig_accept / iter < .3) {
        c_sig = c_sig / c_sig_adjust
      } else if (sig_accept / iter > .6) {
        c_sig = c_sig * c_sig_adjust
      }
    }
  }
  
  
  # ## Griddy-Gibbs step for rho
  # V = solve(beta_prior_var_inv + 1/curr.sigma * XpX )
  # b0 = V %*% (beta_prior_var_inv%*%beta_prior_mean + 1/curr.sigma * t(X) %*% Y )
  # bd = V %*% (beta_prior_var_inv%*%beta_prior_mean + 1/curr.sigma * t(X) %*% curr.WY)
  # e0 = Y - X %*% b0
  # ed = curr.WY - X %*% bd
  # epe0 = as.double(t(e0) %*% e0)
  # eped = as.double(t(ed) %*% ed)
  # epe0d = as.double(t(ed) %*% e0)
  # z = epe0  - 2 * rrhos * epe0d + rrhos^2 * eped
  # z = -(n-k)/2 * log(z)
  # den = rowSums(all_logdets[,curr.mix]) + z + log(beta_prob(rrhos,rho_a))
  # y = rrhos
  # adj = max(den)
  # den = den - adj
  # x = exp(den)
  # isum = sum((y[-1] + y[-length(y)])*(x[-1]  - x[-length(x)])/2)
  # z = abs(x/isum)
  # den = cumsum(z)
  # rnd = runif(1) * sum(z)
  # ind = max(which(den <= rnd))
  # if (is.integer(ind) && ind <= length(rrhos)) {
  #   curr.rho = rrhos[ind]
  #   curr.Ay = Y - curr.rho * curr.WY
  #   #curr.ai_diag = ai_diags[ind]
  #   #curr.ai_tot = ai_tots[ind]
  # }
  
  ## Metropolis-Hastings step for rhos
  curr.llh = sum(all_logdets[cbind(curr.inds, curr.mix)]) - as.double(crossprod(curr.Ay - curr.xb))/ (2*curr.sigma) 
  curr.llh = curr.llh + sum(log(beta_prob(curr.rhos,rho_a)))
  for (tt in 1:T) {
    from_to = ((tt-1)*smalln+1):(tt*smalln)
    accept = 0;
    while (accept!=1) {
      prop.rho = curr.rhos[tt] + ccs[tt]*rnorm(1,0,1)
      if (prop.rho<1 && prop.rho>-1) {
        accept = 1
      }
    }
    prop.rhos = curr.rhos; prop.rhos[tt] = prop.rho
    prop.ind = which(rrhos>=prop.rho)[1]; if (is.na(prop.ind)) {prop.ind = griddy_n}
    prop.inds = curr.inds; prop.inds[tt] = prop.ind
    prop.Ay = curr.Ay
    prop.Ay[from_to] = Y[from_to] - prop.rho*curr.WY[from_to]
    prop.ESS = crossprod(prop.Ay - curr.xb)
    prop.logdet = sum(all_logdets[cbind(prop.inds, curr.mix)])
    prop.llh = prop.logdet - as.double(prop.ESS)/ (2*curr.sigma)
    prop.llh = prop.llh + sum(log(beta_prob(prop.rhos,rho_a)))
    acc_prob = min(1,exp(prop.llh - curr.llh))
    if (rbinom(1,1,acc_prob) == 1) {
      curr.rhos = prop.rhos
      curr.inds = prop.inds
      rhos_accept[tt] = rhos_accept[tt] + 1
      curr.logdet = prop.logdet
      curr.Ay = prop.Ay
    }
    # # adjust candidate distribution based on acceptance probability
    # rtmp[tt,iter] = rhos_accept[tt]/iter
    # if (iter < ndiscard/2) {
    #   if (rtmp[tt,iter]<.3) {
    #     ccs[tt] <- ccs[tt]/c_adjust
    #   } else if (rtmp[tt,iter]>.7) {
    #     ccs[tt] <- ccs[tt]*c_adjust
    #   }
    # }
  }
  
  # sample curr.mix
  curr.xb = X %*% curr.beta
  curr.llh = sum(all_logdets[cbind(curr.inds,curr.mix)]) - as.double(crossprod(curr.Ay - curr.xb))/ (2*curr.sigma)
  for(tt in 1:(T)) {
    sss = (1:ddd)[-curr.mix[tt]]
    prop.mix = curr.mix
    prop.mix[tt] = sample(sss,1)

    prop.WY = curr.WY
    from_to = ((tt-1)*smalln+1):(tt*smalln)
    prop.WY[from_to] = W.list[[prop.mix[tt]]] %*% Y[from_to]
    prop.Ay = curr.Ay
    prop.Ay[from_to] = Y[from_to] - curr.rhos[tt]*prop.WY[from_to]

    prop.llh = sum(all_logdets[cbind(curr.inds,prop.mix)]) - as.double(crossprod(prop.Ay - curr.xb))/ (2*curr.sigma)

    acc_prob = min(1,exp(prop.llh - curr.llh))

    if (rbinom(1,1,acc_prob) == 1) {
      curr.mix <- prop.mix
      curr.WY = prop.WY
      curr.llh = prop.llh
    }

  }
  
  # we are past the burn-in, save the draws
  if (iter > ndiscard) {
    s = iter - ndiscard
    postb[,s] = as.matrix(curr.beta)
    posts[s] = curr.sigma
    postr[,s] = curr.rhos
    acc.rates[s,] = curr.mix
    
    # calculate summary spatial effects
    #post.direct[,s] = curr.ai_diag/n * curr.beta
    #post.total[,s] = curr.ai_tot/n * curr.beta
    #post.indirect[,s] = post.total[,s] - post.direct[,s]
    #AI = solve(diag(n) - curr.rho * W)
    # for (rr in 1:k) {
    #   SW = AI %*% (diag(n) * curr.beta[rr])
    #   post.direct[rr,s] = sum(diag(SW))/n
    #   post.total[rr,s] = sum(SW)/n
    #   post.indirect[rr,s] = post.total[rr,s] - post.direct[rr,s]
    # }
  }
}
# yhat
thinning = 10
thin_chain = seq(1,nretain,by = thinning)
posty = matrix(0,n,length(thin_chain))
cat("posterior calculations\n")
for (ss in 1:length(thin_chain)) {
  cat(ss,"from",length(thin_chain),"\n")
  s = thin_chain[ss]
  curr.xb = X %*% postb[,s]
  curr.ee = rnorm(n,0,sqrt(posts[s]))
  for (tt in 1:T) {
    from_to = ((tt-1)*smalln+1):(tt*smalln)
    curr.W = W.list[[acc.rates[s,tt]]]
    posty[from_to,ss] = as.matrix(
      solve(.sparseDiagonal(smalln) - postr[tt,s] * curr.W) %*% (curr.xb[from_to] + curr.ee[from_to]))
  }
}
cat("Done!\n")
yhat = apply(posty,c(1),mean)
# 
# # fit plots
# pdf("./Estimation/fit_plots.pdf",onefile = TRUE)
# par(mfrow=c(4,4))
# for (ccc in 1:length(shp@data$ISO3)) {
#   ccountry = shp@data$ISO3[ccc]
#   Y_ = exp(Y[countries == ccountry])
#   yhat_ = exp(yhat[countries == ccountry])
#   plot(Y_,  main=ccountry,
#        xlab="Time", ylab="Cases",
#        ylim = c(0,max(c(Y_,yhat_))),
#        type="l")
#   lines(yhat_,col="red")
# }
# par(mfrow=c(1,1))
# dev.off()
# 
# r-squared, rbar-squared
ym = Y - mean(Y)
rsqr1 = crossprod(Y - yhat )
rsqr2 = t(ym) %*% ym;
R2 = 1.0-rsqr1/rsqr2;   # r-squared
rsqr1 = rsqr1/(n-k - 2*T);
rsqr2 = rsqr2/(n-1.0);
R2bar = 1 - (rsqr1/rsqr2)
# 
# 
### output table of beta, R^2, sigma^2, sample size
require(openxlsx)
signif = function(x,y) {sign(x) == sign(y)}
rres = rbind(postb,posts)
rres = t(apply(rres,c(1),function (x) c(mean(x),sd(x),quantile(x,c(.025,.975)) )))
rres = cbind(rres,signif(rres[,ncol(rres)],rres[,ncol(rres)-1]))
colnames(rres) = c("mean","sd","p025","p975","significant")
row.names(rres) = c("alpha",substr(colnames(country_dummies),10,13),"sigma^2")
# add statistics about rho
postr_janfeb = postr[1:38,]
postr_march = postr[39:T,]
rres = rbind(rres,
             rho = c(mean(postr),sd(postr),quantile(postr,c(.025,.975)),
                     signif(quantile(postr,0.025),quantile(postr,.975))),
             rho1 = c(mean(postr_janfeb),sd(postr_janfeb),quantile(postr_janfeb,c(.025,.975)),
                     signif(quantile(postr_janfeb,0.025),quantile(postr_janfeb,.975))),
             rho2 = c(mean(postr_march),sd(postr_march),quantile(postr_march,c(.025,.975)),
                     signif(quantile(postr_march,0.025),quantile(postr_march,.975))))
# add other measures
rres = rbind(rres,
             R2bar = c(R2bar,NA,NA,NA,NA),
             N = c(smalln,NA,NA,NA,NA),
             T = c(T,NA,NA,NA,NA))
# no country fixed FX
rres = rres[-which(row.names(rres) %in% substr(colnames(country_dummies),10,13)),]
# # save data
# wb = loadWorkbook("./Latex/results_SAR4.xlsx")
# writeData(wb,1,rres, startCol = 3, startRow = 3,rowNames = FALSE)
# saveWorkbook(wb,"./Latex/results_SAR4.xlsx",overwrite = TRUE)
### done
# 
# 
# ### compute time specific direct/indirect FX  ###
# # thinning = 100
# # thin_chain = seq(1,nretain,by = thinning)
# # 
# # nhor = 14
# # post.direct = array(0,c(length(thin_chain),T,k))
# # post.indirect = array(0,c(length(thin_chain),T,k))
# # post.total = array(0,c(length(thin_chain),T,k))
# # 
# # post.fcast = array(0,c(length(thin_chain),T,smalln,nhor))
# # 
# # require(progress)
# # pb <- progress_bar$new(
# #   format = "  [:bar] :percent in :eta",
# #   total = length(thin_chain), clear = FALSE, width= 60)
# # 
# # for (ss in 1:length(thin_chain)) {
# #   s = thin_chain[ss]
# #   for (tt in 1:T) {
# #     from_to = ((tt-1)*smalln+1):(tt*smalln)
# #     
# #     curr.W = W.list[[acc.rates[s,tt]]]
# #     curr.AI = solve(.sparseDiagonal(smalln) - postr[tt,s] * curr.W)
# #     
# #     for (kk in 1:k) {
# #       SW = curr.AI %*% (diag(smalln) * postb[kk,s])
# #       post.direct[ss,tt,kk] = sum(diag(SW))/smalln
# #       post.total[ss,tt,kk] = sum(SW)/smalln
# #       post.indirect[ss,tt,kk] = post.total[ss,tt,kk] - post.direct[ss,tt,kk]
# #     }
# #     
# #     yhat = X[from_to,]
# #     for (hh in 1:nhor) {
# #       yhat[,1] = as.vector(curr.AI %*% (yhat %*% postb[,s]))
# #       post.fcast[ss,tt,,hh] = exp(yhat[,1])
# #     }
# #   }
# #   pb$tick()
# # }
# # 
# # post.fcast.mean = apply(post.fcast,c(2,3,4),mean)
# # post.fcast.growth_rt = post.fcast.mean
# # for (tt in 1:T) {
# #   post.fcast.growth_rt[tt,,] = post.fcast.mean[tt,,] / matrix(exp(Xraw[,tt]),smalln,nhor)
# # }
# # 
# # avg_fcast = (apply(post.fcast,c(3,4),mean))
# # avg_X = matrix(apply(exp(Xraw),c(1),mean),smalln,nhor)
# #############
# 
# 
# ##### - look at fcast after 12/03/2020
# 
# upto_date = which(dates == "2020-03-06")
# from_date = which(dates == "2020-03-12")
# upto_W.prob = table(acc.rates[,1:upto_date]) / (nretain * upto_date)
# fcast_T = from_date:T
# nhor = length(fcast_T)
# 
# thinning = 10
# thin_chain = seq(1,nretain,by = thinning)
# 
# post.fcast = array(0,c(length(thin_chain),smalln,nhor))
# post.fcast2 = array(0,c(length(thin_chain),smalln,nhor))
# diff.fcast = array(0,c(length(thin_chain),smalln,nhor))
# 
# require(progress)
# pb <- progress_bar$new(
#   format = "  [:bar] :percent in :eta",
#   total = length(thin_chain), clear = FALSE, width= 60)
# for (ss in 1:length(thin_chain)) {
#   s = thin_chain[ss]
#   
#   start_from_to = ((from_date-2)*smalln+1):((from_date-1)*smalln)
#   yhat = yhat2 = X[start_from_to,]
#   for (hh in 1:nhor) {
#     tt = fcast_T[hh]
#     curr.W = W.list[[sample(1:4,1,prob = upto_W.prob)]]
#     curr.AI = solve(.sparseDiagonal(smalln) - postr[tt,s] * curr.W)
# 
#     yhat[,1] = as.vector(curr.AI %*% (yhat %*% postb[,s]))
#     post.fcast[ss,,hh] = exp(yhat[,1]) #* wbpop$value
#     
#     ### counterfactural with previous rho
#     upto_tt = sample(1:upto_date,1)
#     curr.W2 = W.list[[acc.rates[s,upto_tt]]]
#     curr.AI2 = solve(.sparseDiagonal(smalln) - mean(postr[1:upto_date,s]) * curr.W2)
#     
#     yhat2[,1] = as.vector(curr.AI2 %*% (yhat2 %*% postb[,s]))
#     post.fcast2[ss,,hh] = exp(yhat2[,1]) #* wbpop$value
#     
#     diff.fcast[ss,,hh] = exp(yhat2[,1]) - exp(yhat[,1])
#   }
#   pb$tick()
# }
# # output graphics and figures
# pdf("./Latex/diff_rho_plots.pdf",onefile = TRUE)
# par(mfrow=c(4,4))
# for (ccc in 1:length(shp@data$ISO3)) {
#   ccountry = shp@data$ISO3[ccc]
#   to_plot=data.frame(
#     med1 = apply(post.fcast[,ccc,],c(2),median),
#     low1 = apply(post.fcast[,ccc,],c(2),quantile,.16),
#     high1 = apply(post.fcast[,ccc,],c(2),quantile,.84),
#     med2 = apply(post.fcast2[,ccc,],c(2),median),
#     low2 = apply(post.fcast2[,ccc,],c(2),quantile,.16),
#     high2 = apply(post.fcast2[,ccc,],c(2),quantile,.84)
#   )
#   
#   plot(to_plot$med1,  main=ccountry,
#        xlab="Days", ylab="Cases",
#        ylim = c(0,max(to_plot[,c("med1","med2")])),
#        type="l")
#   lines(to_plot$low1,lty=2)
#   lines(to_plot$high1,lty=2)
#   lines(to_plot$med2,col="red")
#   lines(to_plot$low2,lty=2,col="red")
#   lines(to_plot$high2,lty=2,col="red")
# }
# par(mfrow=c(1,1))
# dev.off()
# 
# ### visualize selected countries
# sel_countries = c("AUT")
# dates <- seq(as.Date("2020-03-10"), length.out = nhor, by="days")
# for (jj in 1:length(sel_countries)) {
#   country = sel_countries[jj]
#   ccc = which(shp@data$ISO3 == country)
#   
#   to_plot=data.frame(
#     days = dates,
#     med3 = apply(diff.fcast[,ccc,],c(2),median),
#     mea3 = apply(diff.fcast[,ccc,],c(2),mean),
#     low3 = apply(diff.fcast[,ccc,],c(2),quantile,.05),
#     high3 = apply(diff.fcast[,ccc,],c(2),quantile,.95),
#     lowest3 = apply(diff.fcast[,ccc,],c(2),quantile,.16),
#     highest3 = apply(diff.fcast[,ccc,],c(2),quantile,.84),
#     sd3 = apply(diff.fcast[,ccc,],c(2),sd)
#   )
#   to_plot = to_plot[1:16,]
#   #to_plot[,-1] = to_plot[,-1] / 1000
#   
#   textsize = 14
#   p = ggplot(to_plot,aes(days,med3)) + geom_line() +
#     theme_minimal() + xlab("") + ylab("Additional infections \nwithout border closure policies") +
#     coord_cartesian(ylim=c(0,floor10(max(to_plot$med3) * 3))) +
#     #coord_cartesian(ylim=c(0,2000)) +
#     ggtitle(countrycode(country,"iso3c","country.name")) +
#     geom_ribbon(aes(ymin=low3,ymax = high3), alpha = .2) +
#     geom_ribbon(aes(ymin=lowest3,ymax = highest3),alpha = .15) +
#     scale_colour_manual("",values="black")+
#     scale_fill_manual("",values=c("grey12", "darkgrey")) +
#     scale_x_date(date_breaks = "3 days", date_labels =  "%b %d") +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 75 , hjust = 1, size = textsize),
#           axis.text.y = element_text(size = textsize, angle = 15),
#           axis.title = element_text(size = textsize -1),
#           axis.title.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           panel.grid.minor.x = element_blank(),
#           panel.grid.major.x = element_line(colour = 'white', size = 0.4),
#           panel.grid.minor.y = element_blank(),
#           panel.grid.major.y = element_line(colour = 'lightgrey', size = 0.4),
#           #axis.title.y = element_text(hj),
#           legend.title = element_text(face = 'bold', size = textsize + 1),
#           legend.position = 'bottom',
#           legend.text = element_text(size = textsize))
#   
#   ggsave(paste0("./Latex/final_country_bw_",country,".pdf"), 
#          p, width = 7, height = 4.25)
#   ggsave(paste0("./Latex/final_country_bw_",country,".jpg"), 
#          p, width = 7, height = 4.25)
# }
# 
# ############
# 
acc.rates_old = acc.rates
res = matrix(0,ddd,T)
for (d in 1:ddd) {
  res[d,] = colSums(acc.rates == d)  / nretain
}
dates <- seq(as.Date("2020-01-23"), length.out = T, by="days")
dates2 <- seq(as.Date("2020-01-28"), length.out = (T-5) , by="days")
# res <- t(res/nsave)
res <- cbind(dates, data.frame(t(res)))

res2 <- cbind(dates2, data.frame(apply(res[,2:ncol(res)],2,rollmean,k=6, align = "right")))
res2[,2:ncol(res2)] <- res2[,2:ncol(res2)]/rowSums(res2[,2:ncol(res2)])

wmats = names(W.list)
colnames(res2) <- c("dates",wmats)
res2 <- gather(res2, "Wmat", "Prob" , -dates)


j <- ggplot(res2, aes(dates, Prob))
j = j + geom_line(aes(color=Wmat), size=1) +
  theme_minimal() +
  xlab("Dates") + ylab("Probability of W") +
  ylim(0,1) +
  scale_x_date(date_breaks = "7 days", date_labels =  "%b %d") #+
# scale_color_discrete(name="Spatial specification", labels=c("Common Currency","Flight intensity","FTA","Language ties", "Common Borders"))


rho_res = apply(postr,c(1),quantile,c(.025,.16,.5,.84,.975))
rho_res = data.frame(dates = dates,t(rho_res))
colnames(rho_res)[-1] = c("lowest","low","mid","high","highest")
rho_range = ceiling(max(abs(rho_res[,-1]))*10) / 10
p = ggplot(rho_res,aes(dates,mid)) + geom_line() +
      theme_minimal() + xlab("Dates") + ylab("Rho") +
      coord_cartesian(ylim=c(-rho_range,rho_range)) +
    geom_ribbon(aes(ymin=low,ymax = high), alpha = .2,fill = "blue") +
    geom_ribbon(aes(ymin=lowest,ymax = highest),alpha = .1,fill = "blue") +
   geom_hline(yintercept = 0,alpha = .5) +
  scale_x_date(date_breaks = "7 days", date_labels =  "%b %d")
# 
# ggsave("./Latex/Wspecfig_SAR.pdf",j,width = 7, height = 5)
# ggsave("./Latex/rho_plot.pdf",p,width = 7, height = 5)
# 
 # save.image(paste0("./Estimation/SAR4",Sys.Date(),".RData"))
 # 
 # load(paste0("./Estimation/SAR42020-06-29.RData"))

 # #### Mcfadden
 # llh.thin <- 10
 # llh.seq <- seq(1,nsave,llh.thin)
 # 
 # llh <- rep(NA, length(llh.seq))
 # 
 # a <- 1
 # for(i in llh.seq){
 #   
 #   
 #   
 #   err <- Y[i]-posty[a]
 #   
 #   ldet <- 0
 #   j <- 1
 #   for(j in 1:T){
 #     curr.W.mat <- acc.rates[i,j]
 #     ldet <- ldet + log(det(diag(N)-postr[i]*W.list[[curr.W.mat]]))
 #   }
 #   llh.2 = - (N*T)/2 * log(pi*posts[i])
 #   llh.3 = - t(err) %*% err / (2*posts[i])
 #   llh[a]= ldet + llh.2 + llh.3
 #   a <- a+1
 # }
 # llh.sof <- mean(llh)
 # 
 
 
 X.ols <- rep(1,length(Y))
 beta.ols <- solve(t(X.ols) %*% X.ols) %*% (t(X.ols) %*% Y)
 resi <- Y-X.ols%*%beta.ols
 
 se2 <- sum(resi ^ 2) / (N*T - 1)
 
 llh.ols <- sum(dnorm(Y, mean=X.ols %*% beta.ols, sd=sqrt(se2), log = TRUE))
 
 
 McFadden <- 1-(llh.sof/llh.ols)


