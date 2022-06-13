###############################################
################ Load libraries ###############
###############################################
library("scales")
library("KFAS")
library("latex2exp")
set.seed(12) #12 #123 #1234 #234 #2345

###############################################
######### Measurement noise scenarios #########
###############################################
iid <- function(n, se, lambda, t, speed) {
  # n : scalar
  # se : vector 
  # lambda : vector
  # t : vector
  # speed : scalar
  
  p <- ncol(se) # number of measurement signals
  
  es <- matrix(nrow=n, ncol=p) # initial noise matrix
  for (i in seq(1,n)) { # iid noise
    es[i,1] <- rnorm(1, sd=se[1])
    es[i,2] <- rnorm(1, sd=se[2])
    es[i,3] <- rnorm(1, sd=se[3])
  }
  
  as <- cos(t) + speed  # true values vector
  ys <- rbind(lambda)%x%cbind(as) + es # measurements vector
  
  return(list("as"=as,"ys"=ys,"es"=es))
}

AR1 <- function(n, se, theta, lambda, t, speed) {
  # n : scalar
  # se : vector 
  # theta : vector
  # lambda : vector
  # t : vector
  # speed : scalar
  
  p <- ncol(se) # number of measurement signals
  
  es <- matrix(nrow=n, ncol=p) # initial noise matrix
  es[1,] <- rnorm(p) # Initial conditions of the first row
  
  for (i in seq(2,n)) { # AR1 noise
    es[i,1] <- theta[1]*es[i-1,1] + rnorm(1, sd=se[1])
    es[i,2] <- theta[2]*es[i-1,2] + rnorm(1, sd=se[2])
    es[i,3] <- theta[3]*es[i-1,3] + rnorm(1, sd=se[3])
  }
  
  as <- cos(t) + speed  # true values vector
  ys <- rbind(lambda)%x%cbind(as) + es # measurements vector
  
  return(list("as"=as,"ys"=ys,"es"=es))
}

nse <- function(n, se, lambda, t, speed) {
  # n : scalar
  # se : vector 
  # lambda : vector
  # t : vector
  # speed : scalar
  
  p <- ncol(se) # number of measurement signals
  
  es <- matrix(nrow=n, ncol=p) # initial noise matrix
  for (i in seq(1,n)) { # iid noise
    es[i,1] <- rnorm(1, sd=se[1])
    es[i,2] <- rnorm(1, sd=se[2])
    es[i,3] <- rnorm(1, sd=se[3])
  }
  
  as <- cos(t) + speed  # true values vector
  ys <- rbind(lambda)%x%cbind(as) + es # measurements vector
  
  ys[500:600,1] <- sin(ys[500:600,1]) + speed + rnorm(1, sd=se[1])
  ys[1300:1700,2] <- sin(ys[1300:1700,2]) + speed + 2
  ys[1100:1300,3] <- (ys[1100:1300,3]) - 6
  ys[1300:1500,3] <- (ys[1300:1500,3]) - 9
  ys[1500:1800,3] <- (ys[1500:1800,3]) + 3
  
  return(list("as"=as,"ys"=ys,"es"=es))
}

###############################################
############ Dynamic Factor Models ############
###############################################

ssm1 <- function(ys){
  ##### Model Design ######
  Z <- rbind(c(1),            # H matrix in paper's notation
             c(NA),
             c(NA))
  p <- ncol(Z)
  q <- nrow(Z)
  H <- matrix(0,q,q)          # R matrix in paper's notation

  Q <- cbind(rep(NA,p))       # Q matrix in paper's notation
  T <- cbind(rep(1,p))        # F matrix in paper's notation
  a0 <- cbind(rep(0,p))
  P0 <- cbind(rep(1e3,p))

  ############ Set the param optimization function ###########
  objf <- function(pars, model, estimate = TRUE) {
    model$Q[,,1] <- exp(pars[1])
    model$H[,,1] <- diag(exp(pars[2:4]))
    ##browser()
    model$Z[,,1] <- c(1,pars[5:6])
    if (estimate) {
      -logLik(model)
    } else {
      model
    }
  }

  ############### Define parameter number ###################
  p0 <- rep(.5,6)

  ########## Set the state space model for each leg #########
  ssm <- SSModel(ys ~ -1 + SSMcustom(Z = Z, T = T, Q = Q, a1 = a0, P1 = P0), H=H)

  ################## Parameter optimization #################
  opt <- optim(p0, objf, model=ssm, method="BFGS") # BFGS, Nelder-Mead, CG, L-BFGS-B, SANN, Brent
  fit <- objf(opt$par, model=ssm, estimate=FALSE)

  ################## Estimate ###############
  out <- KFS(fit, smoothing="mean")
} # 1-factor model design

ssm2 <- function(ys){
  ##### Model Design ######
  Z <- rbind(c(1,0,0),           # H matrix in paper's notation
             c(NA,1,0),
             c(NA,0,1))
  p <- ncol(Z)
  q <- nrow(Z)
  H <- matrix(0,q,q)             # R matrix in paper's notation
  
  Q <- diag(rep(NA,p))*1         # Q matrix in paper's notation
  T <- diag(c(1,rep(NA,p-1)))    # F matrix in paper's notation
  a0 <- cbind(rep(0,p))
  P0 <- diag(rep(1e3,p))
  
  ############ Set the param optimization function ###########
  objf <- function(pars, model, estimate = TRUE) {
    model$Q[,,1] <- diag(exp(pars[1:3]))
    model$H[,,1] <- diag(exp(pars[4:6]))
    model$T[,,1] <- diag(c(1,pars[7:8]))
    
    Z <- model$Z[,,1]
    Z[seq(q-1)+1,1] <- pars[9:10]
    model$Z[,,1] <- Z
    if (estimate) {
      -logLik(model)
    } else {
      model
    }
  }
  
  ############### Define parameter number ###################
  p0 <- rep(.5,10)
  
  ########## Set the state space model for each leg #########
  ssm <- SSModel(ys ~ -1 + SSMcustom(Z = Z, T = T, Q = Q, a1 = a0, P1 = P0), H=H)
  
  ################## Parameter optimization #################
  opt <- optim(p0, objf, model=ssm, method="BFGS") # BFGS, Nelder-Mead, CG, L-BFGS-B, SANN, Brent 
  fit <- objf(opt$par, model=ssm, estimate=FALSE)
  
  ################## Estimate ###############
  out <- KFS(fit, smoothing="mean")
} # 3-factors model design

ssm3 <- function(ys){
  ##### Model Design ######
  Z <- rbind(c(1,0,0,1,0,0),    # H matrix in paper's notation
             c(NA,1,0,0,1,0),
             c(NA,0,1,0,0,1))
  p <- ncol(Z)
  q <- nrow(Z)
  H <- matrix(0,q,q)            # R matrix in paper's notation
  
  Q <- diag(rep(NA,p))*1        # Q matrix in paper's notation
  T <- diag(p)                  # F matrix in paper's notation
  a0 <- cbind(rep(0,p))
  P0 <- diag(rep(1e3,p))
  
  ############ Set the param optimization function ###########
  objf <- function(pars, model, estimate = TRUE) {
    model$T[,,1] <- diag(c(1,tanh(pars[1:2]),exp(pars[3:5])))
    model$Q[,,1] <- diag(exp(pars[6:11]))
    
    Z <- model$Z[,,1]
    Z[seq(q-1)+1,1] <- pars[12:13]
    model$Z[,,1] <- Z
    if (estimate) {
      -logLik(model)
    } else {
      model
    }
  }
  
  ############### Define parameter number ###################
  p0 <- rep(0.05,13)
  
  ########## Set the state space model for each leg #########
  ssm <- SSModel(ys ~ -1 + SSMcustom(Z = Z, T = T, Q = Q, a1 = a0, P1 = P0), H=H)
  
  ################## Parameter optimization #################
  opt <- optim(p0, objf, model=ssm, method="BFGS") # BFGS, Nelder-Mead, CG, L-BFGS-B, SANN, Brent 
  fit <- objf(opt$par, model=ssm, estimate=FALSE)
  
  ################## Estimate ###############
  out <- KFS(fit, smoothing="mean")
  
  return(out)
} # 3-factors model design with AR(1) meas. noise   

###############################################
############## Plotting function ##############
###############################################
plot.results <- function(res,title, estimate=FALSE, out=NA, model_name=NA, relativeRMSE=FALSE) {
  # Change background 
  par(bg="white")
  
  # Define size
  options(repr.plot.width=18, repr.plot.height=12)
  
  # Configure each graph's size
  par(mar = c(5.1, 5.1, 4.1, 2.1),
      cex.main=1.2, cex.lab=1.2, cex.axis=1.2, cex.sub=1.2)
  
  plot(res$ys[,1], col=alpha('red', 0.4), pch=16, cex=.5, 
       xlab="Time [10 mins]", ylab="Speed [kn]",
       ylim=c(0,24), yaxp=c(0, 24, 24),
       xlim=c(0,2000), xaxp=c(0, 2000, 8))#
  points(res$ys[,2], col=alpha('green', 0.4), pch=16, cex=.5)
  points(res$ys[,3], col=alpha('blue', 0.4), pch=16, cex=.5)
  lines(res$as, col="red", lwd=2, lty=2)
  # grid(7, 13, lwd = 1) # grid only in y-direction
  grid()

  if (estimate) {
    # Calculate RMSE
    residuals <- cbind(res$as) - cbind(out[2:length(out)])
    rmse = sqrt(colMeans(residuals^2))
    
    lines(out, col="purple", lwd=2)   # Plot estimate
    legend(x = 1850, y = 6, legend=c(expression(U[nu*c]), expression(U[nu*e]), expression(U[dvl]), expression(alpha), expression(widehat(alpha))),
           bg = "transparent", box.col = "transparent", col=c("red", "green", "blue","red", "purple"), pch=c(16,16,16,NA,NA),lty=c(NA,NA,NA,2,1),cex=1.2)
    text(1000, 3,  bquote("RMSE " ~ widehat(alpha) ~ ": " ~ .(round(rmse,4)) ~ "[kn]"), cex=1.15, pos=3,col="purple")
    
    title(main=paste(title, ' - ',model_name)) # Change title
  } else {
    # Calculate all residuals
    res_Uvc <- cbind(res$as) - cbind(res$ys[,1])
    res_Uve <- cbind(res$as) - cbind(res$ys[,2])
    res_Udvl <- cbind(res$as) - cbind(res$ys[,3])
    res_ybar <- cbind(res$as) - cbind((res$ys[,1]+res$ys[,2]+res$ys[,3])/3)
    
    # Create a ColSd function equal to ColMeans
    colSd <- function(x)sqrt(rowMeans((t(x)-colMeans(x))^2)*((dim(x)[1])/(dim(x)[1]-1)))
    
    if (relativeRMSE){
      # Calculate relative rmse
      rmse_Uvc <- sqrt(colMeans(res_Uvc^2))/colSd(cbind(res$as))
      rmse_Uve <- sqrt(colMeans(res_Uve^2))/colSd(cbind(res$as))
      rmse_Udvl <- sqrt(colMeans(res_Udvl^2))/colSd(cbind(res$as))
      rmse_ybar <- sqrt(colMeans(res_ybar^2))/colSd(cbind(res$as))
    } else {
      # Calculate rmse
      rmse_Uvc <- sqrt(colMeans(res_Uvc^2))
      rmse_Uve <- sqrt(colMeans(res_Uve^2))
      rmse_Udvl <- sqrt(colMeans(res_Udvl^2))
      rmse_ybar <- sqrt(colMeans(res_ybar^2))
    }
    
    lines((res$ys[,1]+res$ys[,2]+res$ys[,3])/3, col=alpha('pink', 0.4), lwd=2)   # Plot estimate
    legend(x = 1850, y = 12, legend=c(expression(U[nu*c]), expression(U[nu*e]), expression(U[dvl]), expression(alpha), expression(bar(bold(y)))),
           bg = "transparent", box.col = "transparent", col=c("red", "green", "blue","red", "pink"), pch=c(16,16,16,NA,NA),lty=c(NA,NA,NA,2,1),cex=1.2)
    text(1000, 5.5,  bquote("RMSE " ~ U[nu*c] ~ ": " ~ .(round(rmse_Uvc,4)) ~ "[kn]"), cex=1.15, pos=3,col="red") 
    text(1000, 4,  bquote("RMSE " ~ U[nu*e] ~ ": " ~ .(round(rmse_Uve,4)) ~ "[kn]"), cex=1.15, pos=3,col="green")
    text(1000, 2.5,  bquote("RMSE " ~ U[dvl] ~ ": " ~ .(round(rmse_Udvl,4)) ~ "[kn]"), cex=1.15, pos=3,col="blue")
    text(1000, 1,  bquote("RMSE " ~ bar(bold(y)) ~ ": " ~ .(round(rmse_ybar,4)) ~ "[kn]"), cex=1.15, pos=3,col="pink")
    
    title(main=paste(title)) # Set title
  }
}

###############################################
######## Data and estimate generation #########
###############################################
# Generate a simulation with iid noise
iid_res <- iid(n=2e3, se=cbind(1,1,1), lambda=rep(1, 3),  # Change to se=cbind(0.5,1,2) for real data noise scenario
               t=seq(0, 6*pi, length.out=2e3), speed=15)
iid_ssm1_out <- ssm1(iid_res$ys)   # 1-factor model
iid_ssm2_out <- ssm2(iid_res$ys)   # 3-factors model
iid_ssm3_out <- ssm3(iid_res$ys)   # 3-factors model with AR(1) meas. noise  


# Generate a simulation with AR1 noise
AR1_res <- AR1(n=2e3, se=cbind(1,1,1),                    # Change to se=cbind(0.5,1,2) for real data noise scenario
               theta=rep(.1, 3), lambda=rep(1, 3),       
               t=seq(0, 6*pi, length.out=2e3), speed=15)
AR1_ssm1_out <- ssm1(AR1_res$ys)    # 1-factor model
AR1_ssm2_out <- ssm2(AR1_res$ys)    # 3-factors model 
AR1_ssm3_out <- ssm3(AR1_res$ys)    # 3-factors model with AR(1) meas. noise 


# Generate a simulation with non-stationary measurement noise
nse_res <- nse(n=2e3, se=cbind(1,1,1), lambda=rep(1, 3),  # Change to se=cbind(0.5,1,2) for real data noise scenario
               t=seq(0, 6*pi, length.out=2e3), speed=15)
nse_ssm1_out <- ssm1(nse_res$ys)    # 1-factor model
nse_ssm2_out <- ssm2(nse_res$ys)    # 3-factors model
nse_ssm3_out <- ssm3(nse_res$ys)    # 3-factors model with AR(1) meas. noise


###############################################
############## Plotting results ###############
###############################################
# Generate iid noise scenario results from all models
plot.results(iid_res,"Independent and identically distributed (i.i.d.) noise") #relativeRMSE=TRUE
plot.results(iid_res,"Independent and identically distributed (i.i.d.) noise",estimate=TRUE,iid_ssm1_out$a[,1],'SSM #1')
plot.results(iid_res,"Independent and identically distributed (i.i.d.) noise",estimate=TRUE,iid_ssm2_out$a[,1],'SSM #2') 
plot.results(iid_res,"Independent and identically distributed (i.i.d.) noise",estimate=TRUE,iid_ssm3_out$a[,1],'SSM #3') 

# Generate AR1 noise scenario results from all models
plot.results(AR1_res,"AR(1) noise")
plot.results(AR1_res,"AR(1) noise",estimate=TRUE,AR1_ssm1_out$a[,1],'SSM #1')
plot.results(AR1_res,"AR(1) noise",estimate=TRUE,AR1_ssm2_out$a[,1],'SSM #2')
plot.results(AR1_res,"AR(1) noise",estimate=TRUE,AR1_ssm3_out$a[,1],'SSM #3') 

# Generate non stationary noise scenario results from all models
plot.results(nse_res,"Non stationary (n.s.) noise")
plot.results(nse_res,"Non stationary (n.s.) noise",estimate=TRUE,nse_ssm1_out$a[,1],'SSM #1')
plot.results(nse_res,"Non stationary (n.s.) noise",estimate=TRUE,nse_ssm2_out$a[,1],'SSM #2') 
plot.results(nse_res,"Non stationary (n.s.) noise",estimate=TRUE,nse_ssm3_out$a[,1],'SSM #3')


