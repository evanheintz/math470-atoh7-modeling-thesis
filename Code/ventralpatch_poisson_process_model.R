#----------------------------------------------------------------------------
##############################################################################
#@@@@@@@@@@@@@@@ Ventral Patch Poisson Process @@@@@@@@@@@@@@@#
# author: Evan Heintz

# description: defines log likelihood and fits updated ventral patch poisson 
# point process model to expanded sweetrain dataset.

# requirements: 'data_expanded_final_new' file from 'thesis_sweetrain.R' script.

# outputs: MLE estimated parameters for updated poisson process model, simulated
# data and plots from the model's intensity function.
##############################################################################
#----------------------------------------------------------------------------

#------ load libraries ------#
packages <- c("stats","tidyverse","tiff","ggforce","gridExtra","viridis","ggthemes")
lapply(packages, library, character.only = TRUE)
theme_set(theme_minimal())



#log likelihood (with separate ventral patch component)
logl <- function(par, data){
  sigma12 <- par[1]
  sigma22 <- par[2]
  h <- par[3]
  d <- par[4]
  e <- par[5]
  f <- par[6]
  k <- par[7]
  p <- par[8]
  m <- par[9]
  tau1 <- par[10]
  tau2 <- par[11]
  b <- par[12]
  c <- par[13]
  #new parameters
  int <- par[14] #intensity of ventral patch expression
  tau3 <- par[15] #time for which expression moves outside ventral patch
  major <- par[16] #minor coord for ellipse of ventral patch
  minor <- par[17] #major coord for ellipse of ventral patch
  xcent <- par[18] #x center for ellipse of ventral patch
  ycent <- par[19] #y center for ellipse of ventral patch
  beta0 <- par[20] #overall intercept 
  #define parameters
  
  a = 2*pi/(tau1-tau3)
  #a is defined in terms of tau1

  ventral_time <- data %>% filter(t <= tau3)
  no_ventral_time <- data %>% filter(t > tau3) #set to tau3 parameter
  
  xoft <- d+b*cos(a*no_ventral_time$t + f)
  yoft <- e+c*sin(a*no_ventral_time$t + f)
  
  central_indicator <- ((ventral_time$x - xcent)^2)/(major^2) + ((ventral_time$y - ycent)^2)/(minor^2) <= 1
  
  #if points within ventral patch ellipse
  ventral <- ifelse(central_indicator, int, 0) 
  
  #ventral patch, if inside ellipse equation from parameters, then will have homogeneous
  #point process intensity for expression
  ventral <- c(ventral, rep(0,length(no_ventral_time$t)))
  
  lambdas_rotational <- c(rep(0,length(ventral_time$t)),
                          (no_ventral_time$t <= tau1) * (h / (2 * pi * sigma12)) * exp( - ((no_ventral_time$x - xoft)^2 + (no_ventral_time$y - yoft)^2) / (2 * sigma12) ))
  #rotational intensity function component
  
  lambdas_radial <- c(rep(0,length(ventral_time$t)), (no_ventral_time$t <= tau2) * (p / sqrt(2 * pi * sigma22)) * exp( - (sqrt((no_ventral_time$x - d)^2 + (no_ventral_time$y - e)^2/m^2) - k*no_ventral_time$t)^2 / (2 * sigma22)))
  #radial intensity function component
  
  lambda <- lambdas_rotational + lambdas_radial + beta0 + ventral
  #combined by adding both components together
  
  if (h <= 0) {
   # print("h is the problem")
    return(99999)
  }
  if (p <= 0){
   # print("p is the problem")
    return(99999)
  } 
  if (m <= 0){
   # print("m is the problem")
    return(99999)
  } 
  
  if (sigma12 <= 0){
   # print("sigma12 is the problem")
    return(99999)
  }
  if (sigma22 <= 0){
   # print("sigma22 is the problem")
    return(99999)
  }
  #restrict f to less than the period of sin(at)
   if (f > 2*pi/a){
   # print("f is the problem")
    return(99999)
  } 
  if (f < 0){
   # print("negative f is the problem")
    return(99999)
  } 
  
  if (tau1 <= 0 | tau1 >= 77){
    #print("tau1 is the problem")
    return(99999)
  } 
  
  if (tau2 <= 0 | tau2 >= 77){
    #print("tau2 is the problem")
    return(99999)
  } 
  
  if (tau3 <= 1 | tau3 >= 10 | tau3 >= tau1 | tau3 >= tau2){
   # print("probably not ventral patch anymore")
    return(99999)
  } 
  
  if (minor <= 0 | major <= 0 | xcent <= 0 | ycent <= 0){
    #print("problems")
    return(99999) 
  } 

   if (beta0 <= 0){
    #print("beta0 is the problem")
    return(99999)
   }
  
  if (sigma12 > 5000 | sigma22 > 5000){
    #print("sigma12 or sigma22 is the problem")
    return(99999)
  }
  
  
  if (sum(is.nan(lambda)) > 0) cat(a, d, e, f, sigma12, sigma22, m, h, p, k)
  
  #integral of intensity function
  # intlam <- # pi*major*minor*int*tau3 + #ventral patch integral
   intlam <- tau3*pi*major*minor*int + h*(tau1 - tau3) + pi*m*p*k*(tau2-tau3)^2 + # (h*tau3) +  #+ pi*m*p*k*tau3^2
    beta0*max(data$t)*pi*d*e # coordinates of center of retina
  # using hull_df for respective video
  #not ventral patch integral from tau3 to max T
  
  
  loglik = sum(log(lambda)) - intlam
  #log likelihood
  
  return(-1.0*loglik)
  #since optim() is set to minimize, we multiply by -1 to ensure that we
  #are getting the MAXIMUM log likelihood
}

#parameter optimization
set.seed(343)
m3 = function(x) signif(x,3)
#number of iterations
tries = 10
#values from math 343 report for optimized parameters and educated guesses for new parameters


par0 <- c(
  4060, # sigma12
  1772, # sigma22
  8.772, # h
  246.2, # d
  234, # e
  11.64, # f
  4.198, # k
  0.04108, # p
  0.7688, # m
  74.42, # tau1
  76.94, # tau2
  -134.4, # b
  48.1, # c
  0.0015, # int 
  5, # tau3
  50, # major
  50, # minor
  165, # xcent
  175, # ycent
  0.0000001 # beta0  
)


par_est <- as.data.frame(matrix(NA, nrow = tries, ncol = 20))

# manually assign each parameter a different sd and making sure
# no conditions are violated from 'logl' function
# also depending on how robust we think parameters should be
sd <- numeric()
sd[1] <- 100 # sigma12
sd[2] <- 100 # sigma22
sd[3] <- 0.05 # h
sd[4] <- 5 # d
sd[5] <- 5 # e
sd[6] <- 0.05 # f
sd[7] <- 0.25 # k
sd[8] <- 0.0025 # p
sd[9] <- 0.025 # m
sd[10] <- 2 # tau1
sd[11] <- 2 # tau2
sd[12] <- 5 # b
sd[13] <- 2 # c
sd[14] <- 0.0001 # int
sd[15] <- 1 # tau3
sd[16] <- 2 # major
sd[17] <- 2 # minor
sd[18] <- 5 # xcent
sd[19] <- 5 # ycent
sd[20] <-0.0000000001 # beta0  

st <- Sys.time()
for (i in 1:tries){
  set.seed(100*i)
  # sample from expanded data with size 5000 (can be altered)
  sample <- sample_n(sweetrain, 3000)
  # sample <- sample_n(data_expanded_final_new %>% select(x,y,t), 5000)
  
  #%%%%% can sample from expanded data, but apparently the 'logl' function
  # is really sensitive to starting values, so parameter estimates were
  # spiking off to infinity (particularly sigma12 and sigma22)
  
  start <- Sys.time()
  
  diff <- 100
  it <- 0
  
  theta1 <- par0
  
  # now reassign theta1 with slightly different random parameters
  # making sure conditions still hold from 'logl' function
  set.seed(100*i)
  
  for(j in 1:length(theta1)){
    theta1[j] <- rnorm(1,theta1[j],sd[j]) #scaled/similar to poisson mean/sd
    #to explore different regions of parameter space in optimization
    if(j != 12 & theta1[j] < 0){ #except for b, which is negative
      theta1[j] <- abs(theta1[j]) #make sure parameters are positive
    }
  }
  
  # check conditions
  if(theta1[10] >= 77){
    #tau1
    theta1[10] <- 74.42 + rnorm(1,0,0.5)
  }
  if(theta1[11] >= 77){
    #tau2
    theta1[11] <- 76.94 - abs(rnorm(1,0,0.5))
  }
  if(theta1[10] > theta1[11]){
    theta1[10] <- 74.42 + rnorm(1,0,0.5)
    theta1[11] <- 76.94 - abs(rnorm(1,0,0.5))
  }
  if(theta1[15] <= 1 | theta1[15] >= 10){
    #tau3
    theta1[15] <- 5 + sample(c(1,-1),1)
  } 
  
  while(diff > 1e-5){
    
    opt_out <- optim(par = theta1, fn = logl, data = sample)
    # optimize parameters from theta 1
    
    diff <- sum(abs(opt_out$par - theta1))
      
    
    it <- it + 1
    print(it)
    print(round(theta1, 2))
    print(diff)
    print("")
    
    if (it > 100) break
  } 
  
  par_est[i, ] <- m3(opt_out$par) #saving parameter estimates
  
  print(Sys.time() - start) 
  print(i)
}
print(Sys.tim
  #sample <- sample_n(data_expanded_final_new %>% select(x,y,t), 3000)
  sample <- sample_n(sweetrain,3000)
  
  start <- Sys.time()
  
  diff <- 100
  it <- 0
  
  theta_new <- c()
  for(j in 1:length(theta1)){
    sd <- abs(theta1[j])*sqrt(abs(theta1[j]))
    theta_new[j] <- rnorm(1,theta1[j],sd) #scaled/similar to poisson mean/sd
  }
  
  while (diff > 1e-5) {
    opt_out <- optim(par = theta_new, fn = logl, data = sample)
    
    diff <- sum(abs(opt_out$par - theta_new))
    # theta1 <- opt_out$par # optimized parameters from optim()
    
    it <- it + 1
    print(it)
    print(round(theta1, 2))
    print(diff)
    print("")
    
    # if (it > 100) break
  }
  
  par_est[i,] <- m3(opt_out$par) 
  print(Sys.time() - start)
  print(i)
  
  i <- i + 1
}



#%%%%%%%%%% saving parameter estimates and simulating data


#average of parameters from iterations
par_est_mean <- numeric()
for(i in 1:length(theta1)){
  par_est_mean[i] <- mean(par_est[,i],na.rm = TRUE)
}

#simulate data
set.seed(343)

compress_points <- function(x, y, d, e) {
  shrink <- sqrt((x - d)^2 / d^2 + (y - e)^2 / e^2)
  shrink[shrink < 1] <- 1 # points inside the ellipse remain
  x_compressed <- ((x - d) / shrink) + d
  y_compressed <- ((y - e) / shrink) + e
  data.frame(x = x_compressed, y = y_compressed)
}


simulate_data_parameters <- function(sigma12=par_est_mean[1], sigma22=par_est_mean[2], 
                                     h=par_est_mean[3], d=par_est_mean[4], e=par_est_mean[5], 
                                     f=par_est_mean[6], k=par_est_mean[7], p=par_est_mean[8], 
                                     m=par_est_mean[9], tau1=par_est_mean[10], 
                                     tau2=par_est_mean[11], b=par_est_mean[12], 
                                     c=par_est_mean[13], int=par_est_mean[14], # playing around with new parameter values
                                     tau3=par_est_mean[15], major=par_est_mean[16],
                                     minor=par_est_mean[17], xcent=par_est_mean[18],
                                     ycent=par_est_mean[19],
                                     beta0 = par_est_mean[20]){ 
  a = 2*pi/(tau1 - tau3) 
  #a is defined in terms of tau1
  
  max <- beta0 + max(int, h / sqrt(2 * pi * sigma12) + p / sqrt(2 * pi * sigma22))
  starting_lambda <- max * 540 * 440 * 75
  #intensity for poisson distribution
  starting <- rpois(n=1, lambda=starting_lambda)
  #simulating a point that has a poisson distribution with defined intensity
  
  x_vec <- runif(starting, min=-15, max=525)
  y_vec <- runif(starting, min = -15, max = 425)
  t_vec <- runif(starting, min=0, max=75)
  #generating points
  
  indicator_v <- (((x_vec - xcent)^2)/(major^2) + ((y_vec - ycent)^2)/(minor^2) <= 1 & t_vec <= tau3)
  #indicator of points in ventral patch
  
  ventral_time <- data.frame(t = t_vec[which(indicator_v == TRUE)],
                             y = y_vec[which(indicator_v == TRUE)],
                             x = x_vec[which(indicator_v == TRUE)])
  
  # ggplot(ventral_time, aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c()
  
  no_ventral_time <- data.frame(t = t_vec[which(indicator_v == FALSE)],
                                y = y_vec[which(indicator_v == FALSE)],
                                x = x_vec[which(indicator_v == FALSE)]) #set to tau3 parameter
  
  combined <- rbind(ventral_time, no_ventral_time)
  x_vec <- combined$x
  y_vec <- combined$y
  t_vec <- combined$t
  
  
  # ggplot(no_ventral_time, aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c()
  
  xoft <- d+b*cos(a*no_ventral_time$t + f)
  yoft <- e+c*sin(a*no_ventral_time$t + f)
  
  central_indicator <- ((ventral_time$x - xcent)^2)/(major^2) + ((ventral_time$y - ycent)^2)/(minor^2) <= 1
  #indicator of a point being in the ventral patch ellipse
  
  ventral <- ifelse(central_indicator, int, 0) 
  #for points within ventral patch ellipse
  
  ventral <- c(ventral, rep(0,length(no_ventral_time$t)))
  
  lambdas_rotational <- c(rep(0,length(ventral_time$t)),
                          (no_ventral_time$t <= tau1 & no_ventral_time$t>=tau3) * (h / (2 * pi * sigma12)) * exp( - ((no_ventral_time$x - xoft)^2 + (no_ventral_time$y - yoft)^2) / (2 * sigma12) ))
  #rotational intensity function component
  
  lambdas_radial <- c(rep(0,length(ventral_time$t)), (no_ventral_time$t <= tau2 & no_ventral_time$t>=tau3) * (p / sqrt(2 * pi * sigma22)) * exp( - (sqrt((no_ventral_time$x - d)^2 + (no_ventral_time$y - e)^2/m^2) - k*(no_ventral_time$t - tau3))^2 / (2 * sigma22)))
  #radial intensity function component
  
  lambda <- beta0 + ventral + lambdas_rotational + lambdas_radial
  #combined by adding components together
  
  prob <- if (max != 0) lambda / max else break
  
  retain_indices <- runif(starting)<=prob
  
  x <- x_vec[retain_indices] 
  y <- y_vec[retain_indices]
  t <- t_vec[retain_indices]
  
  data_compressed <- cbind(compress_points(x=x, y=y, d=d, e=e), t)
  
  return(data_compressed)
  
}

# simulating data using the function
sim_data <- data.frame(simulate_data_parameters())

ggplot(sim_data, aes(x=x,y=y,color=t)) + geom_point() + scale_color_viridis_c() + xlim(0, 550) + ylim(0, 450)

ggplot(sample_n(sweetrain,3000), aes(x=x,y=y,color=t)) + geom_point(size=2) +scale_color_viridis_c()
# ggplot(sample_n(data_expanded_final_new %>% filter(intensity==1),3000), aes(x=x,y=y,color=t)) + geom_point(size=2) +scale_color_viridis_c()



#------radial component only------#
#---beta0 also
simulate_data_parameters <- function(sigma12=par_est_mean[1], sigma22=par_est_mean[2], 
                                     h=par_est_mean[3], d=par_est_mean[4], e=par_est_mean[5], 
                                     f=par_est_mean[6], k=par_est_mean[7], p=par_est_mean[8], 
                                     m=par_est_mean[9], tau1=par_est_mean[10], 
                                     tau2=par_est_mean[11], b=par_est_mean[12], 
                                     c=par_est_mean[13], int=par_est_mean[14], 
                                     tau3=par_est_mean[15], major=par_est_mean[16],
                                     minor=par_est_mean[17], xcent=par_est_mean[18],
                                     ycent=par_est_mean[19],
                                     beta0 = par_est_mean[20]){ 
  a = 2*pi/(tau1 - tau3) 
  #a is defined in terms of tau1
  
  max <- beta0 + max(int, h / sqrt(2 * pi * sigma12) + p / sqrt(2 * pi * sigma22))
  starting_lambda <- max * 540 * 440 * 75
  #intensity for poisson distribution
  starting <- rpois(n=1, lambda=starting_lambda)
  #simulating a point that has a poisson distribution with defined intensity
  
  x_vec <- runif(starting, min=-15, max=525)
  y_vec <- runif(starting, min = -15, max = 425)
  t_vec <- runif(starting, min=0, max=75)
  #generating points
  
  indicator_v <- (((x_vec - xcent)^2)/(major^2) + ((y_vec - ycent)^2)/(minor^2) <= 1 & t_vec <= tau3)
  #indicator of points in ventral patch
  
  ventral_time <- data.frame(t = t_vec[which(indicator_v == TRUE)],
                             y = y_vec[which(indicator_v == TRUE)],
                             x = x_vec[which(indicator_v == TRUE)])
  
  # ggplot(ventral_time, aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c()
  
  no_ventral_time <- data.frame(t = t_vec[which(indicator_v == FALSE)],
                                y = y_vec[which(indicator_v == FALSE)],
                                x = x_vec[which(indicator_v == FALSE)]) #set to tau3 parameter
  
  combined <- rbind(ventral_time, no_ventral_time)
  x_vec <- combined$x
  y_vec <- combined$y
  t_vec <- combined$t
  
  
  # ggplot(no_ventral_time, aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c()
  
  xoft <- d+b*cos(a*no_ventral_time$t + f)
  yoft <- e+c*sin(a*no_ventral_time$t + f)
  
  central_indicator <- ((ventral_time$x - xcent)^2)/(major^2) + ((ventral_time$y - ycent)^2)/(minor^2) <= 1
  #indicator of a point being in the ventral patch ellipse
  
  ventral <- ifelse(central_indicator, int, 0) 
  #for points within ventral patch ellipse
  
  ventral <- c(ventral, rep(0,length(no_ventral_time$t)))
  
  lambdas_rotational <- c(rep(0,length(ventral_time$t)),
                          (no_ventral_time$t <= tau1 & no_ventral_time$t>=tau3) * (h / (2 * pi * sigma12)) * exp( - ((no_ventral_time$x - xoft)^2 + (no_ventral_time$y - yoft)^2) / (2 * sigma12) ))
  #rotational intensity function component
  
  lambdas_radial <- c(rep(0,length(ventral_time$t)), (no_ventral_time$t <= tau2 & no_ventral_time$t>=tau3) * (p / sqrt(2 * pi * sigma22)) * exp( - (sqrt((no_ventral_time$x - d)^2 + (no_ventral_time$y - e)^2/m^2) - k*(no_ventral_time$t - tau3))^2 / (2 * sigma22)))
  #radial intensity function component
  
  lambda <- beta0 + lambdas_radial
  #combined by adding components together
  
  prob <- if (max != 0) lambda / max else break
  
  retain_indices <- runif(starting)<=prob
  
  x <- x_vec[retain_indices] 
  y <- y_vec[retain_indices]
  t <- t_vec[retain_indices]
  
  data <- cbind(x, y, t)
  return(data)
}

# simulating data using the function
sim_data <- data.frame(simulate_data_parameters())

ggplot(sim_data, aes(x=x,y=y,color=t)) + geom_point() + scale_color_viridis_c() + xlim(0, 550) + ylim(0, 450) 




#------rotational component only------#
#---beta0 also
simulate_data_parameters <- function(sigma12=par_est_mean[1], sigma22=par_est_mean[2], 
                                     h=par_est_mean[3], d=par_est_mean[4], e=par_est_mean[5], 
                                     f=par_est_mean[6], k=par_est_mean[7], p=par_est_mean[8], 
                                     m=par_est_mean[9], tau1=par_est_mean[10], 
                                     tau2=par_est_mean[11], b=par_est_mean[12], 
                                     c=par_est_mean[13], int=par_est_mean[14],
                                     tau3=par_est_mean[15], major=par_est_mean[16],
                                     minor=par_est_mean[17], xcent=par_est_mean[18],
                                     ycent=par_est_mean[19],
                                     beta0 = par_est_mean[20]){ 
  a = 2*pi/(tau1 - tau3) 
  #a is defined in terms of tau1
  
  max <- beta0 + max(int, h / sqrt(2 * pi * sigma12) + p / sqrt(2 * pi * sigma22))
  starting_lambda <- max * 540 * 440 * 75
  #intensity for poisson distribution
  starting <- rpois(n=1, lambda=starting_lambda)
  #simulating a point that has a poisson distribution with defined intensity
  
  x_vec <- runif(starting, min=-15, max=525)
  y_vec <- runif(starting, min = -15, max = 425)
  t_vec <- runif(starting, min=0, max=75)
  #generating points
  
  indicator_v <- (((x_vec - xcent)^2)/(major^2) + ((y_vec - ycent)^2)/(minor^2) <= 1 & t_vec <= tau3)
  #indicator of points in ventral patch
  
  ventral_time <- data.frame(t = t_vec[which(indicator_v == TRUE)],
                             y = y_vec[which(indicator_v == TRUE)],
                             x = x_vec[which(indicator_v == TRUE)])
  
  # ggplot(ventral_time, aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c()
  
  no_ventral_time <- data.frame(t = t_vec[which(indicator_v == FALSE)],
                                y = y_vec[which(indicator_v == FALSE)],
                                x = x_vec[which(indicator_v == FALSE)]) #set to tau3 parameter
  
  combined <- rbind(ventral_time, no_ventral_time)
  x_vec <- combined$x
  y_vec <- combined$y
  t_vec <- combined$t
  
  
  # ggplot(no_ventral_time, aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c()
  
  xoft <- d+b*cos(a*no_ventral_time$t + f)
  yoft <- e+c*sin(a*no_ventral_time$t + f)
  
  central_indicator <- ((ventral_time$x - xcent)^2)/(major^2) + ((ventral_time$y - ycent)^2)/(minor^2) <= 1
  #indicator of a point being in the ventral patch ellipse
  
  ventral <- ifelse(central_indicator, int, 0) 
  #for points within ventral patch ellipse
  
  ventral <- c(ventral, rep(0,length(no_ventral_time$t)))
  
  lambdas_rotational <- c(rep(0,length(ventral_time$t)),
                          (no_ventral_time$t <= tau1 & no_ventral_time$t>=tau3) * (h / (2 * pi * sigma12)) * exp( - ((no_ventral_time$x - xoft)^2 + (no_ventral_time$y - yoft)^2) / (2 * sigma12) ))
  #rotational intensity function component
  
  lambdas_radial <- c(rep(0,length(ventral_time$t)), (no_ventral_time$t <= tau2 & no_ventral_time$t>=tau3) * (p / sqrt(2 * pi * sigma22)) * exp( - (sqrt((no_ventral_time$x - d)^2 + (no_ventral_time$y - e)^2/m^2) - k*(no_ventral_time$t - tau3))^2 / (2 * sigma22)))
  #radial intensity function component
  
  lambda <- beta0 + lambdas_rotational
  #combined by adding components together
  
  prob <- if (max != 0) lambda / max else break
  
  retain_indices <- runif(starting)<=prob
  
  x <- x_vec[retain_indices] 
  y <- y_vec[retain_indices]
  t <- t_vec[retain_indices]
  
  data <- cbind(x, y, t)
  return(data)
}

# simulating data using the function
sim_data <- data.frame(simulate_data_parameters())

ggplot(sim_data, aes(x=x,y=y,color=t)) + geom_point() + scale_color_viridis_c() + xlim(0, 550) + ylim(0, 450) 

